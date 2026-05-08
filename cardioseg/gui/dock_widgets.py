import re
import ast
from qtpy.QtWidgets import (
    QWidget, QVBoxLayout, QPushButton, QLabel,
    QProgressBar, QFileDialog, QGroupBox, QSpinBox,
    QLineEdit, QFrame, QHBoxLayout, QScrollArea,
    QDialog, QListWidget, QListWidgetItem, QAbstractItemView,
    QDialogButtonBox, QInputDialog, QCheckBox, QDoubleSpinBox, QComboBox,
    QMessageBox
)

from qtpy.QtCore import Qt, QTimer
from qtpy.QtCore import QThread, Signal

from qtpy.QtGui import QColor, QPalette

from napari.qt import QtViewer
from napari.utils.notifications import show_info
from napari.utils.colormaps import DirectLabelColormap

from PIL import Image

from czifile import imread as czi_imread

from collections import Counter

from cardioseg.segmentation.cellpose_wrapper import CellposeSegmenter
from cardioseg.io.exporters import export_cell_measurements
from cardioseg.spatial.engine import diagnostic_polygon_check

import tifffile
import time
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm 
import matplotlib.colors as mcolors
import scanpy as sc
import imageio.v3 as iio
import os
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed

class CellData:
    def __init__(self):
        self.measurements = {}
        self.genes = {}
    
    def get_gene_expression(self, cid, gene):
        return self.genes.get(cid, {}).get(gene, 0.0)

    def get_top_gene(self, cid):
        genes = self.genes.get(cid, {})
        if not genes:
            return "Unknown"
        return max(genes, key=genes.get)
class CellposeWorker(QThread):
    progress = Signal(float)
    finished = Signal(dict)

    def __init__(self, image, diameter):
        super().__init__()
        self.image = image
        self.diameter = diameter
        self.segmenter = CellposeSegmenter()
    
    def run(self):
        from time import perf_counter  

        t0 = perf_counter()

        def progress_callback(p):
            self.progress.emit(p)
       
        result = self.segmenter.segment(
            self.image,
            diameter=self.diameter,
            progress_callback=progress_callback
    
        )

        self.finished.emit({
            "masks": result["masks"],
            "elapsed": perf_counter() - t0
        })

class SpatialAnnotationWorker(QThread):

    progress = Signal(float)
    finished = Signal(dict)

    def __init__(self, nuclei_labels, label_ids, engine, annotate_cell_type = True, annotate_genes = True):
        super().__init__()
        self.nuclei_labels = nuclei_labels
        self.label_ids = label_ids
        self.engine = engine
        self.annotate_cell_type = annotate_cell_type
        self.annotate_genes = annotate_genes
    
    def run(self):
        results = {}
        # diagnostic_file = "debug_polygon_diagnostics.txt"

        # if os.path.exists(diagnostic_file):
        #     os.remove(diagnostic_file)

        total = len(self.label_ids)
        offset_y, offset_x = getattr(self, "roi_offset", (0, 0))

        if self._is_windows() and total >= 100:
            max_workers = min(8, max(1, (os.cpu_count() or 1) - 1))
            with ThreadPoolExecutor(max_workers=max_workers) as pool:
                future_map = {
                    pool.submit(self._annotate_single, cid, offset_y, offset_x): cid
                    for cid in self.label_ids
                }
                completed = 0
                for future in as_completed(future_map):
                    cid, payload = future.result()
                    if payload is not None:
                        results[cid] = payload
                    completed += 1
                    self.progress.emit(completed / total)
        else:
            for i, cid in enumerate(self.label_ids):
                _, payload = self._annotate_single(cid, offset_y, offset_x)
                if payload is not None:
                    results[cid] = payload
                self.progress.emit((i + 1) / total)
       
        self.finished.emit(results)

    @staticmethod
    def _is_windows():
        return sys.platform.startswith("win")

    def _annotate_single(self, cid, offset_y, offset_x):
        from skimage.measure import find_contours

        mask = self.nuclei_labels == cid
        contours = find_contours(mask.astype(np.uint8), 0.5)
        if not contours:
                return cid, None

        poly = contours[0]
        poly[:, 0] += offset_y
        poly[:, 1] += offset_x

        poly_final = np.fliplr(poly)
        centroid = np.mean(poly_final, axis=0)
        poly_buffered = centroid + 1.2 * (poly_final - centroid)

        annotation_result = self.engine.annotate_polygon(
            poly_buffered,
            centroid,
            annotate_cell_type=self.annotate_cell_type,
            annotate_genes=self.annotate_genes
        )
        return cid, {
            "cell_type": annotation_result["cell_type"],
            "confidence": annotation_result["confidence"],
            "score": annotation_result["score"],
            "genes": annotation_result.get("genes", {}),
             "genes_raw": annotation_result.get("genes_raw", {}),
        }

class ColorManager: 

    def __init__(self, cell_measurements, cell_type_colors):
        self.measurements = cell_measurements
        self.cell_type_colors = cell_type_colors

    def get_mode(self, cell_ids, key):

        if key.startswith("gene:"):
            return "continuous"

        if key == "top_gene":
            return "categorical"

        val = self._get_sample(cell_ids, key)

        if isinstance(val, str):
            return "categorical"
        elif isinstance(val, (int, float)):
            return "continuous"

        return "unknown"
    
    def _get_sample(self, cell_ids, key):
        """Get a representative value to determine type (str vs float)."""

        for cid in cell_ids:
            meas = self.measurements.get(cid, {})

            if key.startswith("gene:"):
                gene = key.split(":", 1)[1]
                val = self._get_gene_expression(meas, gene)
            else:
                val = meas.get(key, None)

            if val is not None:
                return val

        return None
    
    def categorical(self, cell_ids, key, labels):
        filtered = np.zeros_like(labels, dtype=np.int32)

        cat_to_label = {}
        label_to_color = {0: (0, 0, 0, 0)}
        next_label = 1

        for cid in cell_ids:
            meas = self.measurements.get(cid, {})
            val = meas.get(key, "Unknown")
            if key == "cell_type":
                canonical_val = self._canonical_cell_type(val)
                color = self.cell_type_colors.get(canonical_val, self.cell_type_colors["Unknown"])

            if val not in cat_to_label:
                cat_to_label[val] = next_label

                # assign color ONCE per category
                if key == "cell_type":
                    color = self.cell_type_colors.get(val, self.cell_type_colors["Unknown"])
                else:
                    np.random.seed(hash(val) % 2**32)
                    color = (*np.random.rand(3), 1)

                label_to_color[next_label] = color
                next_label += 1

            filtered[labels == cid] = cat_to_label[val]

        legend = {
            "type": "categorical",
            "mapping": cat_to_label,
            "colors": label_to_color,
            "field": key
        }

        return filtered, label_to_color, legend
    
    def _canonical_cell_type(self, value):
        if value in self.cell_type_colors:
            return value

        if isinstance(value, str):
            v = value.strip()
            if v in self.cell_type_colors:
                return v

            normalized = v.lower()
            for known in self.cell_type_colors:
                if known.lower() == normalized:
                    return known
        return "Unknown"

    def _get_gene_expression(self, meas, gene_name):
        genes = meas.get("genes", {})
        if not isinstance(genes, dict):
            return None
        if gene_name in genes:
            return genes.get(gene_name)
        if isinstance(gene_name, str):
            target = gene_name.strip().lower()
            for g, v in genes.items():
                if isinstance(g, str) and g.lower() == target:
                    return v
        return None

    def _resolve_cell_type_color(self, value):
        if value in self.cell_type_colors:
            return self.cell_type_colors[value]

        if isinstance(value, str):
            v = value.strip()
            if v in self.cell_type_colors:
                return self.cell_type_colors[v]

            normalized = v.lower()
            for known, color in self.cell_type_colors.items():
                if known.lower() == normalized:
                    return color

        np.random.seed(hash(value) % 2**32)
        return (*np.random.rand(3), 1)

    def continuous(self, cell_ids, key, labels):
        filtered = np.zeros_like(labels, dtype=np.int32)

        id_map = {}
        values = {}

        # --- STEP 1: remap IDs → 1..N ---
        for i, cid in enumerate(cell_ids, start=1):
            meas = self.measurements.get(cid, {})

            if key.startswith("gene:"):
                gene = key.split(":", 1)[1]
                val = self._get_gene_expression(meas, gene)
                if val is None:
                    val = 0.0
            else:
                val = meas.get(key, 0.0)

            id_map[cid] = i
            values[i] = float(val)
            filtered[labels == cid] = i

        # --- STEP 2: compute limits ---
        if len(values) > 0:
            all_vals = np.array(list(values.values()), dtype=float)
            if key.startswith("gene:"):
                # Keep gene expression scales fully continuous for the selected view.
                vmin, vmax = float(np.min(all_vals)), float(np.max(all_vals))
            else:
                # Robust scaling for non-gene quantitative metadata.
                vmin, vmax = np.percentile(all_vals, [5, 85])
                vmin, vmax = float(vmin), float(vmax)
        else:
            vmin, vmax = 0.0, 1.0

         # Avoid a flat contrast interval when all selected cells share one value.
        if np.isclose(vmin, vmax):
            vmax = vmin + 1e-6

        legend = {
            "type": "continuous",
            "min": float(vmin),
            "max": float(vmax),
            "colormap": "fire"
        }

        return filtered, values, legend

    def top_gene(self, cell_ids, labels):
        mapping = {}
        filtered = np.zeros_like(labels, dtype=np.int32)
        next_label = 1

        for cid in cell_ids:
            gene = self.measurements.get(cid, {}).get("top_gene", "Unknown")
            mapping[cid] = gene
            filtered[labels == cid] = cid

        color_dict = {0: (0,0,0,0)}
        for cid, gene in mapping.items():
            np.random.seed(hash(gene) % 2**32)
            color_dict[cid] = (*np.random.rand(3), 1)

        legend = {
            "type": "categorical",
            "mapping": mapping,
            "colors": color_dict
        }

        return filtered, color_dict, legend

class ControlPanel(QWidget):
    def __init__(self, viewer, info_panel, table_panel):
        super().__init__()
        self.viewer = viewer
        self.info_panel = info_panel
        self.table_panel = table_panel

        self.cell_labels = None
        self.nuclei_labels = None

        self.cell_type_colors = {
        "Cardiomyocyte": (0.84, 0.15, 0.16, 1.0),          # red
        "Fibroblast": (0.12, 0.47, 0.71, 1.0),             # blue
        "Endothelial": (0.17, 0.63, 0.17, 1.0),            # green
        "Lymphatic endothelial": (0.0, 0.75, 0.75, 1.0),   # cyan
        "Mast": (1.0, 0.6, 0.2, 1.0),                      # orange
        "Mesothelial": (0.6, 0.6, 0.6, 1.0),               # gray
        "Myeloid": (0.58, 0.4, 0.74, 1.0),                 # purple
        "Neuronal": (0.55, 0.34, 0.29, 1.0),               # brown
        "NK_T": (0.89, 0.47, 0.76, 1.0),                   # pink
        "Pericyte": (0.5, 0.5, 0.0, 1.0),                  # olive
        "Smooth muscle": (1.0, 0.85, 0.0, 1.0),            # yellow
        "Unknown": (0.7, 0.7, 0.7, 1.0)
    }

        self.cell_measurements = {}
        self.available_genes = set()

        self.viewer.mouse_drag_callbacks.append(self._on_mouse_click)

        self.add_search_bar()
        self._selected_label = None

        self.spatial_engine = None
        self._spatial_initialized = False

        self.nucleus_annotations = {}  
        self.roi_annotations = {}
        self._csv_metadata_to_restore = None
        self._csv_restore_notice = ""

        self.setLayout(self._build_layout())

        self.measurements = self.cell_measurements
        self.color_manager = ColorManager(
            cell_measurements=self.cell_measurements,
            cell_type_colors=self.cell_type_colors
        )

        self.viewer.layers.events.removed.connect(self._on_layer_removed)

    def _on_layer_removed(self, event):
        layer = getattr(event, "value", None)
        if layer is None:
            return

        if layer.name == "Nuclei (Cellpose)":
            self.nuclei_labels = None
            self.all_label_ids = np.array([], dtype=np.uint32)
            self.cell_measurements = {}
            self.color_manager.measurements = self.cell_measurements
            self.table_panel.set_measurements({})
            self.info_panel.clear()
            if hasattr(self, "search_status"):
                self.search_status.setText("")
            self.update_legend({"type": "none"})

    def _build_layout(self):
        layout = QVBoxLayout()

        self.cellpose_progress = QProgressBar()
        self.cellpose_progress.setRange(0, 100)
        self.cellpose_progress.setVisible(False)

        self.cellpose_timer_label = QLabel("Elapsed: 0.0 s")
        self.cellpose_timer_label.setAlignment(Qt.AlignCenter)
        self.cellpose_timer_label.setVisible(False)

        layout.addWidget(self.cellpose_progress)
        layout.addWidget(self.cellpose_timer_label)

        title = QLabel("CardioSeg")
        title.setAlignment(Qt.AlignCenter)
        layout.addWidget(title)

        layout.addWidget(self._image_group())
        layout.addWidget(self._segmentation_group())
        layout.addWidget(self._editing_group())

        layout.addStretch()
        return layout

    def _image_group(self):
        group = QGroupBox("Image")
        layout = QVBoxLayout()
        btn = QPushButton("Load Image")
        btn.clicked.connect(self.load_image)
        layout.addWidget(btn)
        group.setLayout(layout)
        return group

    def _segmentation_group(self):
        group = QGroupBox("Segmentation")
        layout = QVBoxLayout()

        self.diameter_spin = QSpinBox()
        self.diameter_spin.setRange(5, 200)
        self.diameter_spin.setValue(30)
        self.diameter_spin.setPrefix("Diameter: ")

        run_nuclei = QPushButton("Run Cellpose (Global)")
        run_nuclei.clicked.connect(self.run_cellpose_nuclei)

        run_roi = QPushButton("Segment ROI Only")
        run_roi.clicked.connect(self.run_cellpose_roi)

        # export_binary_mask = QPushButton("Export Nuclei Binary Mask (PNG)")
        # export_binary_mask.clicked.connect(self.export_nuclei_binary_mask_png)

        layout.addWidget(self.diameter_spin)
        layout.addWidget(run_nuclei)
        layout.addWidget(run_roi)
        # layout.addWidget(export_binary_mask)
        group.setLayout(layout)
        return group

    def _editing_group(self):
        group = QGroupBox("Editing")
        layout = QVBoxLayout()

        annotation_btn = QPushButton("Annotate Nuclei")
        annotation_btn.clicked.connect(self.annotate_nuclei)
        layout.addWidget(annotation_btn)

        save_btn = QPushButton("Save Measurements (CSV)")
        save_btn.clicked.connect(self.save_measurements)
        layout.addWidget(save_btn)

        load_btn = QPushButton("Load Measurements (CSV)")
        load_btn.clicked.connect(self.load_measurements_csv)
        layout.addWidget(load_btn)

        # roi_export_btn = QPushButton("Export ROI Measurements (CSV)")
        # roi_export_btn.clicked.connect(self.export_roi_csv)
        # layout.addWidget(roi_export_btn)

        # export_btn = QPushButton("Export Layer(s) as PNG")
        # export_btn.clicked.connect(self.export_png)
        # layout.addWidget(export_btn)

        # ruler_btn = QPushButton("Ruler / Measure")
        # ruler_btn.clicked.connect(self.add_ruler)
        # layout.addWidget(ruler_btn)

        group.setLayout(layout)
        return group
    
    def export_nuclei_binary_mask_png(self):
        if self.nuclei_labels is None:
            show_info("No nuclei segmentation found. Run segmentation first.")
            return

        path, _ = QFileDialog.getSaveFileName(
            self,
            "Save Nuclei Binary Mask",
            "nuclei_binary_mask.png",
            "PNG (*.png)"
        )
        if not path:
            return

        binary_mask = (self.nuclei_labels > 0).astype(np.uint8) * 255

        try:
            iio.imwrite(path, binary_mask)
        except Exception as exc:
            show_info(f"Failed to save binary mask PNG: {exc}")
            return

        show_info(f"Binary nuclei mask saved to {path}")
 
    def load_measurements_csv(self):
        import pandas as pd

        path, _ = QFileDialog.getOpenFileName(
            self,
            "Load Measurements CSV",
            "",
            "CSV (*.csv)"
        )
        if not path:
            return

        try:
            df = pd.read_csv(path)
        except Exception as exc:
            show_info(f"Failed to load CSV: {exc}")
            return

        if df.empty:
            show_info("CSV is empty.")
            return

        if "cell_id" in df.columns:
            cell_ids = df["cell_id"]
        elif "Unnamed: 0" in df.columns:
            cell_ids = df["Unnamed: 0"]
        else:
            show_info("CSV must contain a 'cell_id' column (or exported index column).")
            return

        measurements = {}
        for idx, row in df.iterrows():
            try:
                cid = int(cell_ids.iloc[idx])
            except (TypeError, ValueError):
                continue

            meas = {}
            genes = {}
            genes_raw = {}

            for col, value in row.items():
                if col in ("cell_id", "Unnamed: 0"):
                    continue
                if isinstance(value, float) and np.isnan(value):
                    continue

                if col.startswith("gene_norm_"):
                    genes[col.replace("gene_norm_", "", 1)] = float(value)
                    continue

                if col.startswith("gene_raw_"):
                    genes_raw[col.replace("gene_raw_", "", 1)] = float(value)
                    continue

                meas[col] = self._coerce_csv_value(value)

            if genes:
                meas["genes"] = genes
                if "top_gene" not in meas or not meas.get("top_gene"):
                    if genes:  # Safety check to ensure genes dict is not empty
                        top_gene = max(genes, key=genes.get)
                        meas["top_gene"] = top_gene
                        meas["top_gene_expr"] = float(genes[top_gene])
                    else:
                        meas["top_gene"] = "Unknown"
                        meas["top_gene_expr"] = 0.0
            if genes_raw:
                meas["genes_raw"] = genes_raw

            measurements[cid] = meas

        if not measurements:
            show_info("No valid cell measurements found in CSV.")
            return

        self.cell_measurements = measurements
        self.measurements = self.cell_measurements
        self.color_manager.measurements = self.cell_measurements
        self.table_panel.set_measurements(self.cell_measurements)
        self.available_genes = set()
        self._refresh_available_genes()

        has_perimeter_coords = any(
            isinstance(meas.get("Perimeter coordinates"), (list, tuple)) and len(meas.get("Perimeter coordinates")) > 0
            for meas in measurements.values()
        )
        choices = self._prompt_csv_import_mode(has_perimeter_coords)
        if choices is None:
            return

        labels, valid_centroids = self._build_labels_from_measurements(
            measurements,
            use_perimeter=choices["use_perimeter"],
            use_circles=choices["use_circles"],
        )

        if labels is not None:
            if "Nuclei (Cellpose)" in self.viewer.layers:
                self.viewer.layers.remove("Nuclei (Cellpose)")
            self.viewer.add_labels(labels, name="Nuclei (Cellpose)")
            self.nuclei_labels = labels
            self.all_label_ids = np.unique(labels)
            self.all_label_ids = self.all_label_ids[self.all_label_ids != 0]

            if choices["use_cellpose"]:
                self._csv_metadata_to_restore = {cid: dict(meas) for cid, meas in measurements.items()}
                self._run_cellpose_in_csv_extent(measurements)
                show_info(
                    f"Loaded {len(measurements)} nuclei from CSV. Seeded {valid_centroids} nuclei and running Cellpose in CSV area."
                )
            else:
                mode_name = "perimeter coordinates" if choices["use_perimeter"] else "centroid circles"
                show_info(
                    f"Loaded {len(measurements)} nuclei from CSV. Reconstructed masks for {valid_centroids} nuclei using {mode_name}."
                )
        else:
            show_info(f"Loaded {len(measurements)} nuclei from CSV (no centroid/image available to reconstruct masks).")

    def _coerce_csv_value(self, value):
        if isinstance(value, str):
            stripped = value.strip()
            if stripped.startswith(("(", "[")) and stripped.endswith((")", "]")):
                try:
                    return ast.literal_eval(stripped)
                except (ValueError, SyntaxError):
                    return value
            return value
        return value

    def _prompt_csv_import_mode(self, has_perimeter_coords):
        message = QMessageBox(self)
        message.setWindowTitle("CSV import mode")
        message.setIcon(QMessageBox.Question)
        message.setText("Choose how to initialize nuclei from this CSV.")

        if has_perimeter_coords:
            perimeter_btn = message.addButton("Use perimeter reconstruction", QMessageBox.AcceptRole)
            circles_btn = None
        else:
            perimeter_btn = None
            circles_btn = message.addButton("Use fast centroid circles", QMessageBox.AcceptRole)
        cellpose_btn = message.addButton("Use Cellpose in CSV area", QMessageBox.ActionRole)
        cancel_btn = message.addButton(QMessageBox.Cancel)

        message.exec()
        clicked = message.clickedButton()
        if clicked == cancel_btn:
            return None
        if has_perimeter_coords and clicked == perimeter_btn:
            return {"use_perimeter": True, "use_circles": False, "use_cellpose": False}
        if (not has_perimeter_coords) and clicked == circles_btn:
            return {"use_perimeter": False, "use_circles": True, "use_cellpose": False}
        return {"use_perimeter": has_perimeter_coords, "use_circles": not has_perimeter_coords, "use_cellpose": True}
    
    def _build_labels_from_measurements(self, measurements, use_perimeter=False, use_circles=True):
        if "Image" not in self.viewer.layers:
            return None, 0
        image_layer = self.viewer.layers["Image"]

        image_shape = image_layer.data.shape[:2]
        labels = np.zeros(image_shape, dtype=np.uint32)
        valid_centroids = 0

        for cid, meas in measurements.items():
            if use_perimeter:
                contours = meas.get("Perimeter coordinates")
                if isinstance(contours, (list, tuple)):
                    from skimage.draw import polygon
                    for contour in contours:
                        contour_np = np.asarray(contour, dtype=float)
                        if contour_np.ndim != 2 or contour_np.shape[1] < 2 or contour_np.shape[0] < 3:
                            continue
                        contour_np = np.round(contour_np[:, :2]).astype(int)
                        rr, cc = polygon(contour_np[:, 1], contour_np[:, 0], shape=image_shape)
                        if rr.size == 0:
                            continue
                        labels[rr, cc] = cid
                        valid_centroids += 1
                        break
                continue

            if not use_circles:
                continue

            centroid = meas.get("centroid")
            if not isinstance(centroid, (list, tuple)) or len(centroid) < 2:
                continue

            try:
                y = int(round(float(centroid[0])))
                x = int(round(float(centroid[1])))
            except (TypeError, ValueError):
                continue

            if y < 0 or x < 0 or y >= image_shape[0] or x >= image_shape[1]:
                continue

            area = meas.get("area", 1)
            try:
                radius = max(1, int(round(np.sqrt(float(area) / np.pi))))
            except (TypeError, ValueError):
                radius = 1

            y_min = max(0, y - radius)
            y_max = min(image_shape[0], y + radius + 1)
            x_min = max(0, x - radius)
            x_max = min(image_shape[1], x + radius + 1)

            yy, xx = np.ogrid[y_min:y_max, x_min:x_max]
            circle_mask = (yy - y) ** 2 + (xx - x) ** 2 <= radius ** 2

            target = labels[y_min:y_max, x_min:x_max]
            target[(circle_mask) & (target == 0)] = cid
            valid_centroids += 1

        return labels, valid_centroids
    
    def _run_cellpose_in_csv_extent(self, measurements):
        if "Image" not in self.viewer.layers:
            return

        image = self.viewer.layers["Image"].data
        image_shape = image.shape[:2]
        ys, xs = [], []
        for meas in measurements.values():
            centroid = meas.get("centroid")
            if isinstance(centroid, (list, tuple)) and len(centroid) >= 2:
                try:
                    ys.append(int(round(float(centroid[0]))))
                    xs.append(int(round(float(centroid[1]))))
                except (TypeError, ValueError):
                    pass

        if not ys or not xs:
            return

        pad = max(10, int(self.diameter_spin.value()))
        min_y, max_y = max(0, min(ys) - pad), min(image_shape[0], max(ys) + pad)
        min_x, max_x = max(0, min(xs) - pad), min(image_shape[1], max(xs) + pad)
        if min_y >= max_y or min_x >= max_x:
            return

        self._roi_offset = (min_y, min_x)
        self._current_image = image
        roi_crop = image[min_y:max_y, min_x:max_x]

        self.worker = CellposeWorker(roi_crop, self.diameter_spin.value())
        self.worker.progress.connect(self.update_progress)
        self.worker.finished.connect(self.finish_roi_segmentation)
        self.worker.start()

    def _refresh_available_genes(self):
        self.available_genes = set()
        self.gene_min = {}
        self.gene_max = {}

        for cell in self.cell_measurements.values():
            if isinstance(cell.get("genes"), dict):
                self.available_genes.update(cell["genes"].keys())

        for gene in self.available_genes:
            values = [
                cell.get("genes", {}).get(gene, 0.0)
                for cell in self.cell_measurements.values()
            ]
            if not values:
                continue
            self.gene_min[gene] = float(np.min(values))
            self.gene_max[gene] = float(np.max(values))

    def _initialize_spatial_engine(self):

        path, _ = QFileDialog.getOpenFileName(
            self,
            "Select filtered_feature_bc_matrix.h5",
            "",
            "H5 (*.h5)"
        )

        if not path:
            return False

        from cardioseg.spatial.engine import SpatialExpressionEngine

        self.spatial_engine = SpatialExpressionEngine(path)

        self._spatial_initialized = True
        show_info("Spatial dataset loaded.")

        return True
      
    def add_search_bar(self):
        self.search_container = QWidget()
        self.search_container.setObjectName("search_container")

        layout = QHBoxLayout(self.search_container)
        layout.setContentsMargins(8, 2, 8, 2)
        layout.setSpacing(5)

        #magnifying glass
        self.icon_label = QLabel("🔍") # Unicode magnifying glass
        self.icon_label.setObjectName("search_icon")

        #input field
        self.search_input = QLineEdit()
        self.search_input.setPlaceholderText("Search (e.g. area > 500 or 1,2,5)")
        self.search_input.returnPressed.connect(self.execute_search)

        #status label
        self.search_status = QLabel("") 
        self.search_status.setObjectName("search_status")

        #reset button
        self.reset_button = QPushButton("✕")
        self.reset_button.setFixedSize(24, 24)
        self.reset_button.setCursor(Qt.PointingHandCursor)
        self.reset_button.setStyleSheet("color: #777; border: none; font-weight: bold;")
        self.reset_button.clicked.connect(self.clear_search)

        #adding widgets 
        layout.addWidget(self.icon_label)
        layout.addWidget(self.search_input)
        layout.addWidget(self.search_status)
        layout.addWidget(self.reset_button)

        # Use a modern, clean font stack
        font_family = "'Segoe UI', 'Helvetica Neue', sans-serif"

        self.search_container.setStyleSheet(f"""
            /* The Main Container */
            QWidget#search_container {{
                background-color: rgba(40, 44, 52, 220);
                border: 1px solid rgba(255, 255, 255, 30);
                border-radius: 12px;
            }}
            
            QWidget#search_container:hover {{
                border: 1px solid rgba(100, 150, 255, 100);
                background-color: rgba(45, 50, 60, 240);
            }}

            /* Magnifying Glass Icon (The Label on the left) */
            QLabel#search_icon {{
                color: #888;
                font-size: 14px;
                padding-left: 10px;
                background: transparent;
            }}

            /* The Input Field */
            QLineEdit {{
                background: transparent;
                border: none;
                color: #D1D1D1; /* Soft off-white */
                padding: 8px 5px;
                font-family: {font_family};
                font-size: 13px;
            }}
            
            QLineEdit:focus {{
                color: #FFFFFF;
            }}

            /* The Search Status (Removing the 'Gray Box' look) */
            QLabel#search_status {{
                color: #666; /* Muted text */
                font-family: {font_family};
                font-size: 11px;
                font-weight: 500;
                background: transparent; /* Removes the box look */
                padding-right: 12px;
                letter-spacing: 0.5px;
            }}

            /* The Buttons (Close/Clear) */
            QPushButton {{
                background: transparent;
                color: #777;
                border: none;
                padding-right: 10px;
                font-size: 16px;
            }}

            QPushButton:hover {{
                color: #ff5f57; /* Subtle red for 'close' or blue for others */
            }}
        """)

        self.viewer.window.add_dock_widget(
        self.search_container,
        area = "top",
        name = "Search",
        )

    def execute_search(self):
        query = self.search_input.text().strip()

        if not query:
            self.display_cells()
            return

        if query.lower() == "help":
            pass

        if query.lower() == "ruler":
            self.add_ruler()
            return

        #stats for all parameters
        if query.lower().startswith("stats(") and query.endswith(")"):
            param = query[6:-1].strip()
            self.run_stats(param)
            return
        
        # -------------------------
        # SORT
        # -------------------------
        if query.lower().startswith("sort:"):
            field = query.split(":", 1)[1].strip()

            if field not in self.cell_measurements[next(iter(self.cell_measurements))]:
                self.search_status.setText(f"Unknown field: {field}")
                return

            self.display_cells(cell_ids=None, color_by=field)
            self.search_status.setText(f"Colored all cells by {field}")
            return

        # Guard clauses
        if "Nuclei (Cellpose)" not in self.viewer.layers or not self.cell_measurements:
            self.search_status.setText("Run segmentation first")
            return

        if not query:
            self.clear_search()
            return
       
        try:
            # Split by OR first
            or_groups = re.split(r'\s+OR\s+', query, flags=re.IGNORECASE)

            final_matches = set()

            for group in or_groups:
                # Split each group by AND
                and_conditions = re.split(r'\s+AND\s+', group, flags=re.IGNORECASE)

                group_result = None

                for cond in and_conditions:
                    cond = cond.strip()

                    # Handle NOT
                    is_not = False
                    if cond.upper().startswith("NOT "):
                        is_not = True
                        cond = cond[4:].strip()

                    res = self._evaluate_condition(cond)

                    if is_not:
                        res = set(self.cell_measurements.keys()) - res

                    if group_result is None:
                        group_result = res
                    else:
                        group_result = group_result.intersection(res)

                # OR = union
                final_matches = final_matches.union(group_result)

            matched_cell_ids = final_matches

        except Exception as e:
            self.search_status.setText(f"Error: {str(e)}")
            self.clear_search()
            return

        count = len(matched_cell_ids)
        total = len(self.all_label_ids) if hasattr(self, "all_label_ids") else 0
        percent = (count / total * 100) if total else 0

        self.search_status.setText(f"{count} / {total} cells ({percent:.1f}%)")

        if matched_cell_ids:
            color_by = self._detect_color_by_from_query(query)
            self.display_cells(matched_cell_ids, color_by=color_by)
        else:
            # If the query was valid but found nothing
            self.clear_search()

    def _detect_color_by_from_query(self, query):
        sample_cell = next(iter(self.cell_measurements.values()))

        # Prefer gene continuous coloring when a gene query is present.
        gene_match = re.search(r'gene\s*:\s*([A-Za-z0-9_.-]+)', query, flags=re.IGNORECASE)
        if gene_match:
            return f"gene:{gene_match.group(1)}"

        # Fallback: find string metadata fields used in categorical queries.
        matches = re.findall(r'(\w+)\s*:', query)
        for key in matches:
            if key in sample_cell and isinstance(sample_cell.get(key), str):
                return key

        return None
    
    def _evaluate_condition(self, condition):
        condition = condition.strip()
        matched = set()

        # Get a sample to check available keys/types
        sample_cid = next(iter(self.cell_measurements))
        sample_data = self.cell_measurements[sample_cid]
        keys = sample_data.keys()
        # -------------------------
        # 🧬 GENE (no operator → select all cells with this gene)
        # -------------------------
        if condition.lower().startswith("gene:") and not any(op in condition for op in [">", "<", "="]):
            gene_name = condition.split(":", 1)[1].strip().capitalize()

            for cid, meas in self.cell_measurements.items():
                val = meas.get("genes", {}).get(gene_name, 0.0)

                if val > 0:   # only include cells expressing the gene
                    matched.add(cid)

            return matched

        # 1. Range: area: 200-500
        m_range = re.match(r'^(\w+)\s*:\s*(\d+)\s*-\s*(\d+)$', condition)
        if m_range:
            param, lo, hi = m_range.groups()
            if param not in keys:
                raise ValueError(f"Unknown field: {param}")
            
            for cid, meas in self.cell_measurements.items():
                val = meas.get(param)
                if isinstance(val, (int, float)) and float(lo) <= val <= float(hi):
                    matched.add(cid)
            return matched

        # 2. Categorical: cell_type: Cardiomyocyte
        # This matches "word: word" or "word: word with spaces"
        m_cat = re.match(r'^(\w+)\s*:\s*([a-zA-Z\s_]+)$', condition)
        if m_cat:
            param, val = m_cat.groups()
            param, val = param.strip(), val.strip().lower()
            
            if param in keys:
                for cid, meas in self.cell_measurements.items():
                    # We check if the stored value matches the search string
                    if str(meas.get(param, "")).lower() == val:
                        matched.add(cid)
                return matched

        # 3. Numeric Comparison: area > 500
        parts = re.split(r'(>=|<=|>|<|=)', condition)
        if len(parts) == 3:
            param, op, rhs = [p.strip() for p in parts]

            # -------------------------
            # 🧬 GENE SUPPORT
            # -------------------------
            is_gene = False

            if param.startswith("gene:"):
                gene = param.split(":", 1)[1]
                is_gene = True

            elif param not in keys:
                raise ValueError(f"Unknown field: {param}")

            rhs_val = float(rhs)

            ops = {
                '>': lambda a, b: a > b,
                '<': lambda a, b: a < b,
                '>=': lambda a, b: a >= b,
                '<=': lambda a, b: a <= b,
                '=': lambda a, b: round(a, 3) == b,
            }

            for cid, meas in self.cell_measurements.items():
                if is_gene:
                    val = self.get_gene_value(cid, gene)
                else:
                    val = meas.get(param)

                if isinstance(val, (int, float)) and ops[op](val, rhs_val):
                    matched.add(cid)

            return matched

        # 4. IDs: 1, 2, 5
        # Only treat as IDs if there are no special operators/colons
        if not any(char in condition for char in ":><="):
            ids = re.findall(r'\d+', condition)
            if ids:
                return set(int(i) for i in ids)

        raise ValueError(f"Invalid format: {condition}")

    def run_stats(self, param):
        if not self.cell_measurements:
            self.search_status.setText("No data")
            return

        # --- 1. Advanced Parsing ---
        # Default values
        raw_params = param
        filter_key = None
        filter_val = None
        bin_count = 10 # Default

        # A. Extract bin count first (e.g., "area for type:A | bins:20")
        if "|" in raw_params:
            parts = raw_params.split("|")
            raw_params = parts[0].strip()
            if "bins" in parts[1]:
                try:
                    bin_count = int(parts[1].split(":")[1].strip())
                except (IndexError, ValueError):
                    bin_count = 10

        # B. Extract Filter (e.g., "area, intensity for type: Cardiomyocyte")
        if " for " in raw_params:
            parts = raw_params.split(" for ")
            raw_params = parts[0].strip() # This is now "area, intensity"
            
            if ":" in parts[1]:
                f_key, f_val = parts[1].split(":")
                filter_key = f_key.strip()
                filter_val = f_val.strip()

        # C. Split the target parameters (e.g., "area, intensity" -> ["area", "intensity"])
        if raw_params == "":
            keys = set()
            for meas in self.cell_measurements.values():
                keys.update(meas.keys())
        else:
            # Split by comma and clean whitespace
            keys = [p.strip() for p in raw_params.split(",")]
        
        keys = list(keys)

        # --- 2. Save File Dialog ---
        path, _ = QFileDialog.getSaveFileName(self, "Save Stats", "stats.txt", "Text Files (*.txt)")
        if not path:
            self.search_status.setText("Stats cancelled")
            return

        output_lines = []
        skip_keys = {"centroid", "score"}

        roi_size_px = self._get_current_roi_size_px()

        # --- 3. Process each Key ---
        for key in (k for k in keys if k not in skip_keys):
            values = []
            for meas in self.cell_measurements.values():
                # Apply Filter
                if filter_key and filter_val:
                    if str(meas.get(filter_key, "")) != filter_val:
                        continue
                
                if key == "ROI size (px)":
                    val = roi_size_px
                else:
                    val = meas.get(key)

                if val is not None:
                    values.append(val)

            if not values:
                continue

            header = f"\n=== {key} ==="
            if filter_key: header += f" (Filter: {filter_key}={filter_val})"
            output_lines.append(header)

            # --- Numeric Analysis ---
            if all(isinstance(v, (int, float)) for v in values):
                arr = np.array(values)
                mean = np.mean(arr)
                
                output_lines.extend([
                    f"Count: {len(arr)}",
                    f"Mean: {mean:.3f}",
                    f"Median: {np.median(arr):.3f}",
                    f"Min/Max: {np.min(arr):.3f} / {np.max(arr):.3f}",
                    f"\nHistogram ({bin_count} bins):"
                ])

                # Use the adjustable bin_count
                hist, bins = np.histogram(arr, bins=bin_count)
                for i in range(len(hist)):
                    output_lines.append(f"{bins[i]:.2f} - {bins[i+1]:.2f}: {hist[i]}")

                # Status update for the first/main param
                if key == keys[0]:
                    self.search_status.setText(f"Analyzed {len(keys)} metrics. {key} mean: {mean:.2f}")

            # --- Categorical Analysis ---
            else:
                counts = Counter(str(v) for v in values)
                total = sum(counts.values())
                for k, c in counts.items():
                    output_lines.append(f"{k}: {c} ({(c/total)*100:.1f}%)")

        # --- 4. Final Write ---
        with open(path, "w") as f:
            f.write("\n".join(output_lines))

    def add_ruler(self):
        if "Ruler" in self.viewer.layers:
            self.viewer.layers.remove("Ruler")      

        ruler = self.viewer.add_shapes(
            name="Ruler",
            shape_type="line",
            edge_color="red",
            face_color="transparent",
            edge_width=2,
        )

        self.viewer.layers.selection.active = ruler
        ruler.mode = "add_line"

        ruler.events.data.connect(self.update_ruler)
    
    def update_ruler(self, event):
        layer = event.source
        if not layer.data:
            return
        
        p1, p2 = layer.data[-1]
        dist = np.linalg.norm(p2 - p1)

        self.info_panel.show_measurement(dist)

    def load_image(self):
        path, _ = QFileDialog.getOpenFileName(
            self,
            "Load Image",
            "",
            "Images (*.tif *.tiff *.png *.jpg *.jpeg *.czi)",
        )
        if not path:
            return

        lower_path = path.lower()

        if lower_path.endswith(".czi"):
            image = czi_imread(path)
        elif lower_path.endswith((".png", ".jpg", ".jpeg")):
            image = iio.imread(path)
        else:
            image = tifffile.imread(path)

        if image.ndim == 4:
            image = image[0]

        # RGB image
        if image.ndim == 3 and image.shape[-1] in (3, 4):
            self.viewer.add_image(image, name="Image", rgb=True, contrast_limits=None)
        else:
            self.viewer.add_image(image, name="Image", contrast_limits=None)

    def export_png(self):
        from qtpy.QtWidgets import QMessageBox

        if not self.viewer.layers:
            return

        choice = QMessageBox.question(
            self,
            "Export PNG",
            "Export selected nucleus only?",
            QMessageBox.Yes | QMessageBox.No | QMessageBox.Cancel
        )

        if choice == QMessageBox.Cancel:
            return

        path, _ = QFileDialog.getSaveFileName(
            self, "Save PNG", "", "PNG (*.png)"
        )
        if not path:
            return

        if choice == QMessageBox.Yes:
            self._export_selected_nucleus_png(path)
        else:
            self._export_full_viewer_png(path)
    
    def _export_full_viewer_png(self, path):
        self.viewer.screenshot(path=path, canvas_only=True)
    

    def _export_selected_nucleus_png(self, path):
        # Determine which layer to pull the selection from
        if "Search Results" in self.viewer.layers:
            labels_layer = self.viewer.layers["Search Results"]
        elif "Nuclei (Cellpose)" in self.viewer.layers:
            labels_layer = self.viewer.layers["Nuclei (Cellpose)"]
        else:
            return

        selected = labels_layer.selected_label
        if selected == 0:
            show_info("Please select a nucleus (label) first.")
            return

        # Backup visibility of all layers
        visibility = {l.name: l.visible for l in self.viewer.layers}

        try:
            # Hide all layers
            for layer in self.viewer.layers:
                layer.visible = False

            # Show the main Image
            if "Image" in self.viewer.layers:
                self.viewer.layers["Image"].visible = True

            # Show ONLY the selected label
            # We create a temporary layer data where only the 'selected' ID exists
            labels_layer.visible = True
            original_data = labels_layer.data
            temp_data = (original_data == selected).astype(np.uint16) * selected
            labels_layer.data = temp_data

            # Take the screenshot
            self.viewer.screenshot(path=path, canvas_only=True)

        finally:
            # Restore original visibility and data
            for layer in self.viewer.layers:
                layer.visible = visibility.get(layer.name, True)
            
            labels_layer.data = original_data

    def export_roi_csv(self):
        from qtpy.QtWidgets import QFileDialog, QMessageBox
        import numpy as np
        import pandas as pd
        from skimage.measure import find_contours

        # Ensure nuclei exist
        if self.nuclei_labels is None or self.cell_measurements is None:
            QMessageBox.warning(self, "No nuclei", "No nuclei segmentation found.")
            return

        # Find the most recently drawn shapes layer (ROI)
        from napari.layers import Shapes

        shapes_layers = [
            layer for layer in self.viewer.layers
            if isinstance(layer, Shapes)
        ]

        if not shapes_layers:
            QMessageBox.warning(self, "No ROI", "Draw a shape on the viewer first.")
            return

        roi_layer = shapes_layers[-1]

        # New version (correct argument name and using actual image shape)
        roi_masks = roi_layer.to_masks(mask_shape=self.nuclei_labels.shape) 

        if len(roi_masks) == 0 or roi_masks[0].size == 0:
            QMessageBox.warning(self, "Empty ROI", "ROI mask is empty.")
            return

        roi_mask = roi_masks[0].astype(bool)  # Use first shape only for now

        # Get nuclei labels that intersect ROI
        labels = self.nuclei_labels
        labels_in_roi = np.unique(labels[roi_mask])
        labels_in_roi = labels_in_roi[labels_in_roi != 0]  # remove background

        if len(labels_in_roi) == 0:
            QMessageBox.information(self, "No nuclei", "No nuclei found inside ROI.")
            return

        # Ask user where to save CSV
        path, _ = QFileDialog.getSaveFileName(
            self, "Export ROI Measurements", "", "CSV (*.csv)"
        )
        if not path:
            return

        # Build export data
        export_rows = []
        for cid in labels_in_roi:
            meas = self.cell_measurements[cid].copy()

            # Compute perimeter (ROI) coordinates only for this nucleus
            mask = labels == cid
            contours = find_contours(mask.astype(np.uint8), 0.5)
            meas["roi_coords"] = [np.fliplr(c) for c in contours]

            meas["cell_id"] = cid
            export_rows.append(meas)

        # Save to CSV
        df = pd.DataFrame(export_rows)
        df.to_csv(path, index=False)

        QMessageBox.information(self, "Exported", f"Exported {len(export_rows)} nuclei to CSV.")


    def run_cellpose_nuclei(self):
        image_layer = self.viewer.layers.selection.active
        if image_layer is None:
            return

        image = image_layer.data
        self._current_image = image

        diameter = self.diameter_spin.value()

        self.cellpose_progress.setVisible(True)
        self.cellpose_progress.setValue(0)

        self.cellpose_timer_label.setVisible(True)
        self.cellpose_timer_label.setText("Elapsed: 0.0 s")

        # Start worker thread
        self.worker = CellposeWorker(
            image,
            diameter)
    
        self.worker.progress.connect(self.update_progress)
        self.worker.finished.connect(self.finish_segmentation)

        # Start elapsed timer
        self._start_time = time.perf_counter()
        self._elapsed_timer = QTimer()
        self._elapsed_timer.timeout.connect(self._update_elapsed)
        self._elapsed_timer.start(100)  # update every 100ms

        self.worker.start()
  
    def annotate_nuclei(self):

        if self.nuclei_labels is None:
            show_info("Run segmentation first.")
            return

        if not self._spatial_initialized:
            if not self._initialize_spatial_engine(): return

        # If search results exist → annotate only those
        if "Search Results" in self.viewer.layers:
            layer = self.viewer.layers["Search Results"]
            ids = np.unique(layer.data)
            ids = ids[ids != 0]
        else:
            ids = self.all_label_ids

        if len(ids) == 0:
            show_info("No nuclei to annotate.")
            return
        
        annotation_mode, ok = QInputDialog.getItem(
            self,
            "Annotation mode",
            "Choose annotation output:",
            ["Genes only (faster)", "Genes + cell type labels", "Cell type labels only"],
            1,
            False
        )
        if not ok:
            return

        annotate_cell_type = annotation_mode in ("Genes + cell type labels", "Cell type labels only")
        annotate_genes = annotation_mode != "Cell type labels only"
        self._annotation_mode = annotation_mode


        self.cellpose_progress.setVisible(True)
        self.cellpose_progress.setValue(0)

        self.annotation_worker = SpatialAnnotationWorker(
            self.nuclei_labels,
            ids,
            self.spatial_engine,
            annotate_cell_type=annotate_cell_type,
            annotate_genes=annotate_genes

        )

        self.annotation_worker.progress.connect(self.update_progress)
        self.annotation_worker.finished.connect(self.finish_annotation)

        self.annotation_worker.start()
   
    def run_cellpose_roi(self):
        from napari.layers import Image, Shapes
        import numpy as np
        import time
        
        # Search for the first available Image layer regardless of name
        image_layer = next((l for l in self.viewer.layers if isinstance(l, Image)), None)
        
        if image_layer is None:
            print("Error: No Image layer found in the viewer.")
            return

        # Keep track of the name for logging
        print(f"Using image layer: {image_layer.name}")
        
        # Now look for the Shapes layer
        shapes_layer = next((l for l in self.viewer.layers if isinstance(l, Shapes)), None)

        if not shapes_layer or len(shapes_layer.data) == 0:
            print("Error: No shapes found to segment.")
            return
        
        selected_indices = list(shapes_layer.selected_data)
        
        if len(selected_indices) > 0:
            idx = selected_indices[0]  
            print(f"Segmenting Selected Sector {idx}...")
        else:
            idx = -1  
            print("No sector selected, segmenting the latest draw.")

        full_label_mask = shapes_layer.to_labels(labels_shape=image_layer.data.shape[:2])
        
        actual_id = idx + 1 if idx >= 0 else len(shapes_layer.data)
        shape_mask_full = (full_label_mask == actual_id)

        poly = shapes_layer.data[idx]
        
        min_y, min_x = np.min(poly, axis=0).astype(int)
        max_y, max_x = np.max(poly, axis=0).astype(int)
        
        # Ensure we stay within image boundaries
        img_h, img_w = image_layer.data.shape[:2]
        min_y, max_y = max(0, min_y), min(img_h, max_y)
        min_x, max_x = max(0, min_x), min(img_w, max_x)

        roi_crop = image_layer.data[min_y:max_y, min_x:max_x]
        mask_crop = shape_mask_full[min_y:max_y, min_x:max_x]

        roi_crop_masked = roi_crop.copy()
        roi_crop_masked[~mask_crop] = 0

        self._roi_offset = (min_y, min_x)

        self._current_image = image_layer.data

        self.cellpose_progress.setVisible(True)
        self.cellpose_progress.setValue(0)
        self.cellpose_timer_label.setVisible(True)
        self.cellpose_timer_label.setText("Elapsed: 0.0 s")

        # 2. Start worker thread (Crucial change: use roi_crop_masked)
        diameter = self.diameter_spin.value()

        self.worker = CellposeWorker(
            roi_crop_masked,
            diameter
        )
        self.worker.progress.connect(self.update_progress)
        
        # 3. Use a specific "finish" function for ROI to stitch the image back
        self.worker.finished.connect(self.finish_roi_segmentation)

        # 4. Start elapsed timer
        self._start_time = time.perf_counter()
        self._elapsed_timer = QTimer()
        self._elapsed_timer.timeout.connect(self._update_elapsed)
        self._elapsed_timer.start(100)

        self.worker.start()
  
    def _update_elapsed(self):
        elapsed = time.perf_counter() - self._start_time
        self.cellpose_timer_label.setText(f"Elapsed: {elapsed:.1f} s")
    

    def update_progress(self, value):
        self.cellpose_progress.setValue(int(value*100))
        self.cellpose_timer_label.setText(f"Progress: {int(value*100)} %")

    def finish_segmentation(self, result):
        masks = result["masks"]

        if "Nuclei (Cellpose)" in self.viewer.layers:
            self.viewer.layers.remove("Nuclei (Cellpose)")

        self.viewer.add_labels(masks, name="Nuclei (Cellpose)")
        self.nuclei_labels = masks
        self.all_label_ids = np.unique(masks)
        self.all_label_ids = self.all_label_ids[self.all_label_ids != 0]  # exclude background
        self.cell_labels = None

        image = getattr(self, "_current_image", None)
        if image.ndim == 3:
            image_gray = image.mean(axis=-1)
        else:
            image_gray = image

        # 🔹 Compute measurements for all nuclei
        from cardioseg.measurements.measurements import compute_cell_measurements
        self.cell_measurements = compute_cell_measurements(masks,image_gray)
        self._restore_non_morphology_metadata_after_cellpose()
        self.color_manager.measurements = self.cell_measurements

        # 🔹 Update Measurements Table
        self.table_panel.set_measurements(self.cell_measurements)

        self.cellpose_timer_label.setText(f"Finished in {result['elapsed']:.2f} s")
        self.cellpose_progress.setValue(100)

        if hasattr(self, "_elapsed_timer"):
            self._elapsed_timer.stop()
        
        if self._csv_restore_notice:
            show_info(self._csv_restore_notice)
            self._csv_restore_notice = ""
    def _restore_non_morphology_metadata_after_cellpose(self):
        source = self._csv_metadata_to_restore
        if not source:
            return

        morphology_keys = {
            "area", "perimeter", "eccentricity", "solidity", "circularity",
            "mean_intensity", "std_intensity", "centroid", "bbox", "major_axis_length",
            "minor_axis_length", "equivalent_diameter", "extent", "orientation",
            "GLCM contrast", "GLCM homogeneity", "Perimeter coordinates", "ROI size (px)"
        }

        target_centroids = {}
        for cid, meas in self.cell_measurements.items():
            centroid = meas.get("centroid")
            if isinstance(centroid, (list, tuple)) and len(centroid) >= 2:
                target_centroids[cid] = np.array([float(centroid[0]), float(centroid[1])], dtype=float)

        used_target_ids = set()
        restored = 0
        for src_meas in source.values():
            src_centroid = src_meas.get("centroid")
            if not isinstance(src_centroid, (list, tuple)) or len(src_centroid) < 2:
                continue

            src_point = np.array([float(src_centroid[0]), float(src_centroid[1])], dtype=float)
            best_cid = None
            best_dist = None
            for cid, tgt_point in target_centroids.items():
                if cid in used_target_ids:
                    continue
                dist = float(np.linalg.norm(src_point - tgt_point))
                if best_dist is None or dist < best_dist:
                    best_dist = dist
                    best_cid = cid

            if best_cid is None:
                continue

            used_target_ids.add(best_cid)
            target = self.cell_measurements.get(best_cid, {})
            for key, value in src_meas.items():
                if key in morphology_keys:
                    continue
                target[key] = value
            restored += 1

        self._csv_metadata_to_restore = None
        self._csv_restore_notice = f"Restored non-morphology metadata for {restored} cells after Cellpose CSV refinement."

    def finish_annotation(self, results):

        self.nucleus_annotations.update(results)

        for cid, annotation in results.items():

            if cid not in self.cell_measurements:
                continue

            # -------------------------
            # EXISTING METADATA
            # -------------------------
            self.cell_measurements[cid]["cell_type"] = annotation.get("cell_type", "Unknown")
            self.cell_measurements[cid]["confidence"] = annotation.get("confidence", 0)
            self.cell_measurements[cid]["score"] = annotation.get("score", 0)

            # GENE STORAGE
            genes = annotation.get("genes", {}) 
            genes_raw = annotation.get("genes_raw", {}) 

            if not isinstance(genes, dict):
                genes = {}
            if not isinstance(genes_raw, dict):
                genes_raw = {}

            self.cell_measurements[cid]["genes"] = genes or {}
            self.cell_measurements[cid]["genes_raw"] = genes_raw or {}

            if genes:
                top_gene = max(genes, key=genes.get)
                top_gene_expr = float(genes.get(top_gene, 0.0))

            else:
                top_gene = "Unknown"
                top_gene_expr = 0.0

            self.cell_measurements[cid]["top_gene"] = top_gene
            self.cell_measurements[cid]["top_gene_expr"] = top_gene_expr

        self._refresh_available_genes()

        self.table_panel.set_measurements(self.cell_measurements)
        self.color_manager.measurements = self.cell_measurements

        mode_label = getattr(self, "_annotation_mode", "Genes + cell type labels")
        show_info(f"Annotated {len(results)} nuclei ({mode_label}).")
  
    def finish_roi_segmentation(self, result):
        roi_masks = result["masks"]
        offset_y, offset_x = self._roi_offset
        
        if self.nuclei_labels is None:
            self.nuclei_labels = np.zeros(self._current_image.shape[:2], dtype=np.uint32)
        
        current_max = np.max(self.nuclei_labels)
        new_labels = np.where(roi_masks > 0, roi_masks + current_max, 0).astype(np.uint32)

        h, w = roi_masks.shape
        self.nuclei_labels[offset_y:offset_y+h, offset_x:offset_x+w] = new_labels

        self.all_label_ids = np.unique(self.nuclei_labels)
        self.all_label_ids = self.all_label_ids[self.all_label_ids != 0]
        
        self.finish_segmentation({"masks": self.nuclei_labels, "elapsed": result["elapsed"]})
        
        show_info(f"ROI Segmentation complete. Added {np.max(roi_masks) - current_max} nuclei.")
    
    def get_gene_value(self, cid, gene):
        return self.cell_measurements.get(cid, {}).get("genes", {}).get(gene, 0.0)

    def save_measurements(self):
        print("START")
        if not self.cell_measurements:
            print("No measurements")
            show_info("No measurements found to export.")    
            return

        # 1. Determine which labels to export
        # Check if a search is currently active
        if "Search Results" in self.viewer.layers:
            search_layer = self.viewer.layers["Search Results"]
            # Get unique IDs from the search layer (excluding background 0)
            ids_to_export = np.unique(search_layer.data)
            ids_to_export = ids_to_export[ids_to_export != 0]
            
            mode_title = f"Save Search Results ({len(ids_to_export)} cells)"
        else:
            # Export everything
            ids_to_export = list(self.cell_measurements.keys())
            mode_title = f"Save All Measurements ({len(ids_to_export)} cells)"

        if len(ids_to_export) == 0:
            show_info("No nuclei found in the current view to export.")
            return
        
        # 2. Get file path
        path, _ = QFileDialog.getSaveFileName(
            self, mode_title, "", "CSV (*.csv)"
        )
        if not path:
            return

        # 3. Filter the dictionary
        filtered_measurements = {
            cid: dict(self.cell_measurements[cid])
            for cid in ids_to_export 
            if cid in self.cell_measurements
        }

        print("\n=== DEBUG: FILTERED SAMPLE ===")

        if filtered_measurements:
            sample_id = next(iter(filtered_measurements))
            sample = filtered_measurements[sample_id]

            print(f"Sample cell ID: {sample_id}")
            print(f"Total fields: {len(sample)}")

            for k, v in sample.items():
                if isinstance(v, dict):
                    print(f"{k}: dict ({len(v)} keys)")
                else:
                    print(f"{k}: {type(v).__name__} → {v}")
        else:
            print("No filtered measurements found.")

        print("\n=== DEBUG: FULL DATASET SAMPLE ===")

        if self.cell_measurements:
            sample_id = next(iter(self.cell_measurements))
            sample = self.cell_measurements[sample_id]

            print(f"Sample cell ID: {sample_id}")
            print(f"Total fields: {len(sample)}")

            for k, v in sample.items():
                if isinstance(v, dict):
                    print(f"{k}: dict ({len(v)} keys)")
                else:
                    print(f"{k}: {type(v).__name__} → {v}")
        else:
            print("No measurements in full dataset.")

        include_shape, include_texture = self._prompt_optional_export_features()
        if include_shape is None:
            return

        self._attach_optional_export_features(
            filtered_measurements,
            include_shape=include_shape,
            include_texture=include_texture
        )

        has_gene_data = any(
            isinstance(meas.get("genes"), dict) and len(meas.get("genes")) > 0
            for meas in filtered_measurements.values()
        ) or any(
            isinstance(meas.get("genes_raw"), dict) and len(meas.get("genes_raw")) > 0
            for meas in filtered_measurements.values()
        )

        selected_fields = self._prompt_fields_for_export(filtered_measurements)
        if selected_fields is None:
            return
       
        if "Perimeter coordinates" in selected_fields:
            self._attach_perimeter_coordinates(filtered_measurements)
        
        if "ROI size (px)" in selected_fields:
            self._attach_roi_size(filtered_measurements)

        if has_gene_data:
            gene_mode, max_genes, raw_as_integer = self._prompt_gene_export_options()
            if gene_mode is None:
                return
        else:
            gene_mode, max_genes, raw_as_integer = "None", 0, False

        flattened_measurements = {}
        gene_columns = self._choose_gene_columns(filtered_measurements, gene_mode, max_genes)

        for cid, meas in filtered_measurements.items():
            row = {}

            for field in selected_fields:

                value = meas.get(field, None)
                if isinstance(value, (list, tuple)):
                    value = str(value)
                row[field] = value

            for gene in gene_columns.get("normalized", []):
                row[f"gene_norm_{gene}"] = meas.get("genes", {}).get(gene, 0.0)

            for gene in gene_columns.get("raw", []):
                raw_value = meas.get("genes_raw", {}).get(gene, 0.0)
                row[f"gene_raw_{gene}"] = self._format_raw_gene_value(raw_value, raw_as_integer)

            flattened_measurements[cid] = row

        export_cell_measurements(flattened_measurements, path)
        show_info(f"Successfully exported to {path}")
    
    def _prompt_fields_for_export(self, filtered_measurements):
        exportable_fields = set()

        for meas in filtered_measurements.values():
            for key, value in meas.items():
                if isinstance(value, dict):
                    continue
                exportable_fields.add(key)

        exportable_fields.add("Perimeter coordinates")
        exportable_fields.add("ROI size (px)")

        ordered_fields = sorted(exportable_fields)
        default_fields = {
            "area", "perimeter", "eccentricity", "solidity",
            "circularity", "mean_intensity", "std_intensity",
            "centroid", "cell_type", "confidence", "score", "top_gene", "top_gene_expr"
        }

        dialog = QDialog(self)
        dialog.setWindowTitle("Select parameters to export")
        layout = QVBoxLayout(dialog)

        list_widget = QListWidget(dialog)
        list_widget.setSelectionMode(QAbstractItemView.NoSelection)

        for field in ordered_fields:
            item = QListWidgetItem(field)
            item.setFlags(item.flags() | Qt.ItemIsUserCheckable)
            item.setCheckState(Qt.Checked if field in default_fields else Qt.Unchecked)
            list_widget.addItem(item)

        layout.addWidget(list_widget)
        buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel, parent=dialog)
        buttons.accepted.connect(dialog.accept)
        buttons.rejected.connect(dialog.reject)
        layout.addWidget(buttons)

        if dialog.exec() != QDialog.Accepted:
            return None

        selected = [
            list_widget.item(i).text()
            for i in range(list_widget.count())
            if list_widget.item(i).checkState() == Qt.Checked
        ]

        if not selected:
            show_info("No parameters selected for export.")
            return None

        return selected
    
    def _get_current_roi_size_px(self):
        from napari.layers import Shapes

        if self.nuclei_labels is None:
            return None

        shapes_layers = [layer for layer in self.viewer.layers if isinstance(layer, Shapes)]
        if not shapes_layers:
            return None

        roi_layer = shapes_layers[-1]
        roi_masks = roi_layer.to_masks(mask_shape=self.nuclei_labels.shape)
        if len(roi_masks) == 0:
            return None

        return int(np.count_nonzero(roi_masks[0]))
    
    def _attach_roi_size(self, filtered_measurements):
        from napari.layers import Shapes

        if self.nuclei_labels is None:
            show_info("ROI size export requires nuclei labels; leaving values empty.")
            return

        shapes_layers = [layer for layer in self.viewer.layers if isinstance(layer, Shapes)]
        if not shapes_layers:
            show_info("ROI size export requested but no ROI shape was found; leaving values empty.")
            return

        roi_layer = shapes_layers[-1]
        roi_masks = roi_layer.to_masks(mask_shape=self.nuclei_labels.shape)
        if len(roi_masks) == 0:
            show_info("ROI size export requested but ROI mask is empty; leaving values empty.")
            return

        roi_size = int(np.count_nonzero(roi_masks[0]))
        for meas in filtered_measurements.values():
            meas["ROI size (px)"] = roi_size
   
    def _attach_perimeter_coordinates(self, filtered_measurements):
        if self.nuclei_labels is None:
            show_info("Perimeter coordinates require nuclei labels; skipping this export field.")
            return

        from skimage.measure import find_contours

        for cid, meas in filtered_measurements.items():
            mask = self.nuclei_labels == cid
            contours = find_contours(mask.astype(np.uint8), 0.5)
            meas["Perimeter coordinates"] = [np.fliplr(c).tolist() for c in contours]
    
    def _prompt_optional_export_features(self):
        dialog = QDialog(self)
        dialog.setWindowTitle("Optional export features")
        layout = QVBoxLayout(dialog)

        shape_checkbox = QCheckBox("Include shape extras (circularity, solidity, standard_intensity)")
        texture_checkbox = QCheckBox("Include texture extras (GLCM contrast, GLCM homogeneity)")

        shape_checkbox.setChecked(False)
        texture_checkbox.setChecked(False)

        layout.addWidget(shape_checkbox)
        layout.addWidget(texture_checkbox)

        buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel, parent=dialog)
        buttons.accepted.connect(dialog.accept)
        buttons.rejected.connect(dialog.reject)
        layout.addWidget(buttons)

        if dialog.exec() != QDialog.Accepted:
            return None, None

        return shape_checkbox.isChecked(), texture_checkbox.isChecked()

    def _attach_optional_export_features(self, filtered_measurements, include_shape=False, include_texture=False):
        if not include_shape and not include_texture:
            return

        if self.nuclei_labels is None:
            show_info("Optional features require nuclei labels; skipping optional feature extraction.")
            return

        image = getattr(self, "_current_image", None)
        if image is None and "Image" in self.viewer.layers:
            image = self.viewer.layers["Image"].data
        if image is None:
            show_info("No source image available; skipping optional feature extraction.")
            return

        if image.ndim == 3:
            image_gray = image.mean(axis=-1)
        else:
            image_gray = image

        from cardioseg.measurements.measurements import compute_export_features
        optional_features = compute_export_features(
            self.nuclei_labels,
            image_gray,
            include_shape=include_shape,
            include_texture=include_texture,
        )

        for cid, row in optional_features.items():
            if cid in filtered_measurements:
                filtered_measurements[cid].update(row)

    def _prompt_gene_export_options(self):
        mode, ok = QInputDialog.getItem(
            self,
            "Gene export",
            "Gene values to export:",
            ["None", "Normalized", "Raw", "Both", "All (Normalized)", "All (Raw)", "All (Both)"],
            3,
            False
        )
        if not ok:
            return None, None, None

        if mode == "None":
            return mode, 0, False

        if mode.startswith("All"):
            max_genes = None
        else:
            max_genes, ok = QInputDialog.getInt(
                self,
                "Gene export",
                "Maximum number of genes to export per mode:",
                50,
                1,
                2000,
                1
            )
            if not ok:
                return None, None, None

        mode_map = {
            "All (Normalized)": "Normalized",
            "All (Raw)": "Raw",
            "All (Both)": "Both",
        }
        normalized_mode = mode_map.get(mode, mode)

        raw_as_integer = False
        if normalized_mode in ("Raw", "Both"):
            raw_format, ok = QInputDialog.getItem(
                self,
                "Raw gene format",
                "Raw gene values format:",
                ["Integer counts", "Decimal"],
                0,
                False
            )
            if not ok:
                return None, None, None
            raw_as_integer = raw_format == "Integer counts"

        return normalized_mode, max_genes, raw_as_integer

    def _choose_gene_columns(self, filtered_measurements, mode, max_genes):
        if mode == "None":
            return {"normalized": [], "raw": []}

        def top_genes_from(field_name):
            totals = Counter()
            for meas in filtered_measurements.values():
                genes = meas.get(field_name, {})
                if not isinstance(genes, dict):
                    continue
                for gene, value in genes.items():
                    totals[gene] += float(value)

            ranked_genes = [g for g, _ in totals.most_common()]
            if max_genes is None:
                return ranked_genes
            return ranked_genes[:max_genes]

        cols = {"normalized": [], "raw": []}
        if mode in ("Normalized", "Both"):
            cols["normalized"] = top_genes_from("genes")
        if mode in ("Raw", "Both"):
            cols["raw"] = top_genes_from("genes_raw")

        return cols

    def _format_raw_gene_value(self, value, raw_as_integer):
        if not raw_as_integer:
            return value

        try:
            return int(round(float(value)))
        except (TypeError, ValueError):
            return 0

    def _on_mouse_click(self, viewer, event):
        if event.type != "mouse_press":
            return

        if self.nuclei_labels is None or not self.cell_measurements:
            return
        
        pos = np.round(viewer.cursor.position).astype(int)
        
        y, x = pos[-2], pos[-1]

        if (y < 0 or x < 0 or
            y >= self.nuclei_labels.shape[0] or
            x >= self.nuclei_labels.shape[1]):
            return
        
        label_id = int(self.nuclei_labels[y, x])

        if label_id == 0:
            self.clear_highlight()
            self.info_panel.clear()
            return
        else:
            self.select_nucleus(label_id)
  
    def select_nucleus(self, cell_id):

        self.highlight_cell(cell_id)

        meas = self.cell_measurements.get(cell_id)
        if meas:
            self.info_panel.show_cell_info(cell_id, meas)

    def highlight_cell(self, cell_id, layer_name = "Nuclei (Cellpose)",):
        if layer_name not in self.viewer.layers:
            return
        
        labels_layer = self.viewer.layers[layer_name]
        
        labels_layer.selected_label = cell_id

    def render(self, filtered, color_info, legend):
        # --- Step 1: Remove old layer
        if "Search Results" in self.viewer.layers:
            self.viewer.layers.remove("Search Results")

        # --- Step 2: Dim base layer
        if "Nuclei (Cellpose)" in self.viewer.layers:
            self.viewer.layers["Nuclei (Cellpose)"].opacity = 0

        # --- Step 3: Render based on type
        if legend["type"] == "categorical":

            layer = self.viewer.add_labels(
                filtered,
                name="Search Results"
            )

            layer.color_mode = "direct"
            layer.color = color_info
            layer.opacity = 1.0

        elif legend["type"] == "continuous":

            value_img = np.full_like(filtered, np.nan, dtype=float)
            for lab, val in color_info.items():
                value_img[filtered == lab] = val

            mask = filtered > 0

            layer = self.viewer.add_image(
                value_img,
                name="Search Results",
                colormap=legend.get("colormap", "fire"),
                contrast_limits=(legend["min"], legend["max"]),
                blending="minimum",
            )

            # make background invisible
            layer.data[~mask] = legend["min"]
            layer.opacity = 1.0

        else:
            # plain highlight (no coloring)
            layer = self.viewer.add_labels(
                filtered,
                name="Search Results"
            )
            layer.opacity = 1.0

        # --- Step 4: Bring to front
        self.viewer.layers.move(
            self.viewer.layers.index(layer),
            len(self.viewer.layers) - 1
        )

        layer.refresh()
        self.update_legend(legend)
        
    def display_cells(self, cell_ids=None, color_by=None):
        if cell_ids is None:
            cell_ids = list(self.cell_measurements.keys())
        
        cell_ids = [int(c) for c in cell_ids]

        label_ids = np.unique(self.nuclei_labels)

        if color_by is None:
            filtered = np.isin(self.nuclei_labels, cell_ids) * self.nuclei_labels
            self.render(filtered, None, {"type": "none"})
            return

        mode = self.color_manager.get_mode(cell_ids, color_by)

        if "Nuclei (Cellpose)" in self.viewer.layers:
            self.viewer.layers["Nuclei (Cellpose)"].opacity = 0.0

        if mode == "categorical":
            f, c, l = self.color_manager.categorical(cell_ids, color_by, self.nuclei_labels)
        elif mode == "continuous":
            f, c, l = self.color_manager.continuous(cell_ids, color_by, self.nuclei_labels)
            if color_by.startswith("gene:"):
                l["label"] = color_by
        elif mode == "gene":
            gene = color_by.split(":")[1]
            f, c, l = self.color_manager.gene(cell_ids, gene, self.nuclei_labels)
        elif mode == "top_gene":
            f, c, l = self.color_manager.top_gene(cell_ids, self.nuclei_labels)

        self.render(f, c, l)
   
    def update_legend(self, legend):
        if hasattr(self, "_legend_widget") and self._legend_widget is not None:
            try:
                self.viewer.window.remove_dock_widget(self._legend_widget)
            except Exception:
                pass
            self._legend_widget = None
        
        if legend.get("type") in ("none", None):
            return

        container = QWidget()
        main_layout = QVBoxLayout()

        if legend["type"] == "categorical":
            scroll = QScrollArea()
            scroll.setWidgetResizable(True)

            inner = QWidget()
            layout = QVBoxLayout()

            if legend.get("field") == "cell_type":
                present = set(legend.get("mapping", {}).keys())
                items = [
                    (cell_type, self.cell_type_colors.get(cell_type, self.cell_type_colors["Unknown"]))
                    for cell_type in self.cell_type_colors
                    if cell_type in present
                ]
            else:
                items = []
                for cat, lab in legend["mapping"].items():
                    colors = legend.get("colors", {})
                    color = colors.get(lab, colors.get(cat, (0.7, 0.7, 0.7, 1.0)))
                    items.append((cat, color))

            for cat, color in items:
                if isinstance(color, QColor):
                    color = color.getRgbF()
                
                # Ensure color is in 0-1 range
                if isinstance(color, (tuple, list)) and len(color) >= 3:
                    r, g, b = color[0], color[1], color[2]
                    # Handle both 0-1 and 0-255 ranges
                    if max(r, g, b) > 1:
                        r, g, b = r/255, g/255, b/255
                    
                    # Convert to 0-255 for CSS swatches.
                    r_css = int(round(r * 255))
                    g_css = int(round(g * 255))
                    b_css = int(round(b * 255))
                else:
                    r_css = g_css = b_css = 128

                row = QHBoxLayout()

                swatch = QLabel()
                swatch.setFixedSize(14, 14)
                swatch.setStyleSheet(
                    f"background-color: rgb({r_css},"
                    f"{g_css},"
                    f"{b_css});"
                    "border: 1px solid black;"
                )

                text = QLabel(str(cat))

                row.addWidget(swatch)
                row.addWidget(text)
                row.addStretch()

                layout.addLayout(row)

            inner.setLayout(layout)
            scroll.setWidget(inner)
            main_layout.addWidget(scroll)

        elif legend["type"] == "continuous":
            main_layout.addWidget(QLabel(f"{legend.get('label','Value')}"))

            # simple gradient bar
            gradient = QLabel()
            gradient.setFixedHeight(20)
            gradient.setStyleSheet("""
                background: qlineargradient(
                    x1:0, y1:0, x2:1, y2:0,
                    stop:0 rgb(255, 255, 0),
                    stop:0.5 rgb(255, 80, 0),
                    stop:1 rgb(180, 0, 0)
                );
                border: 1px solid black;
            """)

            main_layout.addWidget(gradient)
            main_layout.addWidget(QLabel(f"{legend['min']:.2f} → {legend['max']:.2f}"))

        container.setLayout(main_layout)

        self._legend_widget = container
        self.viewer.window.add_dock_widget(container, area="right")

    def clear_highlight(self):
        if "Nuclei (Cellpose)" in self.viewer.layers:
            self.viewer.layers["Nuclei (Cellpose)"].selected_label = 0
 
    def clear_search(self):
        self.search_input.clear()

        if "Search Results" in self.viewer.layers:
            self.viewer.layers.remove("Search Results")

        if "Nuclei (Cellpose)" in self.viewer.layers:
            layer = self.viewer.layers["Nuclei (Cellpose)"]
            layer.opacity = 1.0
            layer.selected_label = 0

        self.update_legend({"type": "none"})
        self.search_status.setText("")