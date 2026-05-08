from qtpy.QtWidgets import QWidget, QFormLayout, QLabel

import numpy as np

from qtpy.QtWidgets import QWidget, QFormLayout, QLabel

class InfoPanel(QWidget):
    def __init__(self, viewer):
        super().__init__()
        self.viewer = viewer

        layout = QFormLayout()

        self.id_label = QLabel("—")
        self.area_label = QLabel("—")
        self.centroid_label = QLabel("—")
        self.measure_label = QLabel("Measurement: —")

        layout.addRow(self.measure_label)
        layout.addRow("Cell ID", self.id_label)
        layout.addRow("Area", self.area_label)
        layout.addRow("Centroid", self.centroid_label)

        self.fields = {}
        for name in [
            "Eccentricity",
            "Top Gene",
        ]:
            label = QLabel("—")
            layout.addRow(name, label)
            self.fields[name] = label

        self.setLayout(layout)

    def clear(self):
        self.id_label.setText("—")
        self.area_label.setText("—")
        self.centroid_label.setText("—")
        for label in self.fields.values():
            label.setText("—")
        self.clear_measurement()

    def show_cell_info(self, cell_id, info):
        if info is None:
            self.clear()
            return

        cell_type = info.get("cell_type")
        confidence = info.get("confidence", 0)

        if cell_type:
            self.id_label.setText(f"{cell_id} ({cell_type}, {confidence:.2f})")
        else:
            self.id_label.setText(str(cell_id))

        self.area_label.setText(f"{info['area']:.1f} px²")
        y, x = info["centroid"]
        self.centroid_label.setText(f"({y:.1f}, {x:.1f})")

        self.fields["Eccentricity"].setText(f"{info['eccentricity']:.3f}")

        genes = info.get("genes") if isinstance(info, dict) else None
        if isinstance(genes, dict) and genes:
            top_gene = max(genes, key=genes.get)
            top_expr = genes.get(top_gene, 0.0)
            self.fields["Top Gene"].setText(f"{top_gene} ({top_expr:.4f})")
        else:
            top_gene = info.get("top_gene", "Unknown") if isinstance(info, dict) else "Unknown"
            top_expr = info.get("top_gene_expr", None) if isinstance(info, dict) else None
            if top_expr is not None and top_gene not in (None, "Unknown"):
                self.fields["Top Gene"].setText(f"{top_gene} ({float(top_expr):.4f})")
            else:
                self.fields["Top Gene"].setText(str(top_gene))

    def show_measurement(self, dist_px):
        self.measure_label.setText(f"Measurement: {dist_px:.1f} px")

    def clear_measurement(self):
        self.measure_label.setText("Measurement: —")

