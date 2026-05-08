from qtpy.QtWidgets import QWidget, QVBoxLayout, QTableWidget, QTableWidgetItem
from qtpy.QtCore import Qt
import numpy as np

class MeasurementTable(QWidget):
    def __init__(self, viewer):
        super().__init__()
        self.viewer = viewer

        layout = QVBoxLayout()

        self.table = QTableWidget()
        self.table.setColumnCount(7)
        self.table.setHorizontalHeaderLabels([
            "Cell ID",
            "Area",
            "Perimeter",
            "Eccentricity",
            "Mean Intensity",
            "Centroid (Y, X)",
            "Cell Type"
        ])

        self.table.setEditTriggers(QTableWidget.NoEditTriggers)
        self.table.setSelectionBehavior(QTableWidget.SelectRows)
        self.table.setSelectionMode(QTableWidget.SingleSelection)
        self.table.cellClicked.connect(self._on_row_selected)

        layout.addWidget(self.table)
        self.setLayout(layout)

        # mapping cell IDs → row indices
        self.cell_id_to_row = {}

        self._last_id = None

    def set_measurements(self, measurements: dict):
        self.table.setRowCount(0)
        self.cell_id_to_row.clear()

        for row, (cell_id, info) in enumerate(measurements.items()):
            self.table.insertRow(row)
            self.cell_id_to_row[cell_id] = row

            values = [
                cell_id,
                f"{info['area']:.1f}",
                f"{info['perimeter']:.1f}",
                f"{info['eccentricity']:.3f}",
                f"{info['mean_intensity']:.2f}",
                f"({info['centroid'][0]:.1f}, {info['centroid'][1]:.1f})",
                info.get("cell_type", "—"),
                info.get("confidence", None)
            ]

            for col, value in enumerate(values):
                item = QTableWidgetItem(str(value))
                item.setTextAlignment(Qt.AlignCenter)
                self.table.setItem(row, col, item)

        self.table.resizeColumnsToContents()

    def highlight_cell(self, cell_id):
        """Select and scroll to a row programmatically."""
        if cell_id not in self.cell_id_to_row:
            return

        row = self.cell_id_to_row[cell_id]
        self.table.selectRow(row)
        self.table.scrollToItem(
            self.table.item(row, 0),
            QTableWidget.PositionAtCenter
        )

    def _on_row_selected(self, row, column):
        cell_id = int(self.table.item(row, 0).text())

        if self._last_id == cell_id:
            self.table.clearSelection()
            self._last_id = None
            
            if "Nuclei (Cellpose)" in self.viewer.layers:
                self.viewer.layers["Nuclei (Cellpose)"].selected_label = 0
                self.viewer.reset_view()
            
            return
        
        self._last_id = cell_id
        self._highlight_in_viewer(cell_id)

    def _highlight_in_viewer(self, cell_id):
        if "Nuclei (Cellpose)" not in self.viewer.layers:
            return

        labels = self.viewer.layers["Nuclei (Cellpose)"]
        labels.selected_label = cell_id

        # ---- CENTER + ZOOM ON NUCLEUS ----
        mask = labels.data == cell_id
        if not np.any(mask):
            return

        coords = np.column_stack(np.nonzero(mask))
        centroid = coords.mean(axis=0)

        self.viewer.camera.center = tuple(centroid)
        self.viewer.camera.zoom = 4


# class MeasurementTable(QWidget):
#     def __init__(self, viewer):
#         super().__init__()
#         self.viewer = viewer

#         layout = QVBoxLayout()

#         self.table = QTableWidget()
#         self.table.setColumnCount(7)
#         self.table.setHorizontalHeaderLabels([
#             "Cell ID",
#             "Area",
#             "Perimeter",
#             "Eccentricity",
#             "Mean Intensity",
#             "Centroid (Y, X)"
#         ])

#         self.table.setEditTriggers(QTableWidget.NoEditTriggers)
#         self.table.setSelectionBehavior(QTableWidget.SelectRows)
#         self.table.setSelectionMode(QTableWidget.SingleSelection)
#         self.table.cellClicked.connect(self._on_row_selected)

#         layout.addWidget(self.table)
#         self.setLayout(layout)

#         # mapping cell IDs → row indices
#         self.cell_id_to_row = {}

#         self._last_id = None

#     def set_measurements(self, measurements: dict):
#         self.table.setRowCount(0)
#         self.cell_id_to_row.clear()

#         for row, (cell_id, info) in enumerate(measurements.items()):
#             self.table.insertRow(row)
#             self.cell_id_to_row[cell_id] = row

#             values = [
#                 cell_id,
#                 f"{info['area']:.1f}",
#                 f"{info['perimeter']:.1f}",
#                 f"{info['eccentricity']:.3f}",
#                 f"{info['mean_intensity']:.2f}",
#                 f"({info['centroid'][0]:.1f}, {info['centroid'][1]:.1f})"
#             ]

#             for col, value in enumerate(values):
#                 item = QTableWidgetItem(str(value))
#                 item.setTextAlignment(Qt.AlignCenter)
#                 self.table.setItem(row, col, item)

#         self.table.resizeColumnsToContents()

#     def highlight_cell(self, cell_id):
#         """Select and scroll to a row programmatically."""
#         if cell_id not in self.cell_id_to_row:
#             return

#         row = self.cell_id_to_row[cell_id]
#         self.table.selectRow(row)
#         self.table.scrollToItem(
#             self.table.item(row, 0),
#             QTableWidget.PositionAtCenter
#         )

#     def _on_row_selected(self, row, column):
#         cell_id = int(self.table.item(row, 0).text())

#         if self._last_id == cell_id:
#             self.table.clearSelection()
#             self._last_id = None
            
#             if "Nuclei (Cellpose)" in self.viewer.layers:
#                 self.viewer.layers["Nuclei (Cellpose)"].selected_label = 0
#                 self.viewer.reset_view()
            
#             return
        
#         self._last_id = cell_id
#         self._highlight_in_viewer(cell_id)

#     def _highlight_in_viewer(self, cell_id):
#         if "Nuclei (Cellpose)" not in self.viewer.layers:
#             return

#         labels = self.viewer.layers["Nuclei (Cellpose)"]
#         labels.selected_label = cell_id

#         # ---- CENTER + ZOOM ON NUCLEUS ----
#         mask = labels.data == cell_id
#         if not np.any(mask):
#             return

#         coords = np.column_stack(np.nonzero(mask))
#         centroid = coords.mean(axis=0)

#         self.viewer.camera.center = tuple(centroid)
#         self.viewer.camera.zoom = 4
