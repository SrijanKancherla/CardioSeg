from cardioseg.assignment.assign import assign_nuclei_to_cells
from cardioseg.measurements.measurements import compute_cell_measurements
import numpy as np


def connect_layer_events(viewer, info_panel, control_panel, table_panel):

    # -------------------------
    # Helper: ensure measurements exist
    # -------------------------
    def ensure_measurements():
        if control_panel.cell_measurements:
            return

        image = viewer.layers["Image"].data

        if control_panel.nuclei_labels is None:
            control_panel.cell_measurements = compute_cell_measurements(
                control_panel.cell_labels,
                image,
                nuclei_assignment=None,
            )
        else:
            nuclei_map = assign_nuclei_to_cells(
                control_panel.nuclei_labels,
                control_panel.cell_labels,
            )
            control_panel.cell_measurements = compute_cell_measurements(
                control_panel.cell_labels,
                image,
                nuclei_map,
            )

    # -------------------------
    # Mouse click (select cell)
    # -------------------------
    @viewer.mouse_drag_callbacks.append
    def on_click(viewer, event):
        if event.type != "mouse_press":
            return

        if "Cells (Cellpose)" not in viewer.layers:
            return

        labels_layer = viewer.layers["Cells (Cellpose)"]
        y, x = map(int, event.position)

        if (
            y < 0 or x < 0 or
            y >= labels_layer.data.shape[0] or
            x >= labels_layer.data.shape[1]
        ):
            return

        cell_label = int(labels_layer.data[y, x])

        if cell_label == 0:
            labels_layer.selected_label = 0
            info_panel.clear()
            return

        ensure_measurements()

        # VISUAL SELECTION
        labels_layer.selected_label = cell_label

        info = control_panel.cell_measurements.get(cell_label)
        info_panel.show_cell_info(cell_label, info)

        table_panel.highlight_cell(cell_label)

        # -------------------------
        # Ruler handling
        # -------------------------
        def on_shapes_change(event):
            layer = event.source
            if layer.name == "Ruler":
                update_ruler_info(layer, info_panel)

        viewer.layers.events.inserted.connect(
            lambda e: (
                e.value.events.data.connect(on_shapes_change)
                if e.value.name == "Ruler"
                else None
            )
        )
 
def update_ruler_info(layer, info_panel):
    if len(layer.data) == 0:
        info_panel.clear_measurement()
        return

    (y1, x1), (y2, x2) = layer.data[-1]
    dist = np.sqrt((y2 - y1) ** 2 + (x2 - x1) ** 2)

    info_panel.show_measurement(dist)
