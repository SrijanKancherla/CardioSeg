import napari
from qtpy.QtWidgets import QApplication
from pathlib import Path

from .dock_widgets import ControlPanel
from .info_panel import InfoPanel
from .layer_events import connect_layer_events
from .measurements_table import MeasurementTable

def load_stylesheet():
    app = QApplication.instance()
    if app is None:
        return
    style_path = Path(__file__).parent / "styles.qss"
    if style_path.exists():
        with open(style_path) as f:
            app.setStyleSheet(f.read())

def launch_viewer():
    viewer = napari.Viewer(
        title="CardioSeg",
        axis_labels=("Y", "X"),
    )

    load_stylesheet()
    
    info_panel = InfoPanel(viewer)
    table_panel = MeasurementTable(viewer)
    control_panel = ControlPanel(viewer, info_panel, table_panel)

    viewer.window.add_dock_widget(control_panel, area="right", name="Controls")
    viewer.window.add_dock_widget(info_panel, area="right", name="Info")
    viewer.window.add_dock_widget(table_panel, area="right", name="Measurements")

    connect_layer_events(viewer, info_panel, control_panel, table_panel)

    return viewer
