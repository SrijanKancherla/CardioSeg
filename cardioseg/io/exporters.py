import pandas as pd
import numpy as np
import csv
import sys

def export_cell_measurements(measurements, path):

    flat_measurements = {}
    for cell_id, info in measurements.items():
        flat_info = info.copy()
        if "roi_coords" in flat_info:
            flat_info["roi_coords"] = str(flat_info["roi_coords"])
        flat_measurements[cell_id] = flat_info
    
    if _is_windows():
        _export_cell_measurements_windows_fast(flat_measurements, path)
        return

    df = pd.DataFrame.from_dict(flat_measurements, orient="index")
    df.index.name = "cell_id"
    df.to_csv(path)

def _is_windows():
    return sys.platform.startswith("win")


def _export_cell_measurements_windows_fast(flat_measurements, path):
    if not flat_measurements:
        with open(path, "w", newline="", encoding="utf-8") as f:
            f.write("cell_id\n")
        return

    fieldnames = {"cell_id"}
    for row in flat_measurements.values():
        fieldnames.update(row.keys())

    ordered_fields = ["cell_id"] + sorted(f for f in fieldnames if f != "cell_id")

    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=ordered_fields,
            lineterminator="\n",
            extrasaction="ignore"
        )
        writer.writeheader()
        for cell_id, row in flat_measurements.items():
            out_row = {"cell_id": cell_id}
            out_row.update(row)
            writer.writerow(out_row)
