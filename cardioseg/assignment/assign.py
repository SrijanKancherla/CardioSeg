from skimage.measure import regionprops

def assign_nuclei_to_cells(nuclei_labels, cell_labels):
    """
    Assign each nucleus to a cell using centroid lookup.

    Returns:
        dict: nucleus_label -> cell_label (or None)
    """
    assignments = {}

    for prop in regionprops(nuclei_labels):
        y, x = map(int, prop.centroid)

        try:
            cell_label = cell_labels[y, x]
        except IndexError:
            cell_label = 0

        assignments[prop.label] = int(cell_label) if cell_label > 0 else None

    return assignments
