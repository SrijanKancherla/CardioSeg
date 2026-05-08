#polygon logic

import numpy as np
from matplotlib.path import Path


def bins_inside_polygon(polygon, bin_coords):
    """
    polygon: Nx2 array (x, y)
    bin_coords: Mx2 array (x, y)

    Returns indices of bins inside polygon.
    """

    path = Path(polygon)
    inside = path.contains_points(bin_coords)

    return np.where(inside)[0]
