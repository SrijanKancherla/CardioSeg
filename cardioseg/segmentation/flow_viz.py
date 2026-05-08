import numpy as np

def flow_to_rgb(dy, dx):
    """
    Convert flow fields (dy, dx) to an RGB image for visualization.
    """
    magnitude = np.sqrt(dy**2 + dx**2)
    magnitude = np.clip(magnitude / np.percentile(magnitude, 99), 0, 1)

    angle = np.arctan2(dy, dx)  # radians, not degrees!

    rgb = np.zeros((*dy.shape, 3), dtype=np.float32)
    rgb[..., 0] = (np.cos(angle) + 1) / 2 
    rgb[..., 1] = (np.sin(angle) + 1) / 2
    rgb[..., 2] = magnitude

    return rgb
