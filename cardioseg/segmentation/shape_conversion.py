import numpy as np
from skimage.draw import polygon

def shapes_to_labels(shapes, image_shape, label_id):
    mask = np.zeros(image_shape, dtype=np.int32)
    for shape in shapes:
        rr, cc = polygon(shape[:, 0], shape[:, 1])
        mask[rr, cc] = label_id
    return mask
