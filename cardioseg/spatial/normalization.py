#count normalization

import numpy as np


def normalize_vector(counts):
    """
    Total count normalize + log1p
    """

    total = counts.sum()

    if total == 0:
        return counts

    scaled = counts / total * 1e4
    return np.log1p(scaled)
