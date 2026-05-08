#loading AnnData (.h5ad) + coords

import scanpy as sc
import pandas as pd
import numpy as np

def load_spatial_dataset(h5_path):
    """
    Loads Visium HD filtered_feature_bc_matrix.h5
    and attaches full-resolution spatial coordinates.
    """
    adata = sc.read_10x_h5(h5_path)
    adata.var_names_make_unique()
    return adata
