#spatial expression engine

from ast import expr

import numpy as np
from skimage.draw import polygon
from .annotation import (
    MARKER_DICT,
    build_marker_index,
    build_reference_profiles,
    score_markers
)

import os
import scanpy as sc
import pandas as pd
from matplotlib.path import Path

from scipy.spatial import KDTree
class SpatialExpressionEngine:
    def __init__(self, h5_path):
        print("Initializing Spatial Expression Engine...")
        
        # 1. Store path and load the expression matrix
        self.h5_path = h5_path
        self.adata = sc.read_10x_h5(h5_path)
        self.adata.var_names_make_unique()

        self.adata.var_names = [g.strip().capitalize() for g in self.adata.var_names]
        
        self.is_mito = self.adata.var_names.str.startswith(('Mt-', 'MT-'))
        print(f"Detected {self.is_mito.sum()} mitochondrial genes to be excluded from scoring.")

        self.expr_matrix = self.adata.X

        # 2. Locate spatial directory (assuming standard SpaceRanger structure)
        # If h5 is in binned_outputs, we need to go up to get to 'spatial'
        base_dir = os.path.dirname(h5_path)
        spatial_dir = os.path.join(base_dir, "spatial")
        
        # Fallback for nested binned folders
        if not os.path.exists(spatial_dir):
            spatial_dir = os.path.join(os.path.dirname(os.path.dirname(base_dir)), "spatial")

        # 3. Load tissue positions
        parquet_path = os.path.join(spatial_dir, "tissue_positions.parquet")
        csv_path = os.path.join(spatial_dir, "tissue_positions.csv")

        if os.path.exists(parquet_path):
            tissue_positions = pd.read_parquet(parquet_path)
        elif os.path.exists(csv_path):
            tissue_positions = pd.read_csv(csv_path, header=None)
            tissue_positions.columns = [
                "barcode", "in_tissue", "array_row", "array_col",
                "pxl_row_in_fullres", "pxl_col_in_fullres"
            ]
        else:
            raise ValueError(f"No tissue positions found in {spatial_dir}")

        # 4. Align coordinates to AnnData barcodes
        if "barcode" in tissue_positions.columns:
            tissue_positions = tissue_positions.set_index("barcode")
        
        tissue_positions = tissue_positions.loc[self.adata.obs_names]
        
        # Note: SpaceRanger uses (col, row) for (x, y)
        self.bin_coords = tissue_positions[["pxl_col_in_fullres", "pxl_row_in_fullres"]].values
        self.adata.obsm["spatial"] = self.bin_coords

        print("Building spatial index (KDTree)...")
        self.kdtree = KDTree(self.bin_coords)

        # 5. Build marker index (Precomputes which columns match which cell types)
        self.marker_indices = build_marker_index(MARKER_DICT, self.adata.var_names)
        self.reference_profiles = build_reference_profiles(
            self.marker_indices,
            len(self.adata.var_names)
        )
        
        print(f"Spatial engine ready. Loaded {len(self.adata)} spots.")

    def _apply_mito_filter(self, expr_vector):
        """Helper to zero out mito genes in an expression vector."""
        if expr_vector is not None:
            expr_vector[self.is_mito] = 0
        return expr_vector
    
    def _polygon_to_spot_indices(self, polygon):

        centroid = np.mean(polygon, axis=0)
        bbox = polygon.max(axis=0) - polygon.min(axis=0)
        radius = np.linalg.norm(bbox) / 2

        candidate = self.kdtree.query_ball_point(centroid, r=radius + 10)

        if not candidate:
            return []

        coords = self.bin_coords[candidate]

        path = Path(polygon)

        mask = path.contains_points(coords)

        return np.array(candidate)[mask]

    def _get_weighted_expression(self, polygon, k=20, max_dist=200):
        """
        Finds the k-nearest spots to the nucleus centroid and returns 
        a distance-weighted sum of their expression.
        """
        # 1. Get the center of the nucleus
        centroid = np.mean(polygon, axis=0).reshape(1, -1)
        
        # 2. Query the KDTree for the k nearest Visium spots
        dists, indices = self.kdtree.query(centroid, k=k)
        dists = dists.flatten()
        indices = indices.flatten()

        # 3. Filter out spots that are too far away (e.g., outside the tissue)
        valid = dists <= max_dist
        if not np.any(valid):
            return None
        
        v_indices = indices[valid]
        v_dists = dists[valid]

        # 4. Calculate weights (Inverse Distance Weighting)
        # Add a tiny epsilon to avoid division by zero if centroid matches a spot exactly
        weights = 1.0 / np.maximum(v_dists, 1e-6)
        weights /= weights.sum()

        # 5. Extract and weight the expression
        # Use .X for the raw/normalized counts
        X_subset = self.expr_matrix[v_indices]
        
        if hasattr(X_subset, "toarray"):
            X_subset = X_subset.toarray()
            
        weighted_expr = np.sum(X_subset * weights[:, np.newaxis], axis=0)

        return self._apply_mito_filter(weighted_expr)

    def _extract_polygon_expression(self, polygon):

        indices = self._polygon_to_spot_indices(polygon)

        if len(indices) == 0:
            return None

        X = self.expr_matrix[indices]

        # sparse → dense only for the tiny subset
        if hasattr(X, "toarray"):
            X = X.toarray()

        raw_sum = X.sum(axis=0)
    
        return self._apply_mito_filter(raw_sum)

    def annotate_polygon(self, polygon, centroid=None, annotate_cell_type = True, annotate_genes = True):
        if centroid is None:
            centroid = np.mean(polygon, axis=0)

        expr = self._extract_polygon_expression(polygon)

        if expr is None:
            expr = self._get_weighted_expression(polygon, k=20)

        if expr is None:
            return {"cell_type": "Unknown", "confidence": 0.0, "score": 0.0}

        # --- NEW: MITO FILTERING ---
        # Create a mask for genes that start with 'mt-' or 'MT-'
        # This assumes self.adata.var_names is accessible here
        is_mito = self.adata.var_names.str.startswith(('mt-', 'MT-'))
        
        # Zero out mitochondrial expression so they don't count towards the sum
        expr[is_mito] = 0 

        total_sum = expr.sum()
        raw_expr = expr.copy()
        total_sum = raw_expr.sum()

        if total_sum == 0:
            return {"cell_type": "Unknown", "confidence": 0.0, "score": 0.0}


        expr_norm = raw_expr / (total_sum + 1e-9)

        if annotate_cell_type:
            cell_type, confidence, score = score_markers(
                expr_norm,
                self.reference_profiles
            )
        else:
            cell_type, confidence, score = "Unknown", 0.0, 0.0
            
        gene_dict = {}
        gene_dict_raw = {}
        if annotate_genes:
            top_n = 50
            top_idx = np.argsort(raw_expr)[-top_n:]

            gene_dict = {
                self.adata.var_names[i]: float(expr_norm[i])
                for i in top_idx
            }

            gene_dict_raw = {
                self.adata.var_names[i]: float(raw_expr[i])
                for i in top_idx
            }

        return {
            "cell_type": cell_type,
            "confidence": float(confidence),
            "score": float(score),
            "genes": gene_dict,
            "genes_raw": gene_dict_raw
        }

def diagnostic_polygon_check(engine, polygon, polygon_id=None, filepath=None, print_to_console=True):
    import os
    import numpy as np
    from matplotlib.path import Path
    
    # --- 1. Get Data ---
    path = Path(polygon)
    coords = engine.adata.obsm["spatial"]
    inside_mask = path.contains_points(coords)
    num_inside = inside_mask.sum()
    
    top_genes = "N/A"
    total_counts = 0
    
    if num_inside > 0:
        # Extract and sum expression
        X_slice = engine.adata.X[inside_mask]
        if hasattr(X_slice, "toarray"): X_slice = X_slice.toarray()
        
        gene_sums = np.array(X_slice.sum(axis=0)).flatten()
        total_counts = gene_sums.sum()
        
        # Find Top 3 expressed genes in this specific cell
        top_idx = np.argsort(gene_sums)[-3:][::-1]
        top_genes = "|".join([f"{engine.adata.var_names[i]}:{int(gene_sums[i])}" for i in top_idx])

    # --- 2. Reason Logic ---
    status = "OK" if num_inside > 0 else "EMPTY"
    if num_inside > 0 and total_counts < 5:
        status = "LOW_SIGNAL"

    # --- 3. File Output ---
    if filepath:
        write_header = not os.path.exists(filepath)
        with open(filepath, "a") as f:
            if write_header:
                f.write("id,spots,total_counts,top_genes,status\n")
            f.write(f"{polygon_id},{num_inside},{total_counts:.1f},{top_genes},{status}\n")

    if print_to_console:
        print(f"[{polygon_id}] {status} | Counts: {total_counts:.0f} | Top: {top_genes}")