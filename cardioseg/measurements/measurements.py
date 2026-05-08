import numpy as np
from skimage.measure import regionprops
from skimage.feature import graycomatrix, graycoprops

def compute_cell_measurements(cell_labels, image, nuclei_assignment=None):
    measurements = {}

    for prop in regionprops(cell_labels, intensity_image=image):
        cid = prop.label
        perimeter = float(prop.perimeter)
        area = float(prop.area)

        circularity = 0.0
        if perimeter > 0:
            circularity = (4.0 * np.pi * area) / (perimeter ** 2)
            circularity = float(np.clip(circularity, 0.0, 1.0))

        coords = prop.coords
        region_intensity = image[coords[:, 0], coords[:, 1]]
        standard_intensity = float(np.std(region_intensity))

        measurements[cid] = {
            "area": area,
            "perimeter": perimeter,
            "eccentricity": float(prop.eccentricity),
            "mean_intensity": float(prop.mean_intensity),
            "standard_intensity": standard_intensity,
            "circularity": circularity,
            "solidity": float(prop.solidity),
            "centroid": tuple(map(float, prop.centroid)),  # (y, x)
        }

    return measurements

def compute_export_features(cell_labels, image, include_shape=False, include_texture=False, glcm_levels=32):
    """
    Compute optional export-only features without changing table/info panel defaults.
    Returns {cell_id: {feature_name: value}}.
    """
    features = {}
    if not include_shape and not include_texture:
        return features

    for prop in regionprops(cell_labels, intensity_image=image):
        cid = prop.label
        row = {}

        if include_shape:
            perimeter = float(prop.perimeter)
            area = float(prop.area)
            circularity = 0.0
            if perimeter > 0:
                circularity = (4.0 * np.pi * area) / (perimeter ** 2)
                circularity = float(np.clip(circularity, 0.0, 1.0))

            region_values = prop.intensity_image[prop.image]
            row["standard_intensity"] = float(np.std(region_values))
            row["circularity"] = circularity
            row["solidity"] = float(prop.solidity)

        if include_texture:
            minr, minc, maxr, maxc = prop.bbox
            patch = image[minr:maxr, minc:maxc]
            mask = prop.image
            patch_values = patch[mask]

            if patch_values.size > 1:
                vmin = float(np.min(patch_values))
                vmax = float(np.max(patch_values))
                quantized = np.zeros_like(patch, dtype=np.uint8)
                if vmax > vmin:
                    scaled = (patch - vmin) / (vmax - vmin)
                    quantized = np.clip(np.round(scaled * (glcm_levels - 1)), 0, glcm_levels - 1).astype(np.uint8)
                quantized[~mask] = 0

                glcm = graycomatrix(
                    quantized,
                    distances=[1],
                    angles=[0],
                    levels=glcm_levels,
                    symmetric=True,
                    normed=True,
                )
                row["texture_contrast"] = float(graycoprops(glcm, "contrast")[0, 0])
                row["texture_homogeneity"] = float(graycoprops(glcm, "homogeneity")[0, 0])
            else:
                row["texture_contrast"] = 0.0
                row["texture_homogeneity"] = 1.0

        if row:
            features[cid] = row

    return features
