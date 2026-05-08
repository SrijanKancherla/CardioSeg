from dataclasses import dataclass
from typing import Tuple, List

@dataclass
class Cell:
    id: int
    label: int
    centroid: Tuple[float, float]
    area: float
    perimeter: float
    eccentricity: float
    mean_intensity: float
    nuclei_ids: List[int]
