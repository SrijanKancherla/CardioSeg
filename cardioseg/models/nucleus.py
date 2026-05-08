from dataclasses import dataclass
from typing import Tuple, Optional

@dataclass
class Nucleus:
    id: int
    label: int
    centroid: Tuple[float, float]
    area: float
    parent_cell: Optional[int]
