import numpy as np
from typing import Tuple
from ..bindings import linalg as _F

def saxpy(a: float, x: np.ndarray, y: np.ndarray) -> np.ndarray:
    x = np.ascontiguousarray(x, dtype=np.float64)
    y = np.ascontiguousarray(y, dtype=np.float64)
    if x.shape != y.shape:
        raise ValueError("x and y must have the same shape")
    _F.saxpy(x.size, float(a), x, y)
    return y

def scal(a: float, x: np.ndarray) -> np.ndarray:
    x = np.ascontiguousarray(x, dtype=np.float64)
    _F.scal(x.size, float(a), x)
    return x

def dot(x: np.ndarray, y: np.ndarray) -> float:
    x = np.ascontiguousarray(x, dtype=np.float64)
    y = np.ascontiguousarray(y, dtype=np.float64)
    if x.shape != y.shape:
        raise ValueError("x and y must have the same shape")
    out = np.empty(1, dtype=np.float64)
    _F.dot(x.size, x, y, out)
    return float(out[0])

# --- Python-only functions live happily here too ---
def cosine_similarity(x: np.ndarray, y: np.ndarray) -> float:
    """Pure Python/NumPy helper built on top of dot/scal etc."""
    num = dot(x, y)
    denom = np.sqrt(dot(x, x) * dot(y, y))
    if denom == 0.0:
        return 0.0
    return num / denom

