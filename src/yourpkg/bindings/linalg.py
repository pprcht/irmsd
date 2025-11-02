from numpy.ctypeslib import ndpointer
import ctypes as ct
import numpy as np
from .._lib import LIB

# Fortran: subroutine saxpy(n, a, x, y) bind(C)
LIB.saxpy.argtypes = [
    ct.c_int,
    ct.c_double,
    ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),
    ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),
]
LIB.saxpy.restype = None

def saxpy(n: int, a: float, x: np.ndarray, y: np.ndarray) -> None:
    LIB.saxpy(int(n), float(a), x, y)

# Example additional exported functions:
# subroutine scal(n, a, x) bind(C)
LIB.scal.argtypes = [
    ct.c_int,
    ct.c_double,
    ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),
]
LIB.scal.restype = None

def scal(n: int, a: float, x: np.ndarray) -> None:
    LIB.scal(int(n), float(a), x)

# subroutine dot(n, x, y, out) bind(C)  ! out is scalar length-1
LIB.dot.argtypes = [
    ct.c_int,
    ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),
    ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),
    ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),
]
LIB.dot.restype = None

def dot(n: int, x: np.ndarray, y: np.ndarray, out: np.ndarray) -> None:
    LIB.dot(int(n), x, y, out)

