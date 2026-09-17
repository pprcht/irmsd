from __future__ import annotations

import ctypes as ct

import numpy as np
from numpy.ctypeslib import ndpointer

from .._lib import LIB

# argtypes follow the bind(C) signatures in _fortran/fmods/symmetry_exposed.f90
LIB.get_symmetry_fortran.argtypes = [
    ct.c_int,
    ndpointer(dtype=np.int32, flags="C_CONTIGUOUS"),
    ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),
    ct.c_double,
    ct.c_double,
    ct.c_int,
    ct.c_int,
    ct.c_char_p,
    ct.c_int,
]
LIB.get_symmetry_fortran.restype = None

_SYMLEN = 16


def get_symmetry_fortran_raw(
    natoms: int,
    types: np.ndarray,
    coords_flat: np.ndarray,
    threshold: float,
    primary_threshold: float,
    max_axis_order: int,
    max_opt_cycles: int,
) -> str:
    """Direct call to Fortran ``get_symmetry_fortran``.

    Parameters
    ----------
    natoms : int
        Must match the array lengths.
    types : (natoms,) int32, C-contiguous
        Atomic numbers (or type IDs).
    coords_flat : (3*natoms,) float64, C-contiguous
        [x1, y1, z1, x2, ...] in Angstrom.
    threshold, primary_threshold : float
        Final and atom-pairing symmetry tolerances, in Bohr.
    max_axis_order, max_opt_cycles : int
        Highest axis order searched; optimization cycles per element.

    Returns
    -------
    str
        Schoenflies symbol, e.g. "C2v".

    Raises
    ------
    TypeError
        On wrong dtype or non-contiguous arrays.
    ValueError
        If array lengths disagree with natoms.
    """
    if types.dtype != np.int32 or not types.flags.c_contiguous:
        raise TypeError("types must be int32 and C-contiguous")
    if coords_flat.dtype != np.float64 or not coords_flat.flags.c_contiguous:
        raise TypeError("coords_flat must be float64 and C-contiguous")
    if coords_flat.size != 3 * natoms:
        raise ValueError("coords_flat length must be 3*natoms")
    if types.size != natoms:
        raise ValueError("types length must be natoms")

    buf = ct.create_string_buffer(_SYMLEN)
    LIB.get_symmetry_fortran(
        int(natoms),
        types,
        coords_flat,
        float(threshold),
        float(primary_threshold),
        int(max_axis_order),
        int(max_opt_cycles),
        buf,
        _SYMLEN,
    )
    return buf.value.decode("ascii")


LIB.get_symmetry_elements_fortran.argtypes = [
    ct.c_int,
    ndpointer(dtype=np.int32, flags="C_CONTIGUOUS"),
    ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),
    ct.c_double,
    ct.c_double,
    ct.c_int,
    ct.c_int,
    ct.c_char_p,
    ct.c_int,
    ct.c_int,
    ct.POINTER(ct.c_int),
    ndpointer(dtype=np.int32, flags="C_CONTIGUOUS"),
    ndpointer(dtype=np.int32, flags="C_CONTIGUOUS"),
    ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),
    ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),
    ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),
    ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),
    ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),
    ndpointer(dtype=np.int32, flags="C_CONTIGUOUS"),
]
LIB.get_symmetry_elements_fortran.restype = None

# Ih has 63 elements; larger buffers are only needed for loose thresholds.
_MAXEL = 128


def get_symmetry_elements_fortran_raw(
    natoms: int,
    types: np.ndarray,
    coords_flat: np.ndarray,
    threshold: float,
    primary_threshold: float,
    max_axis_order: int,
    max_opt_cycles: int,
) -> tuple[str, dict[str, np.ndarray]]:
    """Direct call to Fortran ``get_symmetry_elements_fortran``.

    Parameters
    ----------
    natoms, types, coords_flat, threshold, primary_threshold, max_axis_order,
    max_opt_cycles
        As in :func:`get_symmetry_fortran_raw`.

    Returns
    -------
    symbol : str
        Schoenflies symbol.
    elements : dict of ndarray
        Arrays over the nel elements, lengths in Angstrom: ``type`` (nel,)
        int32, 1=mirror, 2=inversion, 3=rotation, 4=improper rotation;
        ``order`` (nel,) int32, 0 for Cinf; ``matrix`` (nel,3,3) and
        ``translation`` (nel,3) with x' = R x + t; ``axis``, ``point``
        (nel,3); ``maxdev`` (nel,); ``permutation`` (nel,natoms) int32,
        0-based image of each atom.
    """
    if types.dtype != np.int32 or not types.flags.c_contiguous:
        raise TypeError("types must be int32 and C-contiguous")
    if coords_flat.dtype != np.float64 or not coords_flat.flags.c_contiguous:
        raise TypeError("coords_flat must be float64 and C-contiguous")
    if coords_flat.size != 3 * natoms:
        raise ValueError("coords_flat length must be 3*natoms")
    if types.size != natoms:
        raise ValueError("types length must be natoms")

    maxel = _MAXEL
    while True:
        buf = ct.create_string_buffer(_SYMLEN)
        nel = ct.c_int(0)
        out = {
            "type": np.zeros(maxel, dtype=np.int32),
            "order": np.zeros(maxel, dtype=np.int32),
            "matrix": np.zeros((maxel, 3, 3), dtype=np.float64),
            "translation": np.zeros((maxel, 3), dtype=np.float64),
            "axis": np.zeros((maxel, 3), dtype=np.float64),
            "point": np.zeros((maxel, 3), dtype=np.float64),
            "maxdev": np.zeros(maxel, dtype=np.float64),
            "permutation": np.zeros((maxel, natoms), dtype=np.int32),
        }
        LIB.get_symmetry_elements_fortran(
            int(natoms),
            types,
            coords_flat,
            float(threshold),
            float(primary_threshold),
            int(max_axis_order),
            int(max_opt_cycles),
            buf,
            _SYMLEN,
            maxel,
            ct.byref(nel),
            out["type"],
            out["order"],
            out["matrix"],
            out["translation"],
            out["axis"],
            out["point"],
            out["maxdev"],
            out["permutation"],
        )
        # Fortran reports the full count but writes only min(nel, maxel)
        if nel.value <= maxel:
            break
        maxel = nel.value

    n = nel.value
    return buf.value.decode("ascii"), {k: v[:n] for k, v in out.items()}
