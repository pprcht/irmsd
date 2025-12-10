from __future__ import annotations

"""
Auto-generated bindings for sorter wrapping
"""
import ctypes as ct
from typing import Optional

import numpy as np
from numpy.ctypeslib import ndpointer

from .._lib import LIB

_sorter_exposed_xyz_fortran: ct._CFuncPtr | None = None
_delta_irmsd_list_fortran: ct._CFuncPtr | None = None
_cregen_exposed_fortran: ct._CFuncPtr | None = None 

def _get_sorter_exposed_xyz_fortran() -> ct._CFuncPtr:
    global _sorter_exposed_xyz_fortran
    if LIB is None:
        raise RuntimeError("Library handle not set; call set_library(...) first.")
    if _sorter_exposed_xyz_fortran is None:
        f = LIB.sorter_exposed_xyz_fortran
        f.argtypes = [
            ct.c_int,
            ct.c_int,
            ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),
            ndpointer(dtype=np.int32, flags="C_CONTIGUOUS"),
            ndpointer(dtype=np.int32, flags="C_CONTIGUOUS"),
            ct.c_double,
            ct.c_int,
            ct.c_bool,
            ct.c_int,
            ct.c_double,
            ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),
        ]
        f.restype = None
        _get_sorter_exposed_xyz_fortran = f
    return _get_sorter_exposed_xyz_fortran


def sorter_exposed_xyz_fortran_raw(
    nat,
    nall,
    xyzall,
    atall,
    groups,
    rthresh,
    iinversion,
    allcanon,
    printlvl,
    ethr,
    energies,
) -> None:
    """Low-level wrapper for C symbol 'sorter_exposed_xyz_fortran'."""
    f = _get_sorter_exposed_xyz_fortran()
    f(
        nat,
        nall,
        xyzall,
        atall,
        groups,
        rthresh,
        iinversion,
        allcanon,
        printlvl,
        ethr,
        energies,
    )


def _get_delta_irmsd_list_fortran() -> ct._CFuncPtr:
    global _delta_irmsd_list_fortran
    if LIB is None:
        raise RuntimeError("Library handle not set; call set_library(...) first.")
    if _delta_irmsd_list_fortran is None:
        f = LIB.delta_irmsd_list_fortran
        f.argtypes = [
            ct.c_int,
            ct.c_int,
            ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),
            ndpointer(dtype=np.int32, flags="C_CONTIGUOUS"),
            ct.c_int,
            ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),
            ct.c_bool,
            ct.c_int,
        ]
        f.restype = None
        _get_delta_irmsd_list_fortran = f
    return _get_delta_irmsd_list_fortran


def delta_irmsd_list_fortran_raw(
    nat, nall, xyzall, atall, iinversion, delta, allcanon, printlvl
) -> None:
    """Low-level wrapper for C symbol 'delta_irmsd_list_fortran'."""
    f = _get_delta_irmsd_list_fortran()
    f(nat, nall, xyzall, atall, iinversion, delta, allcanon, printlvl)


#    subroutine cregen_exposed_fortran( &
#  &                     nat,nall,xyzall_ptr,atall_ptr, &
#  &                     groups_ptr,rthresh,ethr,bthr,printlvl, &
#  &                     energies_ptr &
#  &                   ) bind(C,name="cregen_exposed_fortran")


def _get_cregen_exposed_fortran() -> ct._CFuncPtr:
    global _cregen_exposed_fortran
    if LIB is None:
        raise RuntimeError("Library handle not set; call set_library(...) first.")
    if _cregen_exposed_fortran is None:
        f = LIB.cregen_exposed_fortran
        f.argtypes = [
            ct.c_int,
            ct.c_int,
            ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),
            ndpointer(dtype=np.int32, flags="C_CONTIGUOUS"),
            ndpointer(dtype=np.int32, flags="C_CONTIGUOUS"),
            ct.c_double,
            ct.c_double,
            ct.c_double,
            ct.c_int,
            ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),
        ]
        f.restype = None
        _get_cregen_exposed_fortran = f
    return _get_cregen_exposed_fortran


def cregen_fortran_raw(
    nat, nall, xyzall, atall, groups, rthr, ethr, bthr, printlvl, energies
) -> None:
    """Low-level wrapper for C symbol 'cregen_exposed_fortran'."""
    f = _get_cregen_exposed_fortran()
    f(nat, nall, xyzall, atall, groups, rthr, ethr, bthr, printlvl, energies)
