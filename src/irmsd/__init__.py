from .api.linalg import saxpy
from .api.xyz_bridge import xyz_to_fortran, xyz_to_fortran_pair

# Try to expose ase_to_fortran if ASE is present; otherwise export a stub that errors nicely.
try:
    from .utils.ase_io import ase_to_fortran, ase_to_fortran_pair
except Exception:
    def ase_to_fortran(*args, **kwargs):  # type: ignore
        raise ImportError(
            "ASE-based helper not available. Install optional extra: pip install 'irmsd[ase]'"
        )

__all__ = [
    "saxpy",
    "xyz_to_fortran",
    "xyz_to_fortran_pair",
    "ase_to_fortran",
    "ase_to_fortran_pair",
]

