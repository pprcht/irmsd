from importlib.resources import files
import ctypes as ct
import sys
from pathlib import Path

def _lib_filename() -> str:
    return {
        "win32": "yourpkg_fortran.dll",
        "cygwin": "yourpkg_fortran.dll",
        "darwin": "libyourpkg_fortran.dylib",
    }.get(sys.platform, "libyourpkg_fortran.so")

def _find_lib() -> str:
    # standard location (works with wheels & redirect editables)
    cand = files("yourpkg") / _lib_filename()
    if cand.exists():
        return str(cand)
    # optional: remove fallback if your redirect setup is solid
    repo_root = Path(__file__).resolve().parents[2]
    for p in repo_root.glob(f"build/**/{_lib_filename()}"):
        if p.is_file():
            return str(p)
    raise FileNotFoundError(f"Cannot locate {_lib_filename()}.")

# Singleton CDLL
LIB = ct.CDLL(_find_lib())

