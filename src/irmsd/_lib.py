from importlib.resources import files
import ctypes as ct
import sys
from pathlib import Path

def _lib_filenames() -> list[str]:
    if sys.platform in ("win32", "cygwin"):
        # MSVC's name, and MinGW's (which adds a "lib" prefix)
        return ["irmsd_fortran.dll", "libirmsd_fortran.dll"]
    if sys.platform == "darwin":
        return ["libirmsd_fortran.dylib"]
    return ["libirmsd_fortran.so"]

def _find_lib() -> str:
    names = _lib_filenames()
    # standard location (works with wheels & redirect editables)
    for name in names:
        cand = files("irmsd") / name
        if cand.exists():
            return str(cand)
    # optional: remove fallback if the redirect setup is solid
    repo_root = Path(__file__).resolve().parents[2]
    for name in names:
        for p in repo_root.glob(f"build/**/{name}"):
            if p.is_file():
                return str(p)
    raise FileNotFoundError(f"Cannot locate {' or '.join(names)}.")

# Singleton CDLL
LIB = ct.CDLL(_find_lib())

