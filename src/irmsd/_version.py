from __future__ import annotations
from pathlib import Path

try:
    import tomllib  # Python 3.11+
except ModuleNotFoundError:  # Python ≤3.10: need tomli in deps
    import tomli as tomllib


def _get_version() -> str:
    # .../project/src/irmsd/_version.py → project root
    root = Path(__file__).resolve().parents[2]
    pyproject = root / "pyproject.toml"

    data = tomllib.loads(pyproject.read_text(encoding="utf8"))
    try:
        return data["project"]["version"]
    except KeyError as exc:
        raise RuntimeError(
            "Version not found in [project].version in pyproject.toml"
        ) from exc


__version__ = _get_version()
__all__ = ["__version__"]

