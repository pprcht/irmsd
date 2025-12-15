from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:  # avoid hard dependency at import time
    from ase import Atoms  # type: ignore

__all__ = [
    "require_ase",
]


# -----------------------------------------------------------------------------
# Dependency guard
# -----------------------------------------------------------------------------


def require_ase() -> ModuleType:
    """Import and return the ASE module, or raise a clear error if it is
    missing."""
    try:
        import ase  # type: ignore[import]
    except ImportError as exc:
        raise RuntimeError(
            "This function requires ASE, but it is not installed. "
            "Install it with `pip install ase`."
        ) from exc
    return ase


def require_rdkit() -> None:
    """Ensure rdkit is importable; raise a helpful ImportError otherwise.

    Use this at the *start* of any function that depends on ASE. We keep this
    separate so that importing `irmsd` does not immediately require ASE.
    """
    try:
        import rdkit  # noqa: F401
    except Exception as e:  # pragma: no cover
        raise ImportError(
            "rdkit is required for this function. Install optional extra: pip install 'irmsd[rdkit]'"
        ) from e
