# Installation

iRMSD is a Python package with a compiled Fortran backend that is loaded via `ctypes`.
It requires Python 3.10 or newer.

## From PyPI

```bash
pip install irmsd              # core package and the irmsd command-line tool
pip install "irmsd[ase]"       # + ASE (reads any format ASE supports, ASE interface)
pip install "irmsd[rdkit]"     # + RDKit interface
```

Prebuilt wheels are available for

| Platform | Wheel tag |
|---|---|
| Linux, x86\_64 and aarch64 | `manylinux` |
| macOS, Apple silicon | `macosx_*_arm64` |
| Windows, x86\_64 | `win_amd64` |

The wheels contain the compiled library, so no compiler is needed.
On macOS, the minimum supported OS version follows the `macos-latest` runner used by GitHub Actions at release time.

## Parallel sorting (OpenMP)

The wheels are built with OpenMP, which parallelizes the iRMSD comparisons inside the ensemble sorting (`irmsd sort`, {py:func}`irmsd.sorter_irmsd_molecule`, {py:func}`irmsd.prune`).
All available cores are used by default; set `OMP_NUM_THREADS` to limit this:

```bash
OMP_NUM_THREADS=4 irmsd sort ensemble.xyz
```

The result does not depend on the number of threads.
Whether the installed library was built with OpenMP can be checked from Python:

```python
>>> import irmsd
>>> irmsd.openmp_info()
{'enabled': True, 'max_threads': 8}
```

## Building from source

On platforms without a wheel (e.g. Windows on ARM), or when installing from a checkout, pip builds the package from source via scikit-build-core and CMake.
This requires

- a C compiler and a Fortran compiler that CMake can find (e.g. `gcc`/`gfortran`, or MinGW `gcc`/`gfortran` on Windows),
- CMake ≥ 3.22 (installed by pip automatically if missing).

```bash
pip install .                              # from a checkout
pip install irmsd --no-binary irmsd        # force a source build from PyPI
```

OpenMP is used if the Fortran compiler supports it, and the build falls back to a serial library otherwise.
Two CMake options control this and can be passed through pip:

| Option | Default | Effect |
|---|---|---|
| `IRMSD_USE_OPENMP` | `ON` | Build with OpenMP if it is found |
| `IRMSD_REQUIRE_OPENMP` | `OFF` | Stop with an error if OpenMP is not found (set for the release wheels) |

```bash
pip install . -C cmake.define.IRMSD_USE_OPENMP=OFF     # serial build
```

## Development install

An editable install picks up changes to the Python sources directly; after changes to the Fortran code, run the install again to rebuild the library:

```bash
pip install -e ".[dev,ase,rdkit]"
pytest
```

The documentation is built with Sphinx from the `docs/` directory:

```bash
pip install ".[docs]"
make -C docs html          # output in docs/_build/html
```
