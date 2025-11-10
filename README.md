# iRMSD

<div align="center">

[![Latest Version](https://img.shields.io/github/v/release/pprcht/irmsd?color=khaki)](https://github.com/pprcht/irmsd/releases/latest)
[![DOI](https://img.shields.io/badge/DOI-10.1021%2Facs.jcim.4c02143%20-blue)](http://dx.doi.org/10.1021/acs.jcim.4c02143)
[![License: LGPL v3](https://img.shields.io/badge/license-LGPL_v3-coral.svg)](https://www.gnu.org/licenses/lgpl-3.0) 
<!-- ![CI workflow](https://github.com/pprcht/irmsd/actions/workflows/ci.yaml/badge.svg) -->

</div>


A package for the calculation of the permutation-invariant root-mean-square-deviation (RMSD) of atomic positions.

## Quick start

```bash
python -m pip install -U pip
pip install -e .
pytest -q
```

## Code exposed from Fortran 

- [x] CN calculation (D4/cov CN)
- [x] axis0 to get rotational constants
- [x] "canonical" atom identifiers
- [x] Quaternion RMSD calculator
- [ ] ~~Hungarian/LSAP calculator~~
- [ ] two/multiple structure iRMSD procedure (full)
- [ ] Ensemble sorter

### Other TODOs

- [ ] Clean up & unify printouts
- [ ] Think about functions for the cli
- [ ] Write tests (pytest)
- [ ] Write documentation
- [ ] PyPi and conda(-forge) publishing

