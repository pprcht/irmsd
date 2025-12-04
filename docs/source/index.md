<div align="center"> 


<!-- <h1>iRMSD</h1> -->
# iRMSD
<h3><em>Molecular Structure Comparison and Ensemble Pruning</em></h3>

[![Latest Version](https://img.shields.io/github/v/release/pprcht/irmsd?color=khaki)](https://github.com/pprcht/irmsd/releases/latest)
[![DOI](https://img.shields.io/badge/DOI-10.1021%2Facs.jcim.4c02143%20-blue)](http://dx.doi.org/10.1021/acs.jcim.4c02143)
[![License: LGPL v3](https://img.shields.io/badge/license-LGPL_v3-coral.svg)](https://www.gnu.org/licenses/lgpl-3.0) 
<!-- ![CI workflow](https://github.com/pprcht/irmsd/actions/workflows/ci.yaml/badge.svg) -->
[![Tests & Coverage](https://github.com/pprcht/irmsd/actions/workflows/tests-and-coverage.yml/badge.svg)](https://github.com/pprcht/irmsd/actions/workflows/tests-and-coverage.yml)
[![codecov](https://codecov.io/gh/pprcht/irmsd/graph/badge.svg?token=Q1O7IRNITG)](https://codecov.io/gh/pprcht/irmsd)

</div>


---

## iRMSD Package 📦

**iRMSD** is a utility toolkit for *molecular structure comparison*, *ensemble pruning*, and *symmetry-aware RMSD analysis*.
It combines a clean **Python API** with an optimized **Fortran backend**, providing fast and robust routines for large conformational ensembles and multiscale computational workflows.

The package offers:

*  *Fast RMSD evaluation* including symmetry handling, canonicalization, and optimal superposition
* *Structure grouping and pruning* based on distance thresholds or iRMSD criteria
* *Flexible Molecule class* with XYZ/extXYZ parsing and ASE interoperability
* *Low- and high-level APIs* that expose direct Fortran wrappers as well as convenient Python abstractions
* *Extendable infrastructure* for future shape metrics and ensemble workflows

iRMSD is designed for researchers working in computational chemistry, conformational sampling, machine learning for molecules, and structural bio/chem-informatics, and integrates seamlessly with tools such as **ASE**, **RDKit**, **CREST/xTB**, or custom multiscale simulation pipelines.

```{toctree} 
    :maxdepth: 1
   :caption: iRMSD

scientific_background.md
installation.md
cli.md
package_reference.rst
```
