<div align="center"> 


<h1>iRMSD</h1>
<h2><em>Molecular Structure Comparison and Ensemble Pruning</em></h2>

[![Latest Version](https://img.shields.io/github/v/release/pprcht/irmsd?color=khaki)](https://github.com/pprcht/irmsd/releases/latest)
[![DOI](https://img.shields.io/badge/DOI-10.1021%2Facs.jcim.4c02143%20-blue)](http://dx.doi.org/10.1021/acs.jcim.4c02143)
[![License: LGPL v3](https://img.shields.io/badge/license-LGPL_v3-coral.svg)](https://www.gnu.org/licenses/lgpl-3.0) 
<!-- ![CI workflow](https://github.com/pprcht/irmsd/actions/workflows/ci.yaml/badge.svg) -->

</div>


---

## iRMSD Package 📦

**iRMSD** is a high-performance toolkit for **molecular structure comparison**, **ensemble pruning**, and **symmetry-aware RMSD analysis**.
It combines a clean **Python API** with a optimized **Fortran backend**, providing fast and robust routines for large conformational ensembles and multiscale computational workflows.

The package offers:

*  **Fast RMSD evaluation** using symmetry handling, canonicalization, and optimal superposition
* **Structure grouping and pruning** based on distance thresholds or iRMSD criteria
* **Flexible Molecule class** with XYZ/extXYZ parsing and ASE interoperability
* **Low- and high-level APIs** that expose direct Fortran wrappers as well as convenient Python abstractions
* **Extendable infrastructure** for future shape metrics and ensemble workflows

iRMSD is designed for researchers working in computational chemistry, conformational sampling, machine learning for molecules, and structural bio/chem-informatics, and integrates seamlessly with tools such as **ASE**, **RDKit**, **CREST/xTB**, or custom multiscale simulation pipelines.


## Installation 🔧

iRMSD is available via both **PyPI** and **conda-forge**.

**PyPI:**

```bash
pip install irmsd
```

**Conda:**

```bash
conda install -c conda-forge irmsd
```



# Scientific Background 🔬 
[![DOI](https://img.shields.io/badge/DOI-10.1021%2Facs.jcim.4c02143%20-blue)](http://dx.doi.org/10.1021/acs.jcim.4c02143) 

Structural comparison is at the heart of conformational analysis, docking, trajectory processing, and molecular shape workflows. However, the *classical* Cartesian RMSD fails whenever two structures differ by atom ordering, local symmetry, or rotameric permutations. In such cases, the RMSD becomes artificially large even though the molecular properties (e.g. IR/Raman spectra, NMR shifts, energies) are identical. 

The **permutation invariant RMSD (iRMSD)** implemented in this package solves this problem by:

* Assigning **canonical atom identities** independent of input atom order
* Performing **symmetry-aware alignment** using a divide-and-conquer approach
* Solving the **linear sum assignment problem (LSAP, Hungarian algorithm)** efficiently in Fortran
* Handling **false enantiomers** when appropriate
* Using **cached memory + single precision LSAP kernel** for high speed

Together, these yield a robust, fast, and scalable measure of structural similarity suitable for large ensembles (thousands of conformers), proteins, and noncovalent clusters.

Herein, the RMSD with optimal alignment and permutation is defined as:

```math
\mathrm{iRMSD}(\mathbf{X}, \mathbf{Y}) = \min_{\mathbf{P},,\mathbf{U}} \sqrt{ \frac{1}{N} \sum_{i=1}^{N} \left\lVert \mathbf{X}_i - (\mathbf{P},\mathbf{U},\mathbf{Y})_i \right\rVert^2 }
```

where
* (**U**): rotation matrix
* (**P**): permutation matrix representing atom–atom assignment
* The minimization ensures best superposition *and* best atom correspondence

Details on the iRMSD method are extensively discussed in [***J. Chem. Inf. Model.* **2025**, *65* (9), 4501–4511**](https://pubs.acs.org/doi/abs/10.1021/acs.jcim.4c02143) (publically accessible preprint PDF → [**here**](https://doi.org/10.26434/chemrxiv-2024-qmcz4)).


##  Examples  📊

### Example 1 - Symmetric Rotamers (Pentane C<sub>5</sub>H<sub>12</sub>)

```
<!-- Figure Placeholder A:
     Rotamer examples + RMSD heatmaps (conventional vs. iRMSD) -->
<p align="center"><i>[ insert pentane rotamer images and heatmaps here ]</i></p>
```

Most pentane conformers have *four* rotamers that differ only by permutation or inversion of terminal methyl groups.
As shown in Figure 3 of the paper (page 3) :

* Conventional RMSD between these rotamers ranges from 1.2-1.8 Å
* iRMSD correctly identifies all as identical (≈ 0.0 Å) by accounting for permutations and false enantiomers
* Property calculations (IR spectra, NMR shifts) also confirm they are physically equivalent

This simple example illustrates a key problem with classical (quaternion) RMSD-based conformer comparison and the necessity of addressing *both* the alignment *and* permutation problems for chemical workflows.


### Example 2  Validation on Randomized Atom Order Structures

```
<!-- Figure Placeholder B:
     Polar plots of RMSD for randomized copies of pentane/TPPO/taxol/BPTI -->
<p align="center"><i>[ insert validation polar plots for multiple molecules ]</i></p>
```

A robust alignment algorithm must correctly classify structures as identical even if:

* randomly rotated
* atom order is randomly permuted
* (optional) mirrored ("*false enantiomers*")

Figures 7a–d (pages 6–7) show that iRMSD successfully returns ~0 Å for **every pair** in 100 randomized input coordinates of:

* pentane (12 atoms)
* TPPO (46 atoms)
* taxol (113 atoms)
* BPTI (892 atoms)



### Example 3 - Noncovalent Clusters (LJ<sub>75</sub> and (H<sub>2</sub>O)<sub>21</sub>)

```
<!-- Figure Placeholder C:
     RMSD validation for LJ75 and (H2O)21 clusters -->
<p align="center"><i>[ insert LJ75 and water cluster iRMSD plots ]</i></p>
```

Noncovalent clusters break most RMSD algorithms because:

* the molecular graph is *disconnected*
* atom types are all identical (e.g., LJ<sub>75</sub>)
* graph-isomorphism methods become impossible

iRMSD handles these correctly because:

* canonical atom identity assignment does **not** require a connected graph and automatically falls back to the atom types
* LSAP efficiently handles the remaining permutation

For LJ<sub>75</sub>, the full 75×75 LSAP is solved successfully.



### Example 4 - Conformer-Rotamer Ensemble (CRE) Pruning

```
<!-- Figure Placeholder D:
     Parity plot of quaternion RMSD vs iRMSD for ensemble pruning -->
<p align="center"><i>[ insert parity plots for pentane and taxol CREs ]</i></p>
```

iRMSD excels in distinguishing on a single threshold parameter (`RTHR`):

* rotamers (same conformer → iRMSD ≈ 0)
* different conformers (larger iRMSD values)

This is crucial for automated CRE pruning and is an extension to conventional (quaternion) RMSD pruning, e.g. as used in CREST.


## When and How To Use iRMSD?

Use iRMSD whenever you wish to:

* compare conformers with local symmetry
* prune redundant structures in CRE/CREGEN-style workflows
* merge ensembles from different levels of theory (energy-threshold independent)
* process MD trajectories where atom ordering is stable but symmetry remains
* align structures from different toolchains (RDKit ↔ xTB ↔ ORCA ↔ OpenBabel)
* classify rotamers vs conformers with physical correctness


### Python CLI Usage
TODO

### Python Script Usage
TODO



### License
TODO

---

# OTHER TODOs
- [ ] docstrings and actual docs (GH pages?)
- [ ] conda-forge package
- [ ] ci.yml
- [ ] codecov
- [ ] (implementation) Parallelization via OpenMP
- [ ] (implementation) Optional pass of inter-conformer energy threshold (`ethr`)
- [ ] (implementation) Pre-alignment via quaternion RMSD of unique canonical atoms instead of aligning via rotational constants
- [ ] (implementation) "classical" CREGEN pruner based on energy + rot.const. + quaternion RMSD
