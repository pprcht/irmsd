<div align="center"> 


<h1>iRMSD</h1>
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


## Installation 🔧

iRMSD is available via **PyPI**. <!-- both **PyPI** and **conda-forge**. -->

**PyPI:**

```bash
pip install irmsd
```

Prebuilt wheels are available for Linux (`manylinux`) and macOS (`arm64`), so installation is typically straightforward. On macOS, wheel compatibility follows the version of the `macos-latest` runner
used by GitHub Actions. 

We currently do not provide prebuilt wheels for Windows. On Windows, `pip install irmsd` will build the package *from source*. This requires: **1**) Python 3.10 or newer, **2**) A C/C++ compiler (e.g. Visual Studio Build Tools), **3**) A Fortran compiler (e.g. Intel oneAPI Fortran, or MinGW gfortran configured for CMake). If you don’t have a Fortran toolchain available, the installation will likely fail. As an alternative, you can use WSL (Linux subsystem) and install the Linux wheels there.

<!--
**Conda:**

```bash
conda install -c conda-forge irmsd
```
-->

For basic usage instructions, both via the CLI and in-code, [**see below**](#when-and-how-to-use-irmsd).


<br>
<br>

# Scientific Background 🔬 
[![DOI](https://img.shields.io/badge/DOI-10.1021%2Facs.jcim.4c02143%20-blue)](http://dx.doi.org/10.1021/acs.jcim.4c02143) 

Structural comparison is at the heart of conformational analysis, docking, trajectory processing, and molecular shape workflows. However, the *classical* Cartesian RMSD fails whenever two structures differ by atom ordering, local symmetry, or rotameric permutations. In such cases, the RMSD becomes artificially large even though the molecular properties (e.g. IR/Raman spectra, NMR shifts, energies) are identical. 

The **permutation invariant RMSD (iRMSD)** implemented in this package solves this problem by:

- Assigning **canonical atom identities** independent of input atom order
- Performing **symmetry-aware alignment** using a divide-and-conquer approach
- Solving the **linear sum assignment problem (LSAP, Hungarian algorithm)** efficiently in Fortran
- Handling **false enantiomers** when appropriate
- Using **cached memory + single precision LSAP kernel** for high speed

Together, these yield a robust, fast, and scalable measure of structural similarity suitable for large ensembles (thousands of conformers), proteins, and noncovalent clusters.

Herein, the RMSD with optimal alignment and permutation is defined as:

```math
\mathrm{iRMSD}(\mathbf{X}, \mathbf{Y}) = \min_{\mathbf{P},\mathbf{U}} \sqrt{ \frac{1}{N} \sum_{i=1}^{N} \left\lVert \mathbf{X}_i - (\mathbf{P}\mathbf{U}\mathbf{Y})_i \right\rVert^2 }
```

where
* **X**/**Y** : Cartesian coordinates of the two molecules to be compared
* **U**: rotation matrix
* **P**: permutation matrix representing atom–atom assignment
* The minimization ensures best superposition *and* best atom correspondence

Details on the iRMSD method are extensively discussed in [***J. Chem. Inf. Model.* **2025**, *65* (9), 4501–4511**](https://pubs.acs.org/doi/abs/10.1021/acs.jcim.4c02143) (publically accessible preprint PDF → [**here**](https://doi.org/10.26434/chemrxiv-2024-qmcz4)).

<br>

A simple application example is seen below. Here, we have two copies of the same conformer of a drug molecule. Both have entirely random atom order, and are randomly rotated, as seen in the side-by-side coordinate blocks. iRMSD is capable of aligning the structures (i.e., determine **U**) and restoring the correct atom order (i.e., determine **P**).

<table>
  <tr>
  <td >

<p align="center">
  <img src="https://github.com/pprcht/irmsd/raw/main/assets/images/fluoxetine1.jpg" width="75%">
</p>

`struc1.xyz`:

   ```
40

C  -0.0198   0.2158   0.5308
F  -5.0246  -0.2464   0.5593
H  -2.0207   0.0761   3.2888
H   2.0221   5.0303   0.4892
C  -0.4032   0.2003   1.8748
H   1.4542   5.1890  -1.1937
C   2.1773   2.6686  -0.7626
F  -4.6775   0.8294   2.4172
F  -4.3667  -1.3218   2.3289
N   3.0555   3.8290  -0.9233
C   2.3371   5.0467  -0.5601
H   3.7051   1.1947  -0.3434
O   1.3305   0.2076   0.3364
C  -2.3576   0.1386  -0.0950
C   3.3367  -3.7475  -2.0944
C  -2.7426   0.0481   1.2484
C  -4.1855  -0.1659   1.6292
C   2.4651  -1.1750  -1.3319
H   1.7906   2.6512   0.2660
C   2.8963   1.3473  -1.0679
H   1.3203   2.7752  -1.4405
C   1.8806   0.1908  -0.9905
H  -0.7865   0.2662  -1.5101
H   3.3527   1.3925  -2.0651
H   1.1282   0.3962  -1.7541
H   2.8529  -0.6554  -3.3999
H   2.9547  -4.3108  -0.0530
H   3.8626   3.7250  -0.3037
H   0.3581   0.2294   2.6522
C  -1.7511   0.1165   2.2345
H   2.9886   5.9146  -0.7002
H   2.1980  -2.0723   0.6278
H   3.6193  -2.9009  -4.0530
H  -3.1084   0.1174  -0.8845
C   2.5162  -2.2260  -0.4015
C  -1.0075   0.2216  -0.4508
H   3.6502  -4.7430  -2.3936
C   2.9433  -3.5025  -0.7811
C   2.8862  -1.4373  -2.6475
C   3.3158  -2.7117  -3.0256
   ```

  </td><td>

<p align="center">
  <img src="https://github.com/pprcht/irmsd/raw/main/assets/images/fluoxetine2.jpg" width="75%"> 
</p>

`struc2.xyz`:

   ```
40

C  -0.2470  -1.0143  -0.4213
H   0.0828   4.2564  -4.4517
C  -2.1193   2.6178   0.0144
C   1.2089  -1.8204   1.8160
C   0.0619  -0.0804   0.5643
C   0.9632  -2.7436   0.7925
H  -5.9431   1.8117   2.4085
H  -3.2933   0.8099   0.1366
H   2.4555   4.5828  -3.8102
H  -0.0328  -3.0384  -1.1072
H  -3.9856   3.0777   1.8789
H  -2.4640   1.2332   1.6473
C  -0.4255   3.2037  -2.6383
H   1.7581   2.5331  -0.1112
C   1.3640   2.9121  -1.0521
C   1.5502  -4.1306   0.8575
C  -5.0208   1.3357   2.0622
H  -4.4695   0.9920   2.9447
F   1.1253  -4.8340   1.9430
C   2.2338   3.5914  -1.9111
H  -1.7744   3.4140   0.6847
C   0.4427   3.8850  -3.4947
C   0.7599  -0.5019   1.6995
F   1.2525  -4.8990  -0.2270
H  -1.3360   1.0931  -1.2358
H   1.7763  -2.1151   2.6974
H   0.9785   0.2106   2.4926
H  -5.3029   0.4654   1.4587
H  -0.8129  -0.7679  -1.3112
H  -2.6932   3.0927  -0.7916
C  -0.9181   1.8428  -0.5615
F   2.9112  -4.1100   0.9336
H  -1.4519   3.0470  -2.9556
H   3.2752   3.7325  -1.6299
C   1.7755   4.0721  -3.1352
C   0.0230   2.6924  -1.4078
N  -4.2396   2.2887   1.2797
O  -0.2303   1.2522   0.5528
C  -3.0176   1.6462   0.7923
C   0.1996  -2.3349  -0.3081
   ```       

  </td>
  </tr>
</table>

<p align="center">
  <img src="https://github.com/pprcht/irmsd/raw/main/assets/images/Screenshot_20251212_105945.jpg" width="100%"> 
</p>



<br>
<br>

##  Examples  📊

### Example 1 - Symmetric Rotamers (Pentane C<sub>5</sub>H<sub>12</sub>)

<p align="center">
  <img src="https://github.com/pprcht/irmsd/raw/main/assets/images/example1.jpg" width="100%">
</p>

Most pentane conformers have *four* rotamers that differ only by permutation or inversion of terminal methyl groups.
As shown in Figure 3 of the paper (page 3) :

* Conventional RMSD between these rotamers ranges from 1.2-1.8 Å
* iRMSD correctly identifies all as identical (≈ 0.0 Å) by accounting for permutations and false enantiomers
* Property calculations (IR spectra, NMR shifts) also confirm they are physically equivalent

This simple example illustrates a key problem with classical (quaternion) RMSD-based conformer comparison and the necessity of addressing *both* the alignment *and* permutation problems for chemical workflows.

<br>

### Example 2 - Validation on Randomized Atom Order Structures

<p align="center">
  <img src="https://github.com/pprcht/irmsd/raw/main/assets/images/example2.jpg" width="100%"> 
</p>

A robust permutation-handling alignment algorithm must correctly classify structures as identical even if:

* randomly rotated
* atom order is randomly permuted
* (optional) mirrored ("*false enantiomers*")

Figures 7a–d (pages 6–7) show that iRMSD successfully returns ~0 Å for **every pair** in 100 randomized input coordinates of: pentane (12 atoms), TPPO (46 atoms), taxol (113 atoms), BPTI (892 atoms).


<br>

### Example 3 - Noncovalent Clusters (LJ<sub>75</sub> and (H<sub>2</sub>O)<sub>21</sub>)


<p align="center">
  <img src="https://github.com/pprcht/irmsd/raw/main/assets/images/example3.jpg" width="100%">
</p>

Noncovalent clusters break most RMSD algorithms because:

* the molecular graph is *disconnected*
* atom types are all identical (e.g., LJ<sub>75</sub>)
* graph-isomorphism methods become impossible

iRMSD handles these correctly because:

* canonical atom identity assignment does **not** require a connected graph and automatically falls back to the atom types
* LSAP efficiently handles the remaining permutation

For LJ<sub>75</sub>, the full 75×75 LSAP is solved successfully.


<br>

### Example 4 - Conformer-Rotamer Ensemble (CRE) Pruning

<p align="center">
  <img src="https://github.com/pprcht/irmsd/raw/main/assets/images/example4.jpg" width="100%"> 
</p>

iRMSD excels in distinguishing on a single threshold parameter (`RTHR`):

* rotamers (same conformer → iRMSD ≈ 0 Å)
* different conformers (larger iRMSD values)

This is crucial for automated CRE pruning and is an extension to conventional (quaternion) RMSD pruning, e.g. as used in CREST.<br>
The default `RTHR` threshold in `irmsd` to distinguish to structures as conformers is **0.125 Å**, which was adapted from CREST's CREGEN procedure. Additional thresholds, e.g. for the inter-conformer energy difference (`ETHR`) or rotational constants (`BTHR`) are *not* required, but can be used to achieve more efficient pre-sorting.

<br>
<br>

# When and How To Use iRMSD?

Use iRMSD whenever you wish to:

* compare conformers with *local* symmetry
* prune redundant structures in CRE/CREGEN-style workflows
* merge ensembles from different levels of theory (energy-threshold independent)
* process MD trajectories where atom ordering is stable but symmetry remains
* align structures from different toolchains (RDKit ↔ xTB ↔ ORCA ↔ OpenBabel)
* classify rotamers vs conformers with physical correctness


### Python CLI Usage
The iRSMD package comes with an CLI tool `irmsd`. This tool allows you to read multiple structures (e.g. from an extended `xyz`-format file), an perform operations on them. There are three subcommands that can be chosen: `prop`, `compare`, and `sort`/`prune`:
```
usage: irmsd [-h] [-v] {prop,compare,sort,prune} ...

CLI to read an arbitrary number of structures and run selected analysis commands on them.

positional arguments:
  {prop,compare,sort,prune}
                        Subcommand to run.
    prop                Compute structural properties (CN, rotational constants, canonical IDs).
    compare             Compare structures via iRMSD (default) or quaternion RMSD.
    sort (prune)        Sort, prune or cluster structures based on inter-structure measures. By
                        default, the more expensive iRMSD version is used. The use of the molecules'
                        energies is optional (--ethr) is optional but recommended. To fall back to
                        the quicker, but more empirical CREGEN workflow for ensemble sorting (using
                        energies, quaternion RMSDs and rotational constants), use --classic

options:
  -h, --help            show this help message and exit
  -v, --version         show program's version number and exit
```

The `prop` functionality exists as an utility function and can be used to determine some atomic properties for each structure provided:
```
usage: irmsd prop [-h] [--cn] [--rot] [--canonical] [--heavy] structures [structures ...]

positional arguments:
  structures   Paths to structure files (e.g. .xyz, .pdb, .cif).

options:
  -h, --help   show this help message and exit
  --cn         Calculate coordination numbers for each structure and print them as numpy arrays.
  --rot        Calculate the rotational constants.
  --canonical  Calculate the canonical identifiers.
  --heavy      When calculating canonical atom identifiers, consider only heavy atoms.
```

The `compare` subcommand performs a quaternion RMSD or an iRMSD comparison of the *first two* structures provided. It will return the alisgned structures.
```
usage: irmsd compare [-h] [--quaternion] [--inversion {on,off,auto}] [--heavy] [-o OUTPUT]
                     [--ref-idx REF_IDX] [--align-idx ALIGN_IDX]
                     structures [structures ...]

positional arguments:
  structures            Paths to structure files (e.g. .xyz, .pdb, .cif).

options:
  -h, --help            show this help message and exit
  --quaternion          Use the quaternion-based Cartesian RMSD instead of the invariant RMSD.
  --inversion {on,off,auto}
                        Control coordinate inversion in iRMSD runtypes: 'on', 'off', or 'auto'
                        (default: auto). Used only for iRMSD.
  --heavy               When comparing structures, consider only heavy atoms.
  -o OUTPUT, --output OUTPUT
                        Output file name (optional). If not provided, results are only printed.
  --ref-idx REF_IDX     Index of the reference structure in the provided structure list (default: 0,
                        i.e., the first structure).
  --align-idx ALIGN_IDX
                        Index of the structure to align to the reference structure (default: 1, i.e.,
                        the second structure). Used only for quaternion RMSD.
```

Finally, the `sort` (or `prune`) runtype performs ensemble pruning to remove redundant structures from a given structure list. It also splits the structure list into chemically distinct ensembles, should different molecules be included (currently decided via the sum formula).
```
usage: irmsd sort [-h] [--rthr RTHR] [--ethr [ETHR]] [--bthr BTHR] [--ewin EWIN]
                  [--inversion {on,off,auto}] [--align] [--classic] [--heavy] [--maxprint MAXPRINT]
                  [-o OUTPUT]
                  structures [structures ...]

positional arguments:
  structures            Paths to structure files (e.g. .xyz, .pdb, .cif).

options:
  -h, --help            show this help message and exit
  --rthr RTHR           Inter-structure RMSD threshold for sorting in Angström. Structures closer
                        than this threshold are treated as similar.
  --ethr [ETHR]         Inter-structure energy threshold in Hartree. If set, the default is 8.0e-5 Ha
                        (≈0.05 kcal/mol) or a user-specified value. Optional for iRMSD-based
                        runtypes.
  --bthr BTHR           Inter-structure rotational threshold used in the classical CREGEN sorting
                        procedure. The default is 0.01.
  --ewin EWIN           Energy window specification for CREGEN. Structures higher in energy than this
                        threshold (relative to the lowest energy structure in the ensemble) will be
                        removed. There is no default (all conformers are considered).
  --inversion {on,off,auto}
                        Control coordinate inversion when evaluating RMSDs during sorting: 'on',
                        'off', or 'auto' (default: auto). Only for iRMSD-based runtypes.
  --align               Just sort by energy and align.
  --classic, --cregen   Perform conformer classification with the CREGEN workflow based on a
                        comparison of quaternion RMSD, energy, interatomic distances, and rotational
                        constants. This routine is cheaper but more empirical than iRMSD-based
                        sorting. Does NOT restore mismatching atom order. Does not keep individual
                        rotamers.
  --maxprint MAXPRINT   Printout option; determine how many rows are printed for each sorted ensemble.
  -o OUTPUT, --output OUTPUT
                        Optional output file for sorted / clustered results.
```

For more information refer to the docs.

<br>

### Python Script Usage
A list of the provided functions and types can be found in the docs. A simple function could look like this:
```python
from irmsd import read_structures, sorter_irmsd_molecule

# read an xyz file with multiple frames and return a list of irmsd.Molecule objects
molecules = read_structures('/PATH/TO/YOUR/input.xyz')

# Analyze the structure list and assign each molecule to a group (simple list of integers)
# as well as the aligned molecules
groups,aligned_molecules = sorter_irmsd_molecule(molecules)
```

<br>

## Known Edge-Cases and Technical Limitations

* High-symmetry cases (e.g. C<sub>60</sub>, adamantane, etc.): Rotational axes are degenerate and an initial alignment is not possible this way.
* Interchage of atoms between molecules on different fragments in noncovalent complexes: The canonical assignment of two sub-graphs (two molecules that share no covalent connection) is currently independent. Hence, atoms (of the same element) may exchange across fragments, as for the (H<sub>2</sub>O)<sub>21</sub> or LJ<sub>75</sub> case. 
* Mismatches in canonical atom identifiers (comparing to chemical isomers that share a sum formula but have entirely different connectivity) are currently *not* caught and handled automatically.
* Some quaternion RMSDs may be slightly lower when comparing entirely *different* conformers (e.g. see example 4, figure b): This can occur due to the imperfect alignment+LSAP since rotational axes in different conformers varying orientations. Automatically falling back to the lower quaternion RMSD is an implementation TODO.              

<br>
<br>

## License
This project is licensed under the GNU Lesser General Public License (LGPL), version 3.0 or later.
You are free to use, modify, and redistribute the software under the terms of this license.
See the LICENSE file for the full text.

**Disclaimer:**  
This software is provided *“as is”*, without any warranty of any kind, express or implied,
including but not limited to warranties of merchantability, fitness for a particular purpose,
and noninfringement. In no event shall the authors or contributors be liable for any claim,
damages, or other liability arising from the use of this software.

© 2025 Philipp Pracht, Tobias Kaczun.<br> 
If you use this software in academic work, please acknowledge it and cite the [*associated publication*](https://doi.org/10.1021/acs.jcim.4c02143).

---

# OTHER TODOs
- [ ] docstrings and actual docs (GH pages?)
- [ ] conda-forge package
- [x] ci.yml/wheels.yml
- [x] codecov
- [ ] (implementation) Parallelization via OpenMP
- [x] (implementation) Optional pass of inter-conformer energy threshold (`ethr`)
- [ ] (implementation) Pre-alignment via quaternion RMSD of unique canonical atoms instead of aligning via rotational constants
- [x] (implementation) "classical" CREGEN pruner based on energy + rot.const. + quaternion RMSD
