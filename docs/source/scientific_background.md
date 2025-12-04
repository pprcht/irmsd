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

$$
\mathrm{iRMSD}(\mathbf{X}, \mathbf{Y}) = \min_{\mathbf{P},\mathbf{U}} \sqrt{ \frac{1}{N} \sum_{i=1}^{N} \left\lVert \mathbf{X}_i - (\mathbf{P}\mathbf{U}\mathbf{Y})_i \right\rVert^2 }
$$

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

```{image} assets/images/fluoxetine1.jpg
:width: 50%
:align: center
```

   ```
   40

C       -0.0198       0.2158       0.5308
F       -5.0246      -0.2464       0.5593
H       -2.0207       0.0761       3.2888
H        2.0221       5.0303       0.4892
C       -0.4032       0.2003       1.8748
H        1.4542       5.1890      -1.1937
C        2.1773       2.6686      -0.7626
F       -4.6775       0.8294       2.4172
F       -4.3667      -1.3218       2.3289
N        3.0555       3.8290      -0.9233
C        2.3371       5.0467      -0.5601
H        3.7051       1.1947      -0.3434
O        1.3305       0.2076       0.3364
C       -2.3576       0.1386      -0.0950
C        3.3367      -3.7475      -2.0944
C       -2.7426       0.0481       1.2484
C       -4.1855      -0.1659       1.6292
C        2.4651      -1.1750      -1.3319
H        1.7906       2.6512       0.2660
C        2.8963       1.3473      -1.0679
H        1.3203       2.7752      -1.4405
C        1.8806       0.1908      -0.9905
H       -0.7865       0.2662      -1.5101
H        3.3527       1.3925      -2.0651
H        1.1282       0.3962      -1.7541
H        2.8529      -0.6554      -3.3999
H        2.9547      -4.3108      -0.0530
H        3.8626       3.7250      -0.3037
H        0.3581       0.2294       2.6522
C       -1.7511       0.1165       2.2345
H        2.9886       5.9146      -0.7002
H        2.1980      -2.0723       0.6278
H        3.6193      -2.9009      -4.0530
H       -3.1084       0.1174      -0.8845
C        2.5162      -2.2260      -0.4015
C       -1.0075       0.2216      -0.4508
H        3.6502      -4.7430      -2.3936
C        2.9433      -3.5025      -0.7811
C        2.8862      -1.4373      -2.6475
C        3.3158      -2.7117      -3.0256
   ```

  </td><td>

```{image} assets/images/fluoxetine2.jpg
:width: 50%
:align: center
```

   ```
   40

C       -0.2470      -1.0143      -0.4213
H        0.0828       4.2564      -4.4517
C       -2.1193       2.6178       0.0144
C        1.2089      -1.8204       1.8160
C        0.0619      -0.0804       0.5643
C        0.9632      -2.7436       0.7925
H       -5.9431       1.8117       2.4085
H       -3.2933       0.8099       0.1366
H        2.4555       4.5828      -3.8102
H       -0.0328      -3.0384      -1.1072
H       -3.9856       3.0777       1.8789
H       -2.4640       1.2332       1.6473
C       -0.4255       3.2037      -2.6383
H        1.7581       2.5331      -0.1112
C        1.3640       2.9121      -1.0521
C        1.5502      -4.1306       0.8575
C       -5.0208       1.3357       2.0622
H       -4.4695       0.9920       2.9447
F        1.1253      -4.8340       1.9430
C        2.2338       3.5914      -1.9111
H       -1.7744       3.4140       0.6847
C        0.4427       3.8850      -3.4947
C        0.7599      -0.5019       1.6995
F        1.2525      -4.8990      -0.2270
H       -1.3360       1.0931      -1.2358
H        1.7763      -2.1151       2.6974
H        0.9785       0.2106       2.4926
H       -5.3029       0.4654       1.4587
H       -0.8129      -0.7679      -1.3112
H       -2.6932       3.0927      -0.7916
C       -0.9181       1.8428      -0.5615
F        2.9112      -4.1100       0.9336
H       -1.4519       3.0470      -2.9556
H        3.2752       3.7325      -1.6299
C        1.7755       4.0721      -3.1352
C        0.0230       2.6924      -1.4078
N       -4.2396       2.2887       1.2797
O       -0.2303       1.2522       0.5528
C       -3.0176       1.6462       0.7923
C        0.1996      -2.3349      -0.3081
   ```       

  </td>
  </tr>
</table>

<br>
<br>

##  Examples  📊

### Example 1 - Symmetric Rotamers (Pentane C<sub>5</sub>H<sub>12</sub>)

```{image} assets/images/example1.jpg
:width: 100%
:align: center
```
Most pentane conformers have *four* rotamers that differ only by permutation or inversion of terminal methyl groups.
As shown in Figure 3 of the paper (page 3) :

* Conventional RMSD between these rotamers ranges from 1.2-1.8 Å
* iRMSD correctly identifies all as identical (≈ 0.0 Å) by accounting for permutations and false enantiomers
* Property calculations (IR spectra, NMR shifts) also confirm they are physically equivalent

This simple example illustrates a key problem with classical (quaternion) RMSD-based conformer comparison and the necessity of addressing *both* the alignment *and* permutation problems for chemical workflows.

<br>

### Example 2 - Validation on Randomized Atom Order Structures

```{image} assets/images/example2.jpg
:width: 100%
:align: center
```

A robust permutation-handling alignment algorithm must correctly classify structures as identical even if:

* randomly rotated
* atom order is randomly permuted
* (optional) mirrored ("*false enantiomers*")

Figures 7a–d (pages 6–7) show that iRMSD successfully returns ~0 Å for **every pair** in 100 randomized input coordinates of: pentane (12 atoms), TPPO (46 atoms), taxol (113 atoms), BPTI (892 atoms).


<br>

### Example 3 - Noncovalent Clusters (LJ<sub>75</sub> and (H<sub>2</sub>O)<sub>21</sub>)


```{image} assets/images/example3.jpg
:width: 100%
:align: center
```

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

```{image} assets/images/example4.jpg
:width: 100%
:align: center
```

iRMSD excels in distinguishing on a single threshold parameter (`RTHR`):

* rotamers (same conformer → iRMSD ≈ 0 Å)
* different conformers (larger iRMSD values)

This is crucial for automated CRE pruning and is an extension to conventional (quaternion) RMSD pruning, e.g. as used in CREST.<br>
The default `RTHR` threshold in `irmsd` to distinguish to structures as conformers is **0.125 Å**, which was adapted from CREST's CREGEN procedure. Additional thresholds, e.g. for the inter-conformer energy difference (`ETHR`) or rotational constants (`BTHR`) are *not* required, but can be used to achieve more efficient pre-sorting.


