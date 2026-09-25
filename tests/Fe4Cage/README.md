# Fe4Cage: the equal-mass axis fallback for mass-biased pre-alignment

## Test structure

`Fe4Cage.xyz`: a T-symmetric metal-organic cage, 316 atoms, 4 heavy Fe atoms among otherwise light organic ligands (C/H/N/O). 
Exactly the "fewer than 3 uniquely-ranked atoms" case where `min_rmsd` (`irmsd_module.f90`) falls back to principal-axis 
pre-alignment instead of a quaternion fit on unique atoms.

## Problem

The pre-alignment steps (`axis`/`axis_0`, `axis_module.f90`) weigh its inertia tensor by real isotopic mass. 
A small number of heavy atoms among many light ones can dominate that tensor and pull the resulting principal
axes away from the structure's actual ones.
This is a different failure mode from genuine axis degeneracy (e.g. C<sub>60</sub> and adamantane), 
which is still an open limitation where equial weighting does not help.

`pointgroup_isomers.json` holds 200 atom permutations, each an exact point-group operation of this structure (corresponding to a simple rotation).
Before the fix, `get_irmsd` scored only 170/200 of them correctly, with a worst case of ~2 Å, despite them being exact symmetry operations.

## Solution

Relabeling the 4 Fe atoms as Zn (also heavy) reproduces the failure at similar severity, 
while relabeling them as H (light) resolves the problem. This confirms that the failure stems from a mass mismatch.

`axis_0_equal_mass` / `axis_4_equal_mass` (`axis_module.f90`): 
identical to `axis_0`/`axis_4`, except every atom is weighted equally in the inertia tensor instead of by its real mass.
The only change is replacing `atmass = ams(at(i))` with `atmass = 1.0_wp`.

`min_rmsd` (`irmsd_module.f90`) now runs its existing pre-alignment + 32-way grid search (`min_rmsd_rotcheck_permute`) **twice**: 
once mass-weighted as before, and once with the equal-mass axes. 
It keeps whichever alignment reaches the lower LSAP cost. The grid search itself is untouched.

## Verification

- [`alignment_before_after.txt`](alignment_before_after.txt) -- the full,
  per-isomer before/after comparison for all 200 operations.
  All 30 that failed before are `C3` operations, each at exactly 2.007222 Å before and ~1.6×10⁻⁷ Å after the new fix.
- `test_fe4cage_alignment.py::test_all_pointgroup_isomers_score_near_zero`
  -- automated regression check: all 200/200 known-exact operations score
  < 1e-4 Å (worst case observed ~2×10⁻⁷ Å). Before the fix: 170/200, worst
  case ~2 Å.
- The full existing pytest suite (423 tests) passes unchanged.

