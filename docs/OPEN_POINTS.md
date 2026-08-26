# Open points

What this archive does **not** settle. Stated here so that a reader does not
have to discover it by failing to reproduce something.

Two entries of the previous version of this file have since been closed and are
recorded at the end, under *Closed*, rather than deleted.

## 1. No uncertainty band on the reference for field quantities

The uncertainty budget bounds the finite-element reference on the power
quantities and on those alone. It produces **no** band on the fundamental
air-gap flux density, on the synchronous inductance, or on any other field
quantity, and none is asserted anywhere in the paper.

The consequence is stated rather than left to be inferred: a deviation of a few
tenths of a percent on a field quantity cannot be attributed to the model
rather than to the discretisation of the reference, because the size of the
second is unknown.

A second floor applies to the fundamental gap flux density and is now stated in
the paper: under the locked chain its tiling dispersion is **0.440 %** in the
piecewise-constant basis and 0.496 % in the hat basis, so a deviation below
those figures is smaller than the movement the tiling alone produces. The
0.28 % quoted in the abstract is reported as a bound and not as an agreement.

To close it: re-solve the reference project at two mesh densities and report
the movement of each field quantity between them.

## 2. Quantities reported as brackets, not as agreements

The first slot sideband and the cogging torque are governed by a re-entrant
tooth-tip corner singularity. The corner is sharp in the reference model too,
so the computed value depends on the cell size at the corner in the network and
in the finite-element mesh alike. These quantities are reported as brackets.
They are not validated, and the paper does not claim they are.

## 3. The synchronous inductance is not converged over the sweep examined

Over the angular refinement of `outputs/MEC_BLDC/A2_table5_out.txt` the
synchronous inductance deviates by +0.5 %, +0.2 %, +0.1 % and then **+2.5 %** at
the finest grid. `C2_indmesh_out.txt` shows the movement is carried by the phase
self-inductance itself, which rises from 50.379 to 51.549 mH between the last
two grids against a reference of 50.209 mH, and that it survives a second shoe
layer. It is therefore a property of the angular tiling and not of the
subtraction that forms L_d.

The deviation quoted in the abstract is the value in the configuration declared
for Table 13, not the limit of that sequence. The paper says so.

## 4. Two eddy-current formulations coexist for the iron loss

The published no-load iron loss, 21.3358 W, forms the eddy-current term on the
time derivative of the flux density. The same loss formed on the Bertotti
amplitude expression is 22.7409 W. The two differ by **6.59 %**
(`X2_ironloss_out.txt`). Both appear in the paper, at different places and for
different purposes, and each is now named where it is quoted — but a single
formulation would be better than two.

## 5. No experimental measurement

There is none, and none is claimed. The comparison throughout is between two
numerical models, which bounds their difference and never their common distance
from the physical machine.

## 6. The archive carries no DOI

No deposit has been made and **no identifier has been invented**. The manuscript
carries no unresolvable DOI: its data-availability statement says that the
persistent identifier is registered on acceptance and printed in the published
version.

To close it: mint the DOI, then write it, unchanged, in four places — the
manuscript, `.zenodo.json`, `README.md` and `CITATION.cff`.

---

## Closed

### The 41.7-point bracket at n_sh = 1 — **closed 26 August 2026**

The previous version of this file declared the 41.7-point bracket
non-auditable, on the ground that the Newton column of the n_sh = 1 execution
was not reproduced in the archive. **That was wrong, and the transcript is in
this archive.**

`outputs/MEC_BLDC/X1_table5b_reconcile_out.txt` (5 August 2026) carries the full
sweep, n_sh = 1 to 4 at M_s = 900 under both solvers. Read back at full
precision from `code/MEC_BLDC/X1_table5b.mat` by `RUN_G3_DUMP`:

| n_sh | ν = 8 linear | ν = 8 Newton | width |
|---|---|---|---|
| 1 | −9.717218 % | +31.982123 % | **41.699342 points** |
| 2 | −16.018999 % | +26.183413 % | **42.202412 points** |
| 3 | −17.987571 % | +24.526287 % | **42.513858 points** |
| 4 | −18.809286 % | +23.829427 % | **42.638713 points** |

All four published values are exact, and the argument that the bracket does not
close under radial refinement rests on the whole sweep rather than on a pair.

### The reference solver version — **closed 26 August 2026**

Read from the ANSYS project files and from the solution profiles, not inferred:

- **ANSYS Electronics Desktop 2023 R1, Maxwell 2-D, build 2023.1.0.**
- The no-load and armature-field studies are **magnetostatic**, with adaptive
  refinement at `PercentError = 1` %, `PercentRefinement = 30`,
  `MaximumPasses = 10`, `MinimumPasses = 2`; they converged in **2** and **3**
  passes respectively, both reporting *Adaptive Passes converged*.
- The back-electromotive-force, on-load and speed-sweep studies are
  **transient**. Maxwell 2-D does not adapt the mesh in transient, so the mesh
  they carry is the one written into the project. `NonlinearSolverResidual`
  is 1×10⁻⁴ throughout, on the unsmoothed B–H curve.
- Initial mesh TAU(2D), global surface-approximation setting 5, plus two
  length constraints: 10 mm over the stator and rotor bodies, and 1.75 mm
  over the magnets.

This is now recorded in the paper. **It strengthens the reservation of §7.2
rather than weakening it**: the transient reference was not mesh-converged by a
criterion, so what it can resolve is bounded by a mesh the author chose.
