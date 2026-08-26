# Provenance — published item to script to transcript

This is the file a referee should open first. It answers the only
question that matters to them: *where does this number come from, and how
do I reproduce it?*

Numbering follows the **submitted manuscript**. `MANIFEST.md` carries the
numbering of an earlier state of the text; it is kept as a provenance log,
not as a lookup table, and must not be followed by table number.

**Status: partly filled.** A row is filled only where the attribution has been
checked by hand against the transcript. The rest are deliberately left empty.
An empty cell is information; a wrong cell is not.

## Tables

| item | caption | script (`code/MEC_BLDC/`) | transcript (`outputs/`) | block used |
|---|---|---|---|---|
| Table 1 | Air-gap closure at fixed iron model. All network rows use the same lumped iron network; the last row is the finite-element reference | `RUN_CARTER_CMP.m` | `MEC_BLDC/A6_n2b_out.txt` | last complete |
| Table 2 | Positioning of the present formulation against the air-gap closures in use | — | — | narrative table, no chain |
| Table 3 | Data of the test machine, read directly from the finite-element project | `RUN_U2_MACHINE_TABLE.m` | `MEC_BLDC/U2_machine_table_out.txt` | last complete |
| Table 4 | Rank of the bore operator, Ms=360, PMSM data | | | |
| Table 5 | Hat basis: residual tail of the energy series beyond truncation N | | | |
| Table 6 | The increment diagnostic on a quantity known to carry the logarithmic tail. **Lower panel**, truncation unlocked from the tiling | `RUN_R4_UNLOCK.m`, then `RUN_G1_GAMMA.m` for the γ convention | `MEC_BLDC/R4_unlock_out.txt`, `MEC_BLDC/G1_gamma_out.txt` | last complete |
| Table 7 | Upper panel: the divergent bracket of (28), that is ln N + γ + ln\|2 sin(d/2)\|, along the two chains. Lower panel: truncation released from the tiling | upper: block X3; lower: `RUN_R4_UNLOCK.m` and `RUN_M4_HATBASIS.m` | `annulus_reference/X3_trace_table_out.txt`, `MEC_BLDC/R4_unlock_out.txt`, `MEC_BLDC/M4_hatbasis_out.txt` | last complete |
| Table 8 | Verification of Proposition 1 on a representative thin annulus | block C3 | `annulus_reference/C3_truekernel_out.txt` | last complete |
| Table 9 | Nonlinear solver, PMSM. Tightening the Newton tolerance | | | |
| Table 10 | Angular mesh convergence at one layer per shoe sub-region. The fundamental converges; the first slot sideband does not, and neither does the synchronous inductance | `RUN_A2_TABLE5.m`; L_a and M separately in `RUN_IND_MESH.m` | `MEC_BLDC/A2_table5_out.txt` 2nd block, `MEC_BLDC/C2_indmesh_out.txt` | 2nd complete |
| Table 11 | Radial refinement of the shoe, both solvers on identical meshes | `RUN_X1_TABLE5B_RECONCILE.m` | `MEC_BLDC/X1_table5b_reconcile_out.txt` | last complete |
| Table 12 | PMSM air-gap flux density at no load: complete spatial spectrum | `RUN_C4C_TABLE8.m` | `MEC_BLDC/C4c_table8_out.txt` | last complete |
| Table 13 | PMSM at no load: the mesh model at two shoe layerings, its lumped counterpart, and the reference | `RUN_A1_TABLE7.m` | `MEC_BLDC/A1_table7_out.txt` **last block** | last complete |
| Table 14 | Magnet eddy-current losses at no load. The two rows differ only in the source of the field | `RUN_A3_PMLOSS.m`; provenance of the published value in `RUN_D1_PMLOSS_CHAIN.m` | `MEC_BLDC/A3_pmloss_out.txt`, `MEC_BLDC/D1_pmloss_chain_out.txt` | last complete |
| Table 15 | Where the sweep deviation originates | | | |
| Table 16 | PMSM speed sweep at a 500 V dc bus, thirty points | `RUN_R7_SCORECARD.m` for the efficiency rows | `MEC_BLDC/R7_scorecard_out.txt`, points in `MEC_BLDC/G4_points_out.txt` | last complete |
| Table 17 | PMSM on-load transient, steady state | | | |
| Table 18 | PMSM model deviations placed against the resolution floor of the reference | | | |
| Table 19 | Harmonic comparison of the cogging torque | block C3 σ | `MEC_BLDC/C3_sigma_out.txt` | 1st complete |
| Table 20 | The twenty-nine quantities behind Fig. 10 | `RUN_R7_SCORECARD.m`; the two efficiency rows in points by `RUN_G4_POINTS.m` | `MEC_BLDC/R7_scorecard_out.txt`, `MEC_BLDC/G4_points_out.txt` | last complete |

## Figures

| item | caption | script (`code/MEC_BLDC/`) | transcript (`outputs/`) | block used |
|---|---|---|---|---|
| Figure 1 | Cross-sections of the two machines, to scale | | | |
| Figure 2 | Geometry of the two-layer annulus condensed by the operator | | | |
| Figure 3 | The four-branch cell | | | |
| Figure 4 | Radial layering of the stator mesh, developed along the bore | | | |
| Figure 5 | Solution chain | | | |
| Figure 6 | The two comparison axes, kept separate throughout | | | |
| Figure 7 | Air-gap flux density of the 15/14 PMSM at no load | `BLDC_MEC_COMPLET.m` | `MEC_BLDC/BLDC_MEC_COMPLET_out.txt` | last complete |
| Figure 8 | Internal magnetic circuit: magnetomotive force across each magnet and each gap | `BLDC_MEC_COMPLET.m` | `MEC_BLDC/BLDC_MEC_COMPLET_out.txt` | last complete |
| Figure 9 | Back-electromotive force and flux linkage at 1500 rpm | `BLDC_MEC_COMPLET.m` | `MEC_BLDC/BLDC_MEC_COMPLET_out.txt` | last complete |
| Figure 10 | Synthesis of the PMSM validation, twenty-nine quantities | `RUN_R7_SCORECARD.m` | `MEC_BLDC/R7_scorecard_out.txt` | last complete |
| Figure 11 | Saturation curves Lll(i) and ψab(i) from the map | | | |
| Figure 12 | On-load transient computed entirely by the MEC drive model | | | |
| Figure 13 | On-load transient, conversion and losses | | | |
| Figure 14 | Speed sweep at constant dc-bus voltage | | | |
| Figure 15 | Cogging torque of the PMSM | | | |

---

## How to fill one row

1. Open `MANIFEST.md` and find the quantity by **content**, not by number.
   Section 4 of that file carries the closure of 26 August and takes
   precedence over the earlier sections for the quantities it names.
2. Copy across the script and the output file.
3. Open the transcript and confirm the published value is there **at full
   precision**. The project rule is that the authoritative block is the
   **last complete** one, completeness meaning internal consistency
   checked, not merely position in the file.
4. If the value is not there, do not adjust the row. Report it. An empty
   cell is information; a wrong cell is not.

## Declared non-auditable

**One** published item cannot be re-formed from this archive:

- **the uncertainty band of the reference on field quantities** — it does
  not exist, and the manuscript states that it does not.

The 41.7-point bracket at n_sh = 1 **was** listed here and is no longer: the
full n_sh sweep is archived in `outputs/MEC_BLDC/X1_table5b_reconcile_out.txt`
and is read back at full precision by `RUN_G3_DUMP`. See `OPEN_POINTS.md`,
section *Closed*.

Both are described in `OPEN_POINTS.md`.
