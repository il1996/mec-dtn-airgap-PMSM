# Reproducing the results

## 0. What you need

MATLAB R2024a, base product. No toolbox is required by the published chains.
No finite-element licence is required: the reference is shipped as exported
tables under `reference/ANSYS_750W/`.

## 1. Point the code at this checkout — once

The scripts are archived **unchanged**, byte for byte identical to the ones
that produced the published numbers. They therefore still carry the absolute
paths of the machine they ran on. Editing them before archiving would have
broken traceability, so they are adapted here instead, once, visibly and
reversibly.

```matlab
cd code
SET_REFERENCE_PATH('dryrun')   % shows what would change, writes nothing
SET_REFERENCE_PATH             % applies, and lists every file modified
```

A `.bak` file is written beside each modified file, so the operation can be
undone.

## 2. The two-minute test

```matlab
cd MEC_BLDC
RUN_R5_NSH
```

Expected: `B_g1 = 1.07787 T` at `n_sh = 1`, `1.07908 T` at `n_sh = 2`,
unknown counts 5400 and 6480, and the line `GARDE PASSEE` at the end.

Every script in this archive ends with a **guard** of this kind: a check that
would contradict the script's own result if the result were wrong. If a guard
does not print, do not trust the numbers above it.

## 3. Re-forming a published value

1. Find the item in [`PROVENANCE.md`](PROVENANCE.md).
2. Run the script named there.
3. Compare against the transcript named there, **not** against the paper: the
   transcript carries full precision, the paper carries one rounding.

Deviations in the paper are formed on full-precision values and rounded once,
at display. Re-forming a deviation from two rounded printed values will not in
general reproduce the printed deviation, and that is expected.

## 4. Reading a transcript

Each transcript opens with its configuration: machine, tiling `M_s`, radial
layers `n_sh`, surface basis, solver, fringing constant, and the reference path.

The authoritative block is the **last complete** block — complete meaning its
internal consistency has been checked, not merely that it is last. One archived
transcript contains four solved blocks of which the first three are broken by a
warm-start defect; rank alone is not a criterion.

## 4bis. Two files the scripts look for and do not need

**`kfringe_ident.mat`** is not in this archive, and its absence is not a
defect. Every script that looks for it falls back to `kfr = 0.325`, which is
the identified value used throughout the paper:

```matlab
if isfile('kfringe_ident.mat'), S0=load('kfringe_ident.mat'); kfr=S0.kbest;
else,                           kfr=0.325; end
```

The fringing constant is in any case **inert in the meshed chain**; it acts only
on the lumped one. Section 4.3 of the paper reports the sweep that establishes
this.

**Three `.mat` files are included** in `code/MEC_BLDC/` because the closure
blocks read them rather than re-solve: `A1_table7.mat` (the columns of
Table 13), `R7_scorecard.mat` (the twenty-nine compared quantities) and
`X1_table5b.mat` (the shoe-layering sweep). `RUN_G1_GAMMA` reads none of
them and recomputes its chain from scratch.

Two larger caches, `mec_map.mat` and `R6_satmap.mat`, are **not** included:
they are diagnostic inputs, none of the published values depends on them, and
they weigh 7.5 MB and 75 MB. The scripts that read them are shipped so that the
chain can be inspected, not so that it can be re-run without the cache.

## 5. If a number does not come out

Report it rather than adjusting it. The failure modes seen during assembly, in
order of frequency:

1. the wrong block of a transcript was read (see §4);
2. the deviation was re-formed on rounded values (see §3);
3. the configuration differs — check `M_s`, `n_sh`, `mu_r` and the rotor
   position count against the transcript header. Several tables in the paper
   are deliberately at different configurations, and each declares its own.
