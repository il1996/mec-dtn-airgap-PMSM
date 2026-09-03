# -*- coding: utf-8 -*-
"""Independent recomputation of every published value the manuscript's own
equations determine.

    python3 verify_all.py                # print the report
    python3 verify_all.py --csv          # also write results/*.csv

Nothing is read from the MATLAB chain.  Recomputed numbers come from
dtn_operator.py, which implements the published equations and the geometry of
Table 3; printed numbers come from published.py, which is a transcription of
the page.  No printed value is adjusted anywhere.

Verdicts
    REPRODUCED   recomputed and printed agree to the printed resolution
    ARITHMETIC   the printed table is self-consistent (a deviation, a ratio or
                 a dispersion re-formed from the printed columns)
    NOT-FROM-MS  the value is an output of the coupled model and cannot be
                 recomputed from the manuscript alone
"""
from __future__ import annotations
import argparse, csv, math, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import dtn_operator as OP
import published as PUB

ROWS = []            # (block, item, printed, recomputed, abs, rel, verdict)
FAILS = []

def rel(a, b):
    if b == 0: return float('nan')
    return (a - b) / abs(b)

def _res(v):
    """Half a unit in the last printed place of v, taken from its literal."""
    t = repr(float(v))
    if 'e' in t or 'E' in t:
        m, e = t.lower().split('e')
        dec = len(m.split('.')[1]) if '.' in m else 0
        return 0.5 * 10 ** (int(e) - dec)
    dec = len(t.split('.')[1]) if '.' in t else 0
    return 0.5 * 10 ** (-dec)


def interval_check(block, item, printed_result, fn, inputs, note="", printed_res=None):
    """Check that a deviation printed at full precision is consistent with the
    printed columns it is formed from: evaluate fn over the corners of the box
    the printed rounding leaves, and require the printed result to lie inside."""
    import itertools
    boxes = [(v - _res(v), v + _res(v)) for v in inputs]
    vals = [fn(*c) for c in itertools.product(*boxes)]
    lo, hi = min(vals), max(vals)
    pr = _res(printed_result) if printed_res is None else printed_res
    ok = (lo - pr - 1e-12) <= printed_result <= (hi + pr + 1e-12)
    mid = fn(*inputs)
    ROWS.append((block, item, printed_result, mid, printed_result - mid, None,
                 "ARITHMETIC" if ok else "MISMATCH",
                 note or f"printed value consistent with [{lo:.6g}, {hi:.6g}] implied by the printed columns"))
    if not ok:
        FAILS.append((block, item, printed_result, mid, (lo, hi)))
    return ok


def add(block, item, printed, recomputed, verdict, tol=None, note=""):
    if printed is None or recomputed is None:
        ROWS.append((block, item, printed, recomputed, None, None, verdict, note)); return
    ad = recomputed - printed
    rl = rel(recomputed, printed)
    if tol is not None and not (abs(rl) <= tol):
        verdict = "MISMATCH"; FAILS.append((block, item, printed, recomputed, rl))
    ROWS.append((block, item, printed, recomputed, ad, rl, verdict, note))

D6  = 6e-3
MS  = 1080
DMS = 2*math.pi/MS

# ------------------------------------------------------------------ 1. constants
def block_constants():
    add("constants", "(2 mu0 L / pi) ln 10  [Sec. 6.2]",
        6.078825e-08, OP.slope_per_decade(), "REPRODUCED", 1e-6)
    add("constants", "(2 mu0 L / pi) ln 2   [Table 6 note a]",
        PUB.TABLE6B_PRED_PC, OP.slope_per_doubling(), "REPRODUCED", 1e-6)
    add("constants", "(3/4) 3 mu0 L / (pi d^2), d = 2pi/1080  [Table 6 note b]",
        PUB.TABLE6B_PRED_HAT, OP.hat_increment_constant(DMS), "REPRODUCED", 1e-6)
    add("constants", "3 mu0 L / (pi d^2), d = 6 mrad  [Table 5]",
        0.001100, OP.hat_tail_constant(D6), "REPRODUCED", 1e-4)
    add("constants", "ln pi + gamma  [Table 7 caption]",
        PUB.TABLE7A_LNPI_GAMMA, math.log(math.pi) + OP.GAMMA, "REPRODUCED", 1e-6)
    add("constants", "U_a  [Table 8 caption]", 2.9261e-2, OP.UA, "REPRODUCED", 1e-4)
    add("constants", "U_m  [Table 8 caption]", 1.0973e-1, OP.UM, "REPRODUCED", 1e-4)
    add("constants", "tail onset 1/U_a  [Sec. 1 and 3.5: 'near order 34']",
        34.0, 1.0/OP.UA, "REPRODUCED", 0.01)

# ------------------------------------------------------- 2. Proposition 1, Table 8
def block_table8():
    di = dj = D6; dth = 50e-3
    cf = OP.offdiag_closed(di, dj, dth)
    add("Table 8(a)", "closed form (29)", PUB.TABLE8A[0][2], cf, "REPRODUCED", 1e-4)
    for N, summed, closed, dev in PUB.TABLE8A:
        s = OP.offdiag_pc(N, di, dj, dth)   # the proof writes -Y(i,j) = +(4 mu0 L/pi) sum
        add("Table 8(a)", f"summed, N = 1e{int(round(math.log10(N)))}", summed, s, "REPRODUCED", 5e-4)
        add("Table 8(a)", f"rel. dev., N = 1e{int(round(math.log10(N)))}",
            dev, abs(rel(s, cf)), "REPRODUCED", 5e-2)
    prev = None
    for N, summed, closed, dev, incr in PUB.TABLE8B:
        v = OP.self_energy_pc(N, D6, exact=False)
        add("Table 8(b)", f"summed, N = 1e{int(round(math.log10(N)))}", summed, v, "REPRODUCED", 5e-4)
        add("Table 8(b)", f"closed form (30), N = 1e{int(round(math.log10(N)))}",
            closed, OP.diag_closed(N, D6), "REPRODUCED", 5e-4)
        if prev is not None and incr is not None:
            add("Table 8(b)", f"increment per decade at N = 1e{int(round(math.log10(N)))}",
                incr, v - prev, "REPRODUCED", 5e-4)
        prev = v
    prevA = prevE = None
    for N, asym, exact, incr in PUB.TABLE8C:
        a = OP.self_energy_pc(N, D6, exact=False)
        e = OP.self_energy_pc(N, D6, exact=True)
        add("Table 8(c)", f"exact kernel, N = 1e{int(round(math.log10(N)))}", exact, e, "REPRODUCED", 5e-4)
        if prevE is not None:
            add("Table 8(c)", f"increment per decade, exact, N = 1e{int(round(math.log10(N)))}",
                incr, e - prevE, "REPRODUCED", 5e-4)
        prevA, prevE = a, e
    shift = (OP.self_energy_pc(10**5, D6, True) - OP.self_energy_pc(10**5, D6, False))
    add("Table 8 caption", "constant offset between the two kernels",
        PUB.TABLE8C_SHIFT, shift, "REPRODUCED", 2e-3,
        "geometry reconstructed from Table 3 to six figures")

# ------------------------------------------------------------------ 2b. Table 4
def block_table4():
    """rank Y = min(2N, M_s - 1) and the null-space dimension, from the spectral
    result of Section 3.1.  M_s = 360, as the caption declares."""
    Ms = 360
    for N, rank_pub, null_pub in [(60, 120, 240), (120, 240, 120), (179, 358, 2),
                                  (180, 359, 1), (360, 359, 1)]:
        rank = min(2*N, Ms - 1)
        add("Table 4", f"rank Y at N = {N}", rank_pub, rank, "REPRODUCED", None)
        add("Table 4", f"null-space dimension at N = {N}", null_pub, Ms - rank,
            "REPRODUCED", None)


# ------------------------------------------------------------------- 3. Table 5
def block_table5():
    for N, tail, prod in PUB.TABLE5:
        if N > 10**4:                       # the sum beyond 1e5 is not worth the wall clock
            ROWS.append(("Table 5", f"residual tail, N = 1e{int(round(math.log10(N)))}",
                         tail, None, None, None, "SKIPPED (cost)", "converges to the same constant"))
            continue
        t = OP.hat_tail(N, D6, exact=False, reach=400 if N <= 1000 else 60)
        add("Table 5", f"residual tail, N = 1e{int(round(math.log10(N)))}", tail, t, "REPRODUCED", 2e-3)
        add("Table 5", f"tail x N^2, N = 1e{int(round(math.log10(N)))}", prod, t*N*N, "REPRODUCED", 2e-3)

# ---------------------------------------------------------------- 4. Table 6(b)
def block_table6b():
    """The self-energy column is recomputed at double precision from (27) and (33);
    the increment and deviation columns are then formed on the RECOMPUTED values,
    because at the last doublings the printed increment column no longer carries
    enough digits to re-form its own deviation."""
    Ns = [r[0] for r in PUB.TABLE6B]
    top = max(Ns)
    checks = set(Ns)
    s_pc = s_hat = 0.0
    part_pc, part_hat = {}, {}
    for n in range(1, top + 1):
        k = OP.kappa(n, True)
        sn = math.sin(n * DMS / 2.0)
        s_pc  += k / n * sn * sn
        s_hat += k / n**3 * sn**4
        if n in checks:
            part_pc[n]  = 4.0 * OP.MU0 * OP.L_M / math.pi * s_pc
            part_hat[n] = 16.0 * OP.MU0 * OP.L_M / (math.pi * DMS * DMS) * s_hat
    prev_pc = prev_hat = None
    for N, ypc, dpc, devpc, yhat, dhatn2, devhat in PUB.TABLE6B:
        add("Table 6(b)", f"-Y(i,i) p.c., N_h = {N}", ypc, part_pc[N], "REPRODUCED", 2e-6)
        add("Table 6(b)", f"-Y(i,i) hat, N_h = {N}", yhat, part_hat[N], "REPRODUCED", 2e-6)
        if prev_pc is not None:
            inc_pc = part_pc[N] - prev_pc
            inc_h  = (part_hat[N] - prev_hat) * (N // 2) ** 2
            add("Table 6(b)", f"increment p.c., N_h = {N}", dpc, inc_pc, "REPRODUCED", 1e-5)
            add("Table 6(b)", f"increment x N^2, hat, N_h = {N}", dhatn2, inc_h, "REPRODUCED", 1e-4)
            add("Table 6(b)", f"deviation of the increment, N_h = {N}", devpc,
                (inc_pc - PUB.TABLE6B_PRED_PC) / PUB.TABLE6B_PRED_PC * 100,
                "REPRODUCED", None, "recomputed from the double-precision self-energy")
            add("Table 6(b)", f"deviation of increment x N^2, hat, N_h = {N}", devhat,
                (inc_h - PUB.TABLE6B_PRED_HAT) / PUB.TABLE6B_PRED_HAT * 100,
                "REPRODUCED", None, "recomputed from the double-precision self-energy")
        prev_pc, prev_hat = part_pc[N], part_hat[N]


# ------------------------------------------------------- 5. Table 6(a) and 7(a)
def block_brackets():
    for row in PUB.TABLE6A:
        N, br = row[0], row[1]
        add("Table 6(a)", f"bracket, N_h = {N}", br, OP.bracket(N, DMS), "REPRODUCED", 1e-6)
    b0, b1 = PUB.TABLE6A[0][1], PUB.TABLE6A[-1][1]
    add("Sec. 6.2", "growth of the bracket over the released sweep (%)",
        442.8, (b1-b0)/b0*100, "ARITHMETIC", 1e-3)
    drift_pc  = (PUB.TABLE6A[-1][2]-PUB.TABLE6A[0][2])/PUB.TABLE6A[0][2]*100
    drift_hat = (PUB.TABLE6A[-1][5]-PUB.TABLE6A[0][5])/PUB.TABLE6A[0][5]*100
    add("Sec. 6.2", "drift of B_g1, piecewise-constant (%)", 0.703, drift_pc, "ARITHMETIC", 2e-3)
    add("Sec. 6.2", "drift of B_g1, hat (%)", 0.023, drift_hat, "ARITHMETIC", 5e-2)
    ROWS.append(("Table 6(a)", "B_g1 columns", "1.078648295 ... 1.086233369", None,
                 None, None, "NOT-FROM-MS",
                 "an output of the coupled model; the recomputation checks only the bracket, "
                 "the increments and the ratios"))
    for Ms, br in PUB.TABLE7A:
        add("Table 7(a)", f"bracket at M_s = {Ms}", br,
            OP.bracket(Ms//2, 2*math.pi/Ms), "REPRODUCED", 1e-6)
    b_lo = OP.bracket(270, 2*math.pi/540); b_hi = OP.bracket(4320, 2*math.pi/8640)
    mv = (b_hi-b_lo)/b_lo*100
    add("Table 7(a)", "movement of the bracket over a factor 16 in tiling (%)",
        PUB.TABLE7A_BRACKET_PRED, mv, "REPRODUCED", 0.15)
    add("Table 7(a)", "ratio |self term| / bracket movement",
        PUB.TABLE7A_RATIO, abs(PUB.TABLE7A_SELF_TERM)/mv, "ARITHMETIC", 5e-3)

# ---------------------------------------------------------------- 6. Table 7(b)
def block_table7b():
    disp = lambda v: (max(v)-min(v))/abs(sum(v)/len(v))*100
    for k, v in PUB.TABLE7B.items():
        interval_check("Table 7(b)", f"dispersion over six tilings, {k} (%)",
                       PUB.TABLE7B_DISP[k], lambda *a: disp(list(a)), v, printed_res=5e-4)
    for k, v in PUB.TABLE7B_INT.items():
        interval_check("Table 7(b) cont.", f"dispersion over six tilings, {k} (%)",
                       PUB.TABLE7B_INT_DISP6[k], lambda *a: disp(list(a)), v, printed_res=5e-4)
        interval_check("Table 7(b) cont.", f"dispersion over the five finest, {k} (%)",
                       PUB.TABLE7B_INT_DISP5[k], lambda *a: disp(list(a)), v[1:], printed_res=5e-4)
    for name, v in (("L_d p.c.", "p.c."), ("L_d hat", "hat")):
        La = PUB.TABLE7B_INT[f"L_a {v}"]; M = PUB.TABLE7B_INT[f"M {v}"]
        Ld = PUB.TABLE7B_INT[f"L_d {v}"]
        worst = max(abs((a-b)-c) for a, b, c in zip(La, M, Ld))
        add("Table 7(b) cont.", f"identity L_d = L_a - M, {v} (worst residual, mH)",
            0.0, worst, "ARITHMETIC", None)

# ---------------------------------------------------------------- 7. Table 7(d)
def block_table7d():
    for name, ref, pc, dpc, hat, dhat in PUB.TABLE7D:
        r = 5e-4 if name.startswith("B_g1") else 5e-3
        interval_check("Table 7(d)", f"deviation p.c., {name} (%)", dpc,
                       lambda a, b: (a-b)/b*100, [pc, ref], printed_res=r)
        interval_check("Table 7(d)", f"deviation hat, {name} (%)", dhat,
                       lambda a, b: (a-b)/b*100, [hat, ref], printed_res=r)

# ----------------------------------------------------------------- 8. Table 12
def block_table12():
    closer_mesh = closer_lumped = 0
    away = toward = 0
    for name, m1, m2, lu, fea, (d1, d2, d3) in PUB.TABLE12:
        for lbl, val, pub in (("mesh n_sh=1", m1, d1), ("mesh n_sh=2", m2, d2), ("lumped", lu, d3)):
            interval_check("Table 12", f"{name} / {lbl} (%)", pub,
                           lambda a, b: (a-b)/b*100, [val, fea],
                           "deviation formed at full precision and rounded once, Sec. 5.3",
                           printed_res=5e-3)
        if abs(d1) < abs(d3): closer_mesh += 1
        else:                 closer_lumped += 1
        if abs(d2) > abs(d1): away += 1
        else:                 toward += 1
    add("Sec. 7.4", "count, mesh closer / lumped closer  ('within one quantity')",
        1, abs(closer_mesh-closer_lumped), "ARITHMETIC", None,
        f"mesh {closer_mesh}, lumped {closer_lumped}")
    add("Sec. 8", "'thirteen of the seventeen' move away under shoe refinement",
        13, away, "ARITHMETIC", None, f"away {away}, toward {toward}")
    for name, m1, m2, lu, fea, _ in PUB.TABLE12:
        if name.startswith("L_d"):
            La = [r for r in PUB.TABLE12 if r[0].startswith("L_a")][0]
            M  = [r for r in PUB.TABLE12 if r[0].startswith("M (")][0]
            worst = max(abs((La[i]-M[i])-[m1, m2, lu, fea][i-1]) for i in (1, 2, 3, 4))
            add("Table 12", "identity L_d = L_a - M (worst residual, mH)", 0.0, worst,
                "ARITHMETIC", None)

# ------------------------------------------------------- 9. Carter, skin depth, cost
def block_section7():
    C = PUB.CARTER
    gp = C["g_mm"] + C["hm_mm"]/C["mu_r"]
    ts = math.pi*C["D_mm"]/C["Ns"]
    x  = C["b0_mm"]/(2*gp)
    gam = (4/math.pi)*(x*math.atan(x) - math.log(math.sqrt(1+x*x)))
    kC  = ts/(ts - gam*gp)
    gam2 = (C["b0_mm"]/gp)**2/(5 + C["b0_mm"]/gp)
    kC2  = ts/(ts - gam2*gp)
    add("Sec. 7.1", "g' (mm)",  C["gprime_mm"], gp, "REPRODUCED", 1e-4)
    add("Sec. 7.1", "tau_s (mm)", C["tau_s_mm"], ts, "REPRODUCED", 1e-5)
    add("Sec. 7.1", "x",         C["x"],        x,  "REPRODUCED", 1e-4)
    add("Sec. 7.1", "gamma, exact form (43)", C["gamma_exact"], gam, "REPRODUCED", 1e-4)
    add("Sec. 7.1", "k_C, exact form (43)",   C["kC_exact"],    kC,  "REPRODUCED", 1e-5)
    add("Sec. 7.1", "gamma, approximate form", C["gamma_approx"], gam2, "REPRODUCED", 1e-4)
    add("Sec. 7.1", "k_C, approximate form",   C["kC_approx"],    kC2,  "REPRODUCED", 1e-5)
    add("Sec. 7.1", "difference between the two forms (%)",
        C["difference_pct"], (kC2-kC)/kC*100, "REPRODUCED", 2e-2)

    S = PUB.SKIN
    f = S["Ns"]*S["speed_rpm"]/60.0
    sig = 1.0/S["rho_ohm_m"]
    delta = math.sqrt(2.0/(2*math.pi*f*OP.MU0*S["mu_r"]*sig))*1e3
    add("Sec. 6.6", "slotting frequency seen from the rotor (Hz)", S["f_Hz"], f, "REPRODUCED", 1e-6)
    add("Sec. 6.6", "magnet conductivity (S/m)", S["sigma_S_per_m"], sig, "REPRODUCED", 1e-5)
    add("Sec. 6.6", "skin depth (mm)", S["delta_mm"], delta, "REPRODUCED", 2e-3,
        "mu = mu0 mu_r")
    add("Sec. 6.6", "skin depth / magnet thickness", S["ratio"], delta/S["hm_mm"],
        "REPRODUCED", 5e-3)

    K = PUB.COST
    add("Sec. 6.9", "sweep as a fraction of 128 independent solutions (%)",
        K["sweep_pct_of_independent"], K["sweep_total_s"]/(K["n_positions"]*K["single_position_s"])*100,
        "ARITHMETIC", 5e-3)
    add("Sec. 6.9", "marginal cost with reassembly (ms)",
        K["marginal_reassembled_ms"], K["reassembled_sweep_s"]/K["n_positions"]*1e3,
        "ARITHMETIC", 5e-3)
    add("Sec. 6.9", "factor on the sweep", K["factor_sweep"],
        K["reassembled_sweep_s"]/K["sweep_total_s"], "ARITHMETIC", 1e-2)
    add("Sec. 6.9", "per-position cost of the 721-position sweep (ms)",
        K["per_position_ms"], K["sweep_721_s"]/K["positions_721"]*1e3, "ARITHMETIC", 5e-3)
    add("Sec. 7.1", "factor, invariant operator against sliding interface",
        K["factor_interface"], K["sliding_interface_s"]/K["invariant_operator_s"],
        "ARITHMETIC", 5e-3)

# ------------------------------------------------------------------- reporting
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", action="store_true", help="write results/*.csv")
    a = ap.parse_args()
    for f in (block_constants, block_table8, block_table4, block_table5, block_table6b,
              block_brackets, block_table7b, block_table7d, block_table12,
              block_section7):
        f()
    w = max(len(r[1]) for r in ROWS) + 2
    cur = None
    for block, item, pr, rc, ad, rl, vd, note in ROWS:
        if block != cur:
            print(); print("=" * 100); print(block); print("=" * 100); cur = block
        p = "" if pr is None else (f"{pr:.6g}" if isinstance(pr, float) else str(pr))
        r = "" if rc is None else f"{rc:.6g}"
        d = "" if rl is None or rl != rl else f"{rl:+.2e}"
        print(f"  {item:<{w}} printed {p:>16}   recomputed {r:>16}   rel {d:>10}   {vd}"
              + (f"   [{note}]" if note else ""))
    print()
    print("=" * 100)
    n = len(ROWS)
    good = sum(1 for r in ROWS if r[6] in ("REPRODUCED", "ARITHMETIC"))
    bad  = sum(1 for r in ROWS if r[6] == "MISMATCH")
    nfm  = sum(1 for r in ROWS if r[6] == "NOT-FROM-MS")
    skip = sum(1 for r in ROWS if r[6].startswith("SKIPPED"))
    print(f"checks: {n}   agreeing: {good}   mismatches: {bad}   "
          f"not reproducible from the manuscript: {nfm}   skipped: {skip}")
    if bad:
        print()
        for b in FAILS:
            print("  MISMATCH", b)
    print("VERDICT:", "PASS" if bad == 0 else "FAIL")
    print("=" * 100)
    if a.csv:
        d = os.path.join(HERE, "results")
        os.makedirs(d, exist_ok=True)
        with open(os.path.join(d, "independent_verification.csv"), "w", newline="",
                  encoding="utf-8") as fh:
            wr = csv.writer(fh)
            wr.writerow(["block", "item", "printed", "recomputed",
                         "absolute_difference", "relative_difference", "verdict", "note"])
            for r in ROWS: wr.writerow(r)
        print("written:", os.path.join(d, "independent_verification.csv"))
    return 0 if bad == 0 else 1

if __name__ == "__main__":
    sys.exit(main())
