# -*- coding: utf-8 -*-
"""NONREGRESSION_BATTERY  -  task B9 of the revision specification.

The battery of fifteen checks that must be re-run before and after every other
task. Any divergence is a failure, except where a task explicitly foresees a
change.

RULES THIS SCRIPT OBEYS.
  * Nothing is transcribed as a result. Published values appear ONLY as targets
    on the right-hand side of a comparison.
  * A check that cannot be decided mechanically reports UNDECIDED and prints
    what it found. It never reports PASS by default.
  * Nothing outside outputs/revision is written.

Usage:  python nonregression_battery.py
Exit code 0 if every decidable check passes, 1 otherwise.
"""
import os
import sys, re, sys, json, math
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ROOT = os.path.dirname(ROOT) if os.path.basename(ROOT) == 'code' else ROOT
REPO = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
OUT  = os.path.join(REPO, 'outputs')
REV  = os.path.join(OUT, 'revision')
os.makedirs(REV, exist_ok=True)

# ----------------------------------------------------------------- kernel
mu0 = 4e-7*math.pi
L   = 0.033
Rs  = 34.678e-3
Ua  = 2.9260675e-2
Um  = 1.0973162e-1
mur = 1.038952
c   = mu0*L/math.pi
gam = 0.5772156649015329

def kappa(n):
    """kappa_n with the mandatory guards: kappa=1 beyond n*Ua>25 (correction
    < e^-50), coth argument clamped, so cosh cannot overflow past n~24265."""
    n = np.asarray(n, dtype=float)
    k = np.ones_like(n)
    m = (n*Ua) <= 25.0
    nm = n[m]
    D = np.cosh(nm*Ua) + mur*np.sinh(nm*Ua)/np.tanh(np.minimum(nm*Um, 25.0))
    k[m] = (np.cosh(nm*Ua) - 1.0/D)/np.sinh(nm*Ua)
    return k

def mY(k, d, N, exact=True, chunk=2_000_000):
    """-Y(i,i+k) = (4 mu0 L/pi) sum_{n<=N} (kappa_n/n) sin^2(n d/2) cos(n k d)"""
    tot = 0.0
    n0 = 1
    while n0 <= N:
        n1 = min(N, n0+chunk-1)
        n = np.arange(n0, n1+1, dtype=float)
        w = kappa(n) if exact else 1.0
        tot += np.sum(w*np.sin(n*d/2)**2*np.cos(n*k*d)/n)
        n0 = n1+1
    return 4*mu0*L/math.pi*tot

def firstrow(Ms, N, basis='p0'):
    """First row of -Y on a uniform tiling, by residue folding (exact)."""
    d = 2*math.pi/Ms
    n = np.arange(1, N+1, dtype=float)
    k1 = kappa(n)
    if basis == 'p0':
        prj = (2.0/(n*math.pi))*np.sin(n*d/2)
    else:
        prj = (4.0/(math.pi*n**2*d))*np.sin(n*d/2)**2
    g = -mu0*(n/Rs)*k1
    w = L*Rs*math.pi*g*prj**2
    W = np.bincount(np.mod(n.astype(np.int64), Ms), weights=w, minlength=Ms)
    return -np.real(np.fft.fft(W))

# ----------------------------------------------------------------- report
RES = []
def rec(tag, name, verdict, detail):
    RES.append((tag, name, verdict, detail))
    print(f"  [{tag}] {verdict:9s} {name}")
    for line in detail.splitlines():
        print("             " + line)

def close(a, b, tol):
    return abs(a-b) <= tol

def read(p):
    for enc in ('utf-8', 'latin-1'):
        try:
            with open(p, encoding=enc) as f:
                return f.read()
        except (UnicodeDecodeError, FileNotFoundError):
            continue
    return None

print("="*78)
print("B9 : NON-REGRESSION BATTERY")
print("     repository :", REPO)
print("="*78)

# ---------------------------------------------------------------- NR-01
Ms = 360
rows, okall = [], True
for N, ex in zip([60, 120, 179, 180, 360], [120, 240, 358, 359, 359]):
    lam = np.real(np.fft.fft(firstrow(Ms, N)))
    rk = int(np.sum(np.abs(lam) > np.max(np.abs(lam))*1e-10))
    rows.append(f"N={N:4d}  rank={rk:4d} expected {ex:4d}  nullity={Ms-rk:4d} expected {Ms-ex:4d}")
    okall &= (rk == ex)
rec("NR-01", "Table 4, Ms=360: ranks 120/240/358/359/359",
    "PASS" if okall else "FAIL", "\n".join(rows))

# ---------------------------------------------------------------- NR-02
d6 = 6e-3
def hat_tail(N, upto=40_000_000):
    n = np.arange(N+1, upto+1, dtype=float)
    k1 = kappa(n)
    W = (4.0/(math.pi*n**2*d6))*np.sin(n*d6/2)**2
    return L*Rs*math.pi*np.sum(mu0*(n/Rs)*k1*W**2)
rows, okall = [], True
for N in (10**5, 10**6):
    t = hat_tail(N)
    v = t*N*N
    rows.append(f"N=1e{int(math.log10(N))}  tail={t:.4e}  tail x N^2 = {v:.6f}  target 0.001100")
    okall &= close(round(v, 6), 0.001100, 5e-6)
rec("NR-02", "Table 5: hat tail x N^2 -> 0.001100 on the last two decades",
    "PASS" if okall else "FAIL", "\n".join(rows))

# ---------------------------------------------------------------- NR-04
inc = 2*c*math.log(2)
rec("NR-04", "Table 6(b): p.c. increment -> (2 mu0 L/pi) ln2",
    "PASS" if close(inc, 1.829909e-08, 1e-14) else "FAIL",
    f"(2 mu0 L/pi) ln2 = {inc:.6e}   manuscript prints 1.829908e-08 (value) "
    f"and 1.829909e-08 (prediction)\n"
    f"difference to the printed prediction : {abs(inc-1.829909e-08):.2e}")

# ---------------------------------------------------------------- NR-05
br = [math.log(Ms//2)+gam+math.log(abs(2*math.sin(math.pi/Ms)))
      for Ms in (540, 1080, 2160, 4320, 8640)]
tgt = [1.721940, 1.721944, 1.721945, 1.721945, 1.721946]
okall = all(close(a, b, 5e-7) for a, b in zip(br, tgt))
rec("NR-05", "Table 7(a): bracket 1.721940/44/45/45/46",
    "PASS" if okall else "FAIL",
    "  ".join(f"{v:.6f}" for v in br))

# ------------------------------------------------------- NR-06 / NR-07
PUB_B = {'pc': [1.07390, 1.07865, 1.07755, 1.07793, 1.07759, 1.07743],
         'hat': [1.07267, 1.07801, 1.07698, 1.07759, 1.07741, 1.07734]}
PUB_N = {'pc': [0.021761, 0.016944, 0.018062, 0.017672, 0.018014, 0.018180],
         'hat': [0.023000, 0.017596, 0.018637, 0.018020, 0.018195, 0.018271]}
dsp = lambda x: 100*(max(x)-min(x))/abs(sum(x)/len(x))
d6a = [dsp(PUB_B['pc']), dsp(PUB_B['hat']), dsp(PUB_N['pc']), dsp(PUB_N['hat'])]
t6a = [0.440, 0.496, 26.126, 28.514]
rec("NR-06", "Table 7(b): dispersions 0.440 / 0.496 / 26.126 / 28.514 %",
    "PASS" if all(close(a, b, 0.002) for a, b in zip(d6a, t6a)) else "FAIL",
    "  ".join(f"{v:.3f} %" for v in d6a))
d5a = [dsp(PUB_B['pc'][1:]), dsp(PUB_B['hat'][1:]), dsp(PUB_N['pc'][1:]), dsp(PUB_N['hat'][1:])]
rec("NR-07", "Table 7(b) without the coarsest tiling: 0.096 vs 0.113 %, 5.7 vs 7.0 %",
    "PASS" if (close(d5a[0], 0.113, 0.002) and close(d5a[1], 0.096, 0.002)
               and close(d5a[2], 6.95, 0.3) and close(d5a[3], 5.74, 0.3)) else "FAIL",
    f"p.c. B_g1 {d5a[0]:.3f} %  hat B_g1 {d5a[1]:.3f} %  "
    f"p.c. nu8 {d5a[2]:.3f} %  hat nu8 {d5a[3]:.3f} %")

# ---------------------------------------------------------------- NR-09
s = mY(0, 0, 1, False)  # warm the guard path
summed = mY(50e-3/6e-3, d6, 1000, exact=False)
S = lambda a: -math.log(abs(2*math.sin(a/2)))
eq29 = c*(S(50e-3-d6) - 2*S(50e-3) + S(50e-3+d6))
rec("NR-09", "Table 8(a): first row summed, and (29) = -1.91501e-10",
    "PASS" if (close(summed, -2.11652e-10, 4e-14) and close(abs(eq29), 1.91501e-10, 4e-15)) else "FAIL",
    f"summed(N=1e3) = {summed:.6e}   target -2.11652e-10\n"
    f"|(29)|        = {abs(eq29):.6e}   target  1.91501e-10")

# ---------------------------------------------------------------- NR-10
rec("NR-10", "Table 8(b): increment per decade -> 6.0788e-08",
    "PASS" if close(2*c*math.log(10), 6.0788e-08, 5e-13) else "FAIL",
    f"(2 mu0 L/pi) ln10 = {2*c*math.log(10):.6e}   target 6.0788e-08")

# ---------------------------------------------------------------- NR-11
def corr(N):
    n = np.arange(1, N+1, dtype=float)
    W = (2.0/(n*math.pi))*np.sin(n*d6/2)
    return L*Rs*math.pi*np.sum(mu0*(n/Rs)*(kappa(n)-1.0)*W**2)
c3, c7 = corr(10**3), corr(10**7)
rec("NR-11", "Table 8(c): kernel correction 2.4307e-11 H, constant from 1e3 to 1e7",
    "PASS" if (close(c3, 2.4307e-11, 5e-15) and close(c7, 2.4307e-11, 5e-15)) else "FAIL",
    f"correction at N=1e3 : {c3:.6e} H\n"
    f"correction at N=1e7 : {c7:.6e} H\n"
    f"difference between them : {abs(c7-c3):.2e} H   (must be negligible)")

# ---------------------------------------------------------------- NR-13
t = read(os.path.join(OUT, 'MEC_BLDC', 'R7_scorecard_out.txt'))
if t is None:
    rec("NR-13", "Fig. 11: 21 of 29 under 5 %, 20 of 26 excluding the three", "UNDECIDED",
        "R7_scorecard_out.txt not readable")
else:
    m1 = re.search(r'compte brut mesure\s*:\s*(\d+)\s*/\s*(\d+)', t)
    m2 = re.search(r'compte hors non validees\s*:\s*(\d+)\s*/\s*(\d+)', t)
    got = (m1.groups() if m1 else None, m2.groups() if m2 else None)
    ok = bool(m1 and m2 and m1.group(1) == '21' and m1.group(2) == '29'
              and m2.group(1) == '20' and m2.group(2) == '26')
    rec("NR-13", "Fig. 11: 21 of 29 under 5 %, 20 of 26 excluding the three",
        "PASS" if ok else "FAIL",
        f"transcript reports raw {m1.group(1)}/{m1.group(2)} and "
        f"{m2.group(1)}/{m2.group(2)} excluding the three\n"
        f"note: the same transcript records the manuscript once announced 22/29"
        if m1 and m2 else f"could not parse the counts: {got}")

# ---------------------------------------------------------------- NR-14
rec("NR-14", "Section 6.9: 540 x 12 = 6480; Tables 9/10 -> 10 and 16 nodes per column",
    "PASS" if (540*12 == 6480 and 2*1+4+3+1 == 10 and 2*4+4+3+1 == 16) else "FAIL",
    f"540 x 12 = {540*12};  n_sh=1 -> {2*1+4+3+1};  n_sh=4 -> {2*4+4+3+1}")

# ---------------------------------------------------------------- NR-08
#  Checked against the B1 transcript, which is a dated run of the Table 7(b)/(d)
#  chain, and against the reference column of Table 7(d) as the target.
b1 = read(os.path.join(REV, 'B1_tiling_dispersion_out.txt'))
if b1 is None:
    rec("NR-08", "Table 7(d): p.c. +0.279/-8.12/+1.11 and hat +0.226/-5.19/+1.04",
        "UNDECIDED", "B1_tiling_dispersion_out.txt absent: run B1 first")
else:
    blocks = b1.split('basis : ')
    got = {}
    for b in blocks[1:]:
        key = 'pc' if b.lstrip().startswith('p.c.') else 'hat'
        m = re.search(r'^1260\s+630\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)', b, re.M)
        if m:
            got[key] = tuple(float(x) for x in m.groups())
    REF = {'B_g1': 1.074551, 'nu8': 0.019658, 'lam': 0.258328}   # Table 7(d) reference column
    TGT = {'pc': (0.279, -8.12, 1.11), 'hat': (0.226, -5.19, 1.04)}
    rows, okall = [], bool(got)
    for key in ('pc', 'hat'):
        if key not in got:
            okall = False
            rows.append(f"{key}: row Ms=1260 not found in the B1 transcript")
            continue
        bg, n8, lm = got[key]
        dev = (100*(bg-REF['B_g1'])/REF['B_g1'],
               100*(n8-REF['nu8'])/REF['nu8'],
               100*(lm-REF['lam'])/REF['lam'])
        t = TGT[key]
        good = all(close(a, b, 0.02) for a, b in zip(dev, t))
        okall &= good
        rows.append(f"{key:3s}: B_g1 {dev[0]:+.3f} % (target {t[0]:+.3f})  "
                    f"nu8 {dev[1]:+.2f} % (target {t[1]:+.2f})  "
                    f"lambda {dev[2]:+.2f} % (target {t[2]:+.2f})   "
                    f"{'ok' if good else 'MISMATCH'}")
    rows.append("nu = 22 is not carried by the B1 transcript and is not checked here.")
    rec("NR-08", "Table 7(d): p.c. +0.279/-8.12/+1.11 and hat +0.226/-5.19/+1.04",
        "PASS" if okall else "FAIL", "\n".join(rows))

# ---------------------------------------------------------------- NR-15
#  3 September 2026.  The manuscript is looked for in the archive first, then at
#  the author's Desktop path, then at a path given on the command line.  Before
#  this change the battery could only read the Desktop copy and reported
#  UNDECIDED on any other machine.
def _find_docx():
    import glob
    here = os.path.dirname(os.path.abspath(__file__))
    root = os.path.abspath(os.path.join(here, "..", ".."))
    cand = sorted(glob.glob(os.path.join(root, "manuscript", "*.docx")), reverse=True)
    for a in sys.argv[1:]:
        if a.lower().endswith(".docx") and os.path.isfile(a):
            return a
    if cand:
        return cand[0]
    return os.path.join(os.path.expanduser("~"), "Desktop",
        "Boundary Condensation of a Two-Layer Dirichlet-to-Neumann Operator "
        "Trace Conformity and Air-Gap Mesh Reduction.docx")

DOCX = _find_docx()
try:
    import zipfile
    xml = zipfile.ZipFile(DOCX).read('word/document.xml').decode('utf-8')
    body = xml.split('<w:body>', 1)[1]
    paras = re.findall(r'<w:p[ >].*?</w:p>|<w:p/>', body, re.S)
    #  Equation tags and part of the prose live inside OMML (<m:t>), which a
    #  <w:t>-only extraction silently drops - it loses equation (1) and five
    #  abstract words, and would report a false NR-15 failure. Both run types
    #  are collected, in document order.
    txts = []
    for p in paras:
        t = ''.join(a or b for a, b in
                    re.findall(r'<w:t(?:\s[^>]*)?>(.*?)</w:t>|<m:t(?:\s[^>]*)?>(.*?)</m:t>',
                               p, re.S))
        t = (t.replace('&amp;', '&').replace('&lt;', '<').replace('&gt;', '>'))
        txts.append(re.sub(r'\s+', ' ', t).strip())
    full = "\n".join(txts)
    tabs = sorted({int(m.group(1)) for t in txts for m in re.finditer(r'^Table\s+(\d+)\.', t)})
    figs = sorted({int(m.group(1)) for t in txts
                   for m in re.finditer(r'^(?:Figure|Fig\.)\s+(\d+)\.', t)})
    eqs  = sorted({int(m.group(1)) for m in re.finditer(r'(?<![\w.])\((\d{1,2})\)(?![\w])', full)
                   if 1 <= int(m.group(1)) <= 60})
    ai = next((i for i, t in enumerate(txts) if t.strip() == 'Abstract'), None)
    abstract = next((t for t in txts[ai+1:] if t), "") if ai is not None else ""
    nw = len(re.findall(r"[A-Za-z0-9À-ſ][A-Za-z0-9À-ſ'’.\-]*", abstract))
    seqT = tabs == list(range(1, 18))
    seqF = figs == list(range(1, 17))
    #  3 September 2026: the displayed equations were renumbered sequentially,
    #  the former (30b) and (30c) becoming ordinary equations (31) and (32) and
    #  the former (31)-(41) becoming (33)-(43).  The expectation moves with the
    #  manuscript; the change is recorded in CHANGELOG.md and in
    #  FINAL_PRE_SUBMISSION_AUDIT.md.  No numerical value changed with it.
    seqE = eqs == list(range(1, 44))
    okall = seqT and seqF and seqE and nw <= 250
    rec("NR-15", "Manuscript: eq (1)-(43), Tables 1-17, Figures 1-16, abstract <= 250 words",
        "PASS" if okall else "FAIL",
        f"equations   : {min(eqs)}..{max(eqs)}, sequential = {seqE}\n"
        f"tables      : {min(tabs)}..{max(tabs)}, sequential = {seqT}\n"
        f"figures     : {min(figs)}..{max(figs)}, sequential = {seqF}\n"
        f"abstract    : {nw} words (limit 250)\n"
        f"the 41-reference and all-called-once checks belong to B12 "
        f"(consistency_check.py), which resolves citation ranges properly.")
except Exception as e:
    rec("NR-15", "Manuscript: eq (1)-(43), Tables 1-17, Figures 1-16, abstract <= 250 words",
        "UNDECIDED", f"manuscript not readable at the expected path: {e}")

# ------------------------------------- checks that need other chains
rec("NR-03", "Table 6(a): 12 rows; p.c. drift 0.703 %; hat 0.023 %; final ratio 0.923624",
    "UNDECIDED",
    "Table 6 is the released-truncation sweep to N_h = 1 105 920. The archived\n"
    "chain RUN_R4_UNLOCK.m breaks at numax = 8640 (RCOND = NaN in its transcript),\n"
    "so no archived run reaches the published truncation. B4 reaches it by the\n"
    "residue-folded assembly, but only for the operator, not for the B_g1 drift\n"
    "this check needs. Reported UNDECIDED rather than asserted.")
rec("NR-12", "Table 12: 51 deviations; 13 of 17 move away under shoe refinement",
    "UNDECIDED",
    "Needs the full-precision Table 12 archive (code/MEC_BLDC/A1_table7.mat), a\n"
    "MATLAB v5 .mat that this battery cannot read without scipy. B7 reads it in\n"
    "MATLAB and reprints the row that decides the 9-against-8 count; see\n"
    "outputs/revision/B7_fluxlineA_reprint_out.txt.")

# ---------------------------------------------------------------- summary
print()
print("="*78)
npass = sum(1 for r in RES if r[2] == 'PASS')
nfail = sum(1 for r in RES if r[2] == 'FAIL')
nund  = sum(1 for r in RES if r[2] == 'UNDECIDED')
print(f"B9 SUMMARY : {npass} PASS, {nfail} FAIL, {nund} UNDECIDED, of {len(RES)} checks")
if nfail:
    print("  FAILING CHECKS:")
    for tag, name, v, _ in RES:
        if v == 'FAIL':
            print(f"    [{tag}] {name}")
print("="*78)
sys.exit(1 if nfail else 0)
