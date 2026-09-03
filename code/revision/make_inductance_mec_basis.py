# -*- coding: utf-8 -*-
"""Derive inductance_mec_basis.m from the published inductance_mec.m by four
   surgical substitutions. Byte-preserving (latin-1). Published file is READ ONLY."""
import os, difflib
R = r"C:\Users\hp\Desktop\mec-dtn-airgap-pmsm"
src = os.path.join(R,"code","MEC_BLDC","inductance_mec.m")
dstdir = os.path.join(R,"code","revision")
os.makedirs(dstdir, exist_ok=True); os.makedirs(os.path.join(R,"outputs","revision"), exist_ok=True)
dst = os.path.join(dstdir,"inductance_mec_basis.m")

raw = open(src,'rb').read()
EOL = "\r\n" if b"\r\n" in raw else "\n"
print("detected line ending:", repr(EOL))
orig = raw.decode('latin-1')
lines = orig.split(EOL)

BANNER = [
"%INDUCTANCE_MEC_BASIS  Basis-aware copy of inductance_mec.m (task B1).",
"%",
"%   DERIVED FILE - DO NOT EDIT BY HAND. Generated from",
"%   code/MEC_BLDC/inductance_mec.m by code/revision/make_inductance_mec_basis.py.",
"%   The published original is not modified. The ONLY differences are:",
"%     (1) the function name, and a 5th argument 'basis' defaulting to 'p0';",
"%     (2) that basis is forwarded to airgap_magnet;",
"%     (3) R.basis is recorded on the returned struct.",
"%   With basis='p0' this reproduces inductance_mec exactly; the calling script",
"%   proves that numerically, at every tiling, before using it.",
]

SUBS = [
 ("function R = inductance_mec(M, Nsurf, muI, kfringe)",
  "function R = inductance_mec_basis(M, Nsurf, muI, kfringe, basis)"),
 ("AG=airgap_magnet(M,ths,dths,numax);",
  "AG=airgap_magnet(M,ths,dths,numax,basis);"),
 ("R.Nsurf=Nsurf; R.muI=muI; R.kfringe=kfringe;",
  "R.Nsurf=Nsurf; R.muI=muI; R.kfringe=kfringe; R.basis=basis;"),
]
out=[]
for i,l in enumerate(lines):
    rep=l
    for a,b in SUBS:
        if l.strip()==a:
            rep=l.replace(a,b); break
    out.append(rep)
    if i==0:                                   # banner directly under the function line
        out.extend(BANNER)
        out.append("if nargin<5||isempty(basis), basis='p0'; end")

# assert each substitution fired exactly once
joined=EOL.join(out)
for a,b in SUBS:
    assert joined.count(b)==1, "substitution failed: "+b
    assert joined.count(a)==0 or a in b, "original still present: "+a
assert joined.count("if nargin<5||isempty(basis), basis='p0'; end")==1

open(dst,'wb').write(joined.encode('latin-1'))
print("wrote", dst)
print("\n=== unified diff, published original -> derived copy ===")
for line in difflib.unified_diff(lines, out, fromfile="code/MEC_BLDC/inductance_mec.m",
                                 tofile="code/revision/inductance_mec_basis.m", lineterm='', n=2):
    print(line)
