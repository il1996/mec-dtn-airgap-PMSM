# -*- coding: utf-8 -*-
"""Two-layer bore operator, rebuilt from the published equations alone.

Nothing here is imported from the MATLAB chain.  Every quantity is formed from
the equations of the manuscript and from the geometry of Table 3, so that a
reader can check the published closed forms and tables without the original
chain and without MATLAB.

Equations used, by their number in the submitted manuscript:

    (10)  alpha_n = 1 / [ cosh(n Ua) + mu_r sinh(n Ua) coth(n Um) ]
    (13)  g_n     = - mu0 (n/Rs) (cosh(n Ua) - alpha_n) / sinh(n Ua)
          ghat_n  = pi L Rs g_n                                   (Section 3.1)
    (14)  W^c, W^s : projection of a column indicator
    (26)  W_n(i)  = (2 / n pi) sin(n d / 2) e^{-i n theta_i}
    (29)  off-diagonal closed form, four arguments
    (30)  diagonal closed form
    (31)  first-neighbour closed form
    (32)  -Y_N = c ln N Delta_h + O(1)
    (33)  hat coefficients W_n = 4 sin^2(n d / 2) / (pi n^2 d)

Geometry, Table 3 of the manuscript:
    stator bore diameter D = 69.356 mm      -> Rs  = 34.678 mm
    mechanical air gap    g = 1.000 mm      -> Rro = 33.678 mm
    magnet thickness     hm = 3.500 mm      -> rmi = 30.178 mm
    recoil permeability  mu_r = 1.0390
    active length         L = 33.0 mm
"""

import math

MU0 = 4.0 * math.pi * 1e-7
GAMMA = 0.5772156649015328606          # Euler-Mascheroni

# ---------------------------------------------------------------- geometry
D_BORE_MM = 69.356
GAP_MM    = 1.000
HM_MM     = 3.500
MU_R      = 1.0390
L_M       = 0.033

RS  = D_BORE_MM / 2.0e3
RRO = RS - GAP_MM / 1.0e3
RMI = RRO - HM_MM / 1.0e3
UA  = math.log(RS / RRO)
UM  = math.log(RRO / RMI)

C_CONST = MU0 * L_M / math.pi           # c of Proposition 1


def kappa(n, exact=True):
    """(cosh(n Ua) - alpha_n) / sinh(n Ua).  Tends to 1; = 1 for the
    asymptotic kernel."""
    if not exact:
        return 1.0
    x = n * UA
    if x > 350.0:                        # cosh/sinh have already degenerated
        return 1.0
    y = n * UM
    coth = 1.0 / math.tanh(y) if y < 350.0 else 1.0
    dn = math.cosh(x) + MU_R * math.sinh(x) * coth
    return (math.cosh(x) - 1.0 / dn) / math.sinh(x)


def ghat(n, exact=True):
    """Assembled bore admittance of (13) times pi L Rs, in H.  Negative."""
    return -math.pi * L_M * MU0 * n * kappa(n, exact)


# ------------------------------------------------------- condensed entries
def self_energy_pc(N, d, exact=True):
    """-Y(i,i) under the piecewise-constant basis, by direct summation of (27)."""
    s = 0.0
    for n in range(1, N + 1):
        s += kappa(n, exact) / n * math.sin(n * d / 2.0) ** 2
    return 4.0 * MU0 * L_M / math.pi * s


def self_energy_hat(N, d, exact=True):
    """-Y(i,i) under the hat basis, by direct summation with (33)."""
    s = 0.0
    for n in range(1, N + 1):
        s += kappa(n, exact) / n ** 3 * math.sin(n * d / 2.0) ** 4
    return 16.0 * MU0 * L_M / (math.pi * d * d) * s


def offdiag_pc(N, di, dj, dtheta, exact=True):
    """-Y(i,j) by direct summation, the form the proof of Proposition 1 starts from."""
    s = 0.0
    for n in range(1, N + 1):
        s += (1.0 / n) * math.sin(n * di / 2.0) * math.sin(n * dj / 2.0) * math.cos(n * dtheta)
    return 4.0 * MU0 * L_M / math.pi * s


# ------------------------------------------------------------ closed forms
def offdiag_closed(di, dj, dtheta):
    """(29): the four-argument closed form, valid for |i-j| >= 2."""
    a1 = 0.5 * (di - dj) + dtheta
    a2 = 0.5 * (di - dj) - dtheta
    a3 = 0.5 * (di + dj) + dtheta
    a4 = 0.5 * (di + dj) - dtheta
    lg = lambda a: math.log(abs(2.0 * math.sin(a / 2.0)))
    return C_CONST * (-lg(a1) - lg(a2) + lg(a3) + lg(a4))


def bracket(N, d):
    """ln N + gamma + ln|2 sin(d/2)|, the divergent bracket of (30)."""
    return math.log(N) + GAMMA + math.log(abs(2.0 * math.sin(d / 2.0)))


def diag_closed(N, d):
    """(30): -Y(i,i) = (2 mu0 L / pi) [ ln N + gamma + ln|2 sin(d/2)| ] + o(1)."""
    return 2.0 * MU0 * L_M / math.pi * bracket(N, d)


def neighbour_closed(N, d):
    """(31): the first-neighbour entry on a contiguous tiling."""
    return (-C_CONST * (math.log(N) + GAMMA)
            + C_CONST * (math.log(abs(2.0 * math.sin(d)))
                         - 2.0 * math.log(abs(2.0 * math.sin(d / 2.0)))))


def hat_tail(N, d, exact=False, reach=400):
    """Residual tail of the hat energy series beyond N (Table 5)."""
    s = 0.0
    for n in range(N + 1, N * reach):
        s += kappa(n, exact) / n ** 3 * math.sin(n * d / 2.0) ** 4
    return 16.0 * MU0 * L_M / (math.pi * d * d) * s


# ------------------------------------------------- parameter-free constants
def slope_per_decade():
    """(2 mu0 L / pi) ln 10, the increment of the diagonal per decade."""
    return 2.0 * MU0 * L_M / math.pi * math.log(10.0)


def slope_per_doubling():
    """(2 mu0 L / pi) ln 2, note (a) of Table 6."""
    return 2.0 * MU0 * L_M / math.pi * math.log(2.0)


def hat_tail_constant(d):
    """3 mu0 L / (pi d^2): the limit of (residual tail) x N^2, Table 5."""
    return 3.0 * MU0 * L_M / (math.pi * d * d)


def hat_increment_constant(d):
    """(3/4) 3 mu0 L / (pi d^2): the limit of (increment N -> 2N) x N^2,
    note (b) of Table 6."""
    return 0.75 * hat_tail_constant(d)


if __name__ == "__main__":
    print(f"Rs  = {RS*1e3:.4f} mm     Rro = {RRO*1e3:.4f} mm     rmi = {RMI*1e3:.4f} mm")
    print(f"Ua  = {UA:.6f}          Um  = {UM:.6f}")
    print(f"c   = mu0 L / pi        = {C_CONST:.6e} H")
    print(f"(2 mu0 L / pi) ln 10    = {slope_per_decade():.6e} H")
    print(f"(2 mu0 L / pi) ln 2     = {slope_per_doubling():.6e} H")
