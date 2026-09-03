# -*- coding: utf-8 -*-
"""Values as printed in the submitted manuscript.

This file is a transcription, nothing more.  It exists so that the checks in
verify_all.py compare a recomputation against the printed page rather than
against themselves.  Every entry carries the table and column it was read
from.  No value here has been adjusted.
"""

# ---- Table 5 : hat basis, residual tail beyond N, d = 6 mrad
TABLE5 = [           # N,        residual tail (H),  tail x N^2
    (10**2, 3.4267e-08, 0.000343),
    (10**3, 8.6824e-10, 0.000868),
    (10**4, 1.0837e-11, 0.001084),
    (10**5, 1.1002e-13, 0.001100),
    (10**6, 1.0998e-15, 0.001100),
]

# ---- Table 6(a) : truncation released from the tiling, M_s = 1080
#      N_h, bracket, B_g1 p.c. (T), increment p.c., ratio p.c., B_g1 hat (T), ratio hat
TABLE6A = [
    (    540, 1.721944, 1.078648295, None,       None,     1.078007335, None),
    (   1080, 2.415091, 1.079783044, 1.135e-3,   None,     1.078217653, None),
    (   2160, 3.108239, 1.080846644, 1.064e-3,   0.937299, 1.078246080, 0.135161),
    (   4320, 3.801386, 1.081745009, 8.984e-4,   0.844646, 1.078253756, 0.270048),
    (   8640, 4.494533, 1.082528484, 7.835e-4,   0.872112, 1.078255723, 0.256244),
    (  17280, 5.187680, 1.083223960, 6.955e-4,   0.887681, 1.078256218, 0.251692),
    (  34560, 5.880827, 1.083848282, 6.243e-4,   0.897691, 1.078256342, 0.250433),
    (  69120, 6.573974, 1.084413314, 5.650e-4,   0.905032, 1.078256373, 0.250109),
    ( 138240, 7.267122, 1.084927989, 5.147e-4,   0.910878, 1.078256381, 0.250027),
    ( 276480, 7.960269, 1.085399306, 4.713e-4,   0.915758, 1.078256383, 0.250007),
    ( 552960, 8.653416, 1.085832896, 4.336e-4,   0.919951, 1.078256384, 0.250002),
    (1105920, 9.346563, 1.086233369, 4.005e-4,   0.923624, 1.078256384, 0.249999),
]

# ---- Table 6(b) : self-energy of a column, both bases, M_s = 1080
#      N_h, -Y(i,i) p.c. (H), increment p.c. (H), dev (%), -Y(i,i) hat (H),
#      increment hat x N^2 (H), dev (%)
TABLE6B = [
    (    540, 4.358622e-08, None,         None,        3.102211e-08, None,        None),
    (   1080, 6.437687e-08, 2.079065e-08, +13.615765,  3.575272e-08, 1.379445e-03, +57.202804),
    (   2160, 8.224183e-08, 1.786496e-08,  -2.372371,  3.638082e-08, 7.326222e-04, -16.509733),
    (   4320, 1.004208e-07, 1.817902e-08,  -0.656150,  3.655899e-08, 8.312588e-04,  -5.269024),
    (   8640, 1.186889e-07, 1.826810e-08,  -0.169329,  3.660534e-08, 8.648971e-04,  -1.435567),
    (  17280, 1.369802e-07, 1.829127e-08,  -0.042699,  3.661705e-08, 8.742661e-04,  -0.367865),
    (  34560, 1.552773e-07, 1.829713e-08,  -0.010698,  3.661998e-08, 8.766819e-04,  -0.092561),
    (  69120, 1.735759e-07, 1.829860e-08,  -0.002676,  3.662072e-08, 8.772908e-04,  -0.023178),
    ( 138240, 1.918749e-07, 1.829896e-08,  -0.000669,  3.662090e-08, 8.774433e-04,  -0.005797),
    ( 276480, 2.101740e-07, 1.829905e-08,  -0.000167,  3.662095e-08, 8.774814e-04,  -0.001449),
    ( 552960, 2.284730e-07, 1.829908e-08,  -0.000042,  3.662096e-08, 8.774910e-04,  -0.000362),
    (1105920, 2.467721e-07, 1.829908e-08,  -0.000010,  3.662096e-08, 8.774933e-04,  -0.000091),
]
TABLE6B_PRED_PC  = 1.829909e-08     # note (a): (2 mu0 L / pi) ln 2
TABLE6B_PRED_HAT = 8.774941e-04     # note (b): (3/4) 3 mu0 L / (pi d^2)

# ---- Table 7(a) : the divergent bracket under the locked chain
TABLE7A = [(540, 1.721940), (1080, 1.721944), (2160, 1.721945),
           (4320, 1.721945), (8640, 1.721946)]
TABLE7A_LNPI_GAMMA = 1.721946       # ln pi + gamma, as printed in the caption
TABLE7A_SELF_TERM  = -0.21          # movement of the self term, per cent
TABLE7A_BRACKET_PRED = +0.0003      # movement a stationary bracket predicts, per cent
TABLE7A_RATIO      = 644            # "644 times smaller"

# ---- Table 7(b) : the bore in both bases under the locked chain, six tilings
TABLE7B = {
    "Bg1 p.c.":  [1.07390, 1.07865, 1.07755, 1.07793, 1.07759, 1.07743],
    "Bg1 hat":   [1.07267, 1.07801, 1.07698, 1.07759, 1.07741, 1.07734],
    "nu8 p.c.":  [0.021761, 0.016944, 0.018062, 0.017672, 0.018014, 0.018180],
    "nu8 hat":   [0.023000, 0.017596, 0.018637, 0.018020, 0.018195, 0.018271],
}
TABLE7B_DISP = {"Bg1 p.c.": 0.440, "Bg1 hat": 0.496, "nu8 p.c.": 26.126, "nu8 hat": 28.514}

# Table 7(b) continued : the four integral quantities, six tilings, two bases
TABLE7B_INT = {
    "lambda peak p.c.": [0.259990, 0.261587, 0.261206, 0.261317, 0.261189, 0.261128],
    "L_a p.c.":         [50.07544, 50.97078, 50.77264, 50.85949, 50.80404, 50.77613],
    "M p.c.":           [-2.18502, -2.23023, -2.22015, -2.22458, -2.22177, -2.22035],
    "L_d p.c.":         [52.26046, 53.20100, 52.99278, 53.08407, 53.02581, 52.99648],
    "lambda peak hat":  [0.259531, 0.261379, 0.261018, 0.261206, 0.261132, 0.261099],
    "L_a hat":          [49.83764, 50.83138, 50.65639, 50.79041, 50.76983, 50.75922],
    "M hat":            [-2.17334, -2.22313, -2.21426, -2.22107, -2.22003, -2.21949],
    "L_d hat":          [52.01097, 53.05450, 52.87065, 53.01148, 52.98986, 52.97872],
}
TABLE7B_INT_DISP6 = {"lambda peak p.c.": 0.612, "L_a p.c.": 1.766, "M p.c.": 2.039,
                     "L_d p.c.": 1.777, "lambda peak hat": 0.708, "L_a hat": 1.964,
                     "M hat": 2.251, "L_d hat": 1.976}
TABLE7B_INT_DISP5 = {"lambda peak p.c.": 0.176, "L_a p.c.": 0.390, "M p.c.": 0.453,
                     "L_d p.c.": 0.392, "lambda peak hat": 0.138, "L_a hat": 0.345,
                     "M hat": 0.400, "L_d hat": 0.347}

# ---- Table 7(d) : both bases against the reference, M_s = 1260
TABLE7D = [   # quantity, reference, p.c., dev p.c. (%), hat, dev hat (%)
    ("B_g1 (T)",       1.074551, 1.077547, +0.279, 1.076981, +0.226),
    ("nu = 8 (T)",     0.019658, 0.018062,  -8.12, 0.018637,  -5.19),
    ("nu = 22 (T)",    0.032760, 0.035238,  +7.56, 0.036196, +10.49),
    ("lambda peak (Wb)", 0.258328, 0.261206, +1.11, 0.261018, +1.04),
]

# ---- Table 8(a) : off diagonal, d_i = d_j = 6 mrad, Delta theta = 50 mrad
TABLE8A = [   # N, summed (H), closed form (29) (H), relative deviation
    (10**3, -2.1165e-10, -1.9150e-10, 1.052e-1),
    (10**4, -2.4041e-10, -1.9150e-10, 2.554e-1),
    (10**5, -2.0196e-10, -1.9150e-10, 5.460e-2),
    (10**6, -1.9154e-10, -1.9150e-10, 2.260e-4),
    (10**7, -1.9150e-10, -1.9150e-10, 2.300e-5),
]

# ---- Table 8(b) : diagonal, asymptotic kernel, d = 6 mrad
TABLE8B = [   # N, summed (H), closed form (30) (H), rel dev, increment per decade (H)
    (10**3, 6.4338e-08, 6.2541e-08, 2.874e-02, None),
    (10**4, 1.2346e-07, 1.2333e-07, 1.051e-03, 5.9121e-08),
    (10**5, 1.8412e-07, 1.8412e-07, 9.523e-06, 6.0657e-08),
    (10**6, 2.4491e-07, 2.4491e-07, 7.692e-06, 6.0792e-08),
    (10**7, 3.0569e-07, 3.0569e-07, 1.373e-06, 6.0786e-08),
]

# ---- Table 8(c) : diagonal, exact two-layer kernel, d = 6 mrad
TABLE8C = [   # N, asymptotic (H), exact (H), increment per decade, exact (H)
    (10**4, 1.2346e-07, 1.2348e-07, 5.9121e-08),
    (10**5, 1.8412e-07, 1.8414e-07, 6.0657e-08),
    (10**6, 2.4491e-07, 2.4493e-07, 6.0792e-08),
    (10**7, 3.0569e-07, 3.0572e-07, 6.0786e-08),
]
TABLE8C_SHIFT = 2.4307e-11          # caption: constant offset between the two kernels

# ---- Table 12 : the seventeen compared quantities, three network columns
#      name, mesh n_sh=1, mesh n_sh=2, lumped, FEA, published deviations
TABLE12 = [
    ("B_g1 (T)",                1.07787, 1.07908, 1.07755, 1.07455, (+0.31, +0.42, +0.28)),
    ("|B_r| mean (T)",          0.76182, 0.76273, 0.76147, 0.76091, (+0.12, +0.24, +0.07)),
    ("B_r peak (T)",            0.98051, 0.97536, 0.98600, 0.99073, (-1.03, -1.55, -0.48)),
    ("B_t rms (T)",             0.13755, 0.13478, 0.14317, 0.13326, (+3.22, +1.14, +7.43)),
    ("slot sideband nu=8 (T)",  0.01704, 0.01581, 0.01806, 0.01966, (-13.34, -19.58, -8.12)),
    ("flux line A (Wb/m)",     0.005861, 0.005868, 0.005860, 0.005855, (+0.11, +0.22, +0.09)),
    ("flux linkage, peak (Wb)", 0.25894, 0.25908, 0.26121, 0.25833, (+0.24, +0.29, +1.11)),
    ("MMF per magnet (A)",      761.674, 759.329, 762.912, 769.445, (-1.01, -1.31, -0.85)),
    ("MMF per gap (A)",         711.386, 711.880, 725.064, 721.504, (-1.40, -1.33, +0.49)),
    ("MMF dispersion (%)",        6.369,   5.959,   6.810,   6.585, (-3.27, -9.50, +3.42)),
    ("phase back-EMF, peak (V)",238.813, 238.818, 240.877, 233.031, (+2.48, +2.48, +3.37)),
    ("six-step envelope (V)",   464.556, 464.645, 468.719, 459.918, (+1.01, +1.03, +1.91)),
    ("L_a (mH)",                 50.379,  50.484,  50.773,  50.209, (+0.34, +0.55, +1.12)),
    ("M (mH)",                   -2.016,  -2.021,  -2.220,  -2.133, (-5.48, -5.24, +4.09)),
    ("L_d (mH)",                 52.395,  52.505,  52.993,  52.342, (+0.10, +0.31, +1.24)),
    ("k_r",                     1.07069, 1.06665, 1.05220, 1.08693, (-1.49, -1.87, -3.20)),
    ("k_l",                     0.85326, 0.85368, 0.85371, 0.86723, (-1.61, -1.56, -1.56)),
]

# ---- Section 7.1 : Carter's factor, exact and approximate forms
CARTER = dict(b0_mm=2.000, g_mm=1.000, hm_mm=3.500, mu_r=1.0390,
              D_mm=69.356, Ns=15,
              gprime_mm=4.3686, tau_s_mm=14.5259, x=0.22891,
              gamma_exact=0.033072, kC_exact=1.010046,
              gamma_approx=0.038402, kC_approx=1.011684, difference_pct=0.16)

# ---- Section 6.6 : skin depth in the magnet
SKIN = dict(Ns=15, speed_rpm=1688, f_Hz=422.0, rho_ohm_m=1.80e-6,
            sigma_S_per_m=555556, mu_r=1.0390,
            delta_mm=32.25, hm_mm=3.5, ratio=9.2)

# ---- Sections 6.9 and 7.1 : cost figures
COST = dict(sweep_total_s=0.119, single_position_s=0.111, n_positions=128,
            reassembled_sweep_s=13.79, marginal_reassembled_ms=108,
            marginal_backsub_ms=0.06, factor_sweep=115,
            sweep_pct_of_independent=0.84,
            sliding_interface_s=0.714, invariant_operator_s=0.056, factor_interface=12.7,
            positions_721=721, sweep_721_s=3.0, per_position_ms=4.16)
