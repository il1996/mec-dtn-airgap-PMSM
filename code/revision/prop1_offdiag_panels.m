% prop1_offdiag_panels.m
% -------------------------------------------------------------------------
% TASK B2 - Table 8(a): the case where (29) fails, and the operational
%           truncation.
%
% Panel A  (non-regression) : published Table 8(a), d_i = d_j = 6 mrad,
%                             dtheta = 50 mrad, asymptotic kernel, N = 1e3..1e7.
% Panel I  (new)            : the FIRST-NEIGHBOUR case dtheta = d, d = 6 mrad,
%                             where (29) fails because a_3 = (d_i+d_j)/2 - dtheta
%                             vanishes and one of the four series becomes the
%                             harmonic series.  Compared against (30b).
% Panel II (new)            : the OPERATIONAL truncation.  EXACT two-layer
%                             kernel, Ms = 1080, d = 2*pi/1080, k = 1,2,3,5,10,
%                             N in {540, 1080, 1e4, 1e5, 1e6}; deviation of each
%                             entry from its N -> infinity limit.  k = 1 has no
%                             limit and is reported as DIVERGENT.
%
% Author  : revision chain, task B2
% Date    : 2026-09-02
% Outputs : outputs/revision/B2_prop1_offdiag_out.txt
%
% Governing rules honoured: no invented value, no harmonisation, no tuned
% parameter, no wall-clock benchmark; writes only under code/revision and
% outputs/revision.
% -------------------------------------------------------------------------

clear; clc;

RUN_DATE = '2026-09-02';
ROOT     = 'C:\Users\hp\Desktop\mec-dtn-airgap-pmsm';
OUTDIR   = fullfile(ROOT,'outputs','revision');
OUTFILE  = fullfile(OUTDIR,'B2_prop1_offdiag_out.txt');

if ~exist(OUTDIR,'dir'); mkdir(OUTDIR); end

% ---- diary APPENDS: rotate any previous transcript so the file holds one run
if exist(OUTFILE,'file')
    stamp = datestr(now,'yyyymmdd_HHMMSS'); %#ok<TNOW1,DATST>
    movefile(OUTFILE, fullfile(OUTDIR, ...
        sprintf('B2_prop1_offdiag_out.superseded_%s.txt',stamp)));
end
diary(OUTFILE); diary on;

% =========================================================================
% VALIDATED KERNEL  (verbatim from the shared specification)
% =========================================================================
mu0   = 4*pi*1e-7;
L     = 0.033;            % m,   active length
R_s   = 34.678e-3;        % m,   stator bore (reported, not used below)
U_a   = 2.9260675e-2;     % logarithmic thickness, air layer
U_m   = 1.0973162e-1;     % logarithmic thickness, magnet layer
mu_r  = 1.038952;         % magnet recoil permeability
c     = mu0*L/pi;         % = 1.320000e-08
gamma = 0.5772156649015329;

GUARD_UA = 25;            % n*U_a > 25  ->  kappa_n = 1 exactly
GUARD_UM = 25;            % clamp of the coth argument
n_guard  = floor(GUARD_UA/U_a);   % last harmonic that still carries kappa ~= 1

Ms_op = 1080;             % operational tiling
d_op  = 2*pi/Ms_op;       % operational column width
d6    = 6e-3;             % 6 mrad, the width used by published Table 8

% =========================================================================
% CONFIGURATION BLOCK
% =========================================================================
fprintf('=========================================================================\n');
fprintf(' TASK B2 - Proposition 1, off-diagonal entries: where (29) fails, and\n');
fprintf('           the behaviour at the operational truncation\n');
fprintf(' Date (ISO)      : %s\n', RUN_DATE);
fprintf(' Script          : %s\n', fullfile(ROOT,'code','revision','prop1_offdiag_panels.m'));
fprintf(' Transcript      : %s\n', OUTFILE);
fprintf(' MATLAB          : %s (%s)\n', version, computer);
fprintf('-------------------------------------------------------------------------\n');
fprintf(' CONFIGURATION\n');
fprintf('   Basis         : piecewise-constant (unit) on contiguous arcs of the\n');
fprintf('                   stator bore; condensed two-layer annulus operator.\n');
fprintf('   Kernel        : two variants, both used below and always labelled.\n');
fprintf('                   ASYMPTOTIC : kappa_n replaced by 1 for all n.\n');
fprintf('                   EXACT      : alpha_n = 1/(cosh(n Ua) + mu_r sinh(n Ua) coth(n Um))\n');
fprintf('                                kappa_n = (cosh(n Ua) - alpha_n)/sinh(n Ua)\n');
fprintf('   Guard         : n*U_a > %g  ->  kappa_n = 1 exactly (n > %d);\n', GUARD_UA, n_guard);
fprintf('                   coth argument clamped at min(n*U_m,%g).\n', GUARD_UM);
fprintf('   General entry : -Y(i,j) = (4 mu0 L/pi) SUM_{n=1..N} (kappa_n/n)\n');
fprintf('                              * sin(n d_i/2) sin(n d_j/2) cos(n dtheta)\n');
fprintf('                   (with d_i = d_j = d and dtheta = k d this is the\n');
fprintf('                    condensed uniform-tiling form of the specification).\n');
fprintf('   Closed forms  : (29)  -Y = c[ ln|2 sin(a3/2)| + ln|2 sin(a4/2)|\n');
fprintf('                                 - ln|2 sin(a1/2)| - ln|2 sin(a2/2)| ]\n');
fprintf('                   a1 = (di-dj)/2 - dtheta, a2 = (di-dj)/2 + dtheta,\n');
fprintf('                   a3 = (di+dj)/2 - dtheta, a4 = (di+dj)/2 + dtheta.\n');
fprintf('                   (30b) -Y(i,i+-1) = -c(ln N + gamma)\n');
fprintf('                                      + c[ ln|2 sin d| - 2 ln|2 sin(d/2)| ]\n');
fprintf('   Tilings       : d = 6 mrad (published Table 8) and the OPERATIONAL\n');
fprintf('                   uniform tiling Ms = %d, d = 2*pi/Ms = %.6e rad.\n', Ms_op, d_op);
fprintf('   Truncations   : N = 1e3..1e7 (panels A and I); N in {540, 1080, 1e4,\n');
fprintf('                   1e5, 1e6} (panel II, the operational range).\n');
fprintf('   Summation     : direct, vectorised, chunked; no closed form is used\n');
fprintf('                   inside any summed column.\n');
fprintf('   Constants     : mu0 = %.15e   L = %.6f m   R_s = %.6e m\n', mu0, L, R_s);
fprintf('                   U_a = %.7e  U_m = %.7e  mu_r = %.6f\n', U_a, U_m, mu_r);
fprintf('                   c = mu0*L/pi = %.6e   gamma = %.16f\n', c, gamma);
fprintf('   No tuned parameter, no fitted constant, no wall-clock measurement.\n');
fprintf('=========================================================================\n\n');

% =========================================================================
% CONTROL VALUES - printed BEFORE any new result
% =========================================================================
fprintf('=========================================================================\n');
fprintf(' CONTROL-VALUE COMPARISON (must pass before any new result is used)\n');
fprintf('=========================================================================\n');
fprintf('%-58s %16s %16s %12s  %s\n','control','reference','reproduced','deviation','verdict');

ok_all = true;

% ---- C1 : kappa_1
k1 = kappa_vec(1, U_a, U_m, mu_r, GUARD_UA, GUARD_UM);
ok_all = report_rel('kappa_1', 7.4607050, k1, 1e-6, '%16.10f') && ok_all;

% ---- C2 : first n with |kappa_n - 1| < 1e-6
nn   = (1:5000).';
kv   = kappa_vec(nn, U_a, U_m, mu_r, GUARD_UA, GUARD_UM);
nfirst = find(abs(kv-1) < 1e-6, 1, 'first');
ok = (nfirst == 181); ok_all = ok_all && ok;
fprintf('%-58s %16d %16d %12s  %s\n','first n with |kappa_n - 1| < 1e-6', 181, nfirst, '-', verdict(ok));

% ---- C3 : min / max kappa_n over n <= 5000
ok_all = report_rel('min kappa_n, n <= 5000', 1.000000, min(kv), 1e-6, '%16.6f') && ok_all;
ok_all = report_rel('max kappa_n, n <= 5000', 7.460705, max(kv), 1e-6, '%16.6f') && ok_all;

% ---- C4 : -Y(i,i) at Ms = 1080, N = 34560, asymptotic
Yii = offdiag_sum(34560, d_op, d_op, 0.0, false, mu0, L, U_a, U_m, mu_r, GUARD_UA, GUARD_UM, n_guard);
ok_all = report_rel('-Y(i,i)    Ms=1080, N=34560, asymptotic', 1.55254492e-07, Yii, 1e-8, '%16.8e') && ok_all;

% ---- C5 : -Y(i,i+-1) at Ms = 1080, N = 34560, asymptotic
Yi1 = offdiag_sum(34560, d_op, d_op, d_op, false, mu0, L, U_a, U_m, mu_r, GUARD_UA, GUARD_UM, n_guard);
ok_all = report_rel('-Y(i,i+-1) Ms=1080, N=34560, asymptotic', -6.84780037e-08, Yi1, 1e-5, '%16.8e') && ok_all;

% ---- C6 : (30b) relative deviation at N = 34560
b30b = eq30b(34560, d_op, c, gamma);
dev  = (Yi1 - b30b)/abs(b30b);
ok   = abs(dev - (-8.343e-06)) <= abs(3*(-8.343e-06));   % "factor 3" tolerance
ok_all = ok_all && ok;
fprintf('%-58s %16.3e %16.3e %12s  %s\n','(30b) rel. deviation at N=34560', -8.343e-06, dev, '-', verdict(ok));

% ---- C7 : published Table 8(a), the four/five rows must stay identical
fprintf('\n  Non-regression, published Table 8(a): d_i = d_j = 6 mrad, dtheta = 50 mrad,\n');
fprintf('  asymptotic kernel.  Manuscript values as rendered.\n');
fprintf('  %-6s %14s %14s %14s %14s %12s %12s  %s\n', ...
        'N','summed(ms)','summed(here)','(29)(ms)','(29)(here)','reldev(ms)','reldev(here)','verdict');
ms_sum  = [-2.1165e-10, -2.4041e-10, -2.0196e-10, -1.9154e-10, -1.9150e-10];
ms_c29  = [-1.9150e-10, -1.9150e-10, -1.9150e-10, -1.9150e-10, -1.9150e-10];
ms_dev  = [ 1.052e-01,   2.554e-01,   5.460e-02,   2.260e-04,   2.300e-05];
c29_a   = eq29(d6, d6, 50e-3, c);
for e = 3:7
    N = 10^e; ii = e-2;
    s = offdiag_sum(N, d6, d6, 50e-3, false, mu0, L, U_a, U_m, mu_r, GUARD_UA, GUARD_UM, n_guard);
    dv = abs(s - c29_a)/abs(c29_a);
    ok = (abs(s - ms_sum(ii))/abs(ms_sum(ii)) < 5e-4) && ...
         (abs(c29_a - ms_c29(ii))/abs(ms_c29(ii)) < 5e-4) && ...
         (abs(dv - ms_dev(ii))/abs(ms_dev(ii)) < 5e-3);
    ok_all = ok_all && ok;
    fprintf('  1e%-4d %14.4e %14.4e %14.4e %14.4e %12.3e %12.3e  %s\n', ...
            e, ms_sum(ii), s, ms_c29(ii), c29_a, ms_dev(ii), dv, verdict(ok));
end

% ---- C8 : increment per decade of -Y(i,i+-1), d = 6 mrad, asymptotic
%           (computed in panel I below; the reference is -c*ln10)
fprintf('\n%-58s %16.4e %16s %12s  %s\n','-c*ln(10) (reference increment per decade)', ...
        -c*log(10), '(see panel I)', '-', 'INFO');

fprintf('\n  ALL CONTROLS ABOVE: %s\n', verdict(ok_all));
fprintf('  (the k=2,3,5,10 deviation controls are checked inside panel II)\n');
fprintf('=========================================================================\n\n');

% =========================================================================
% PANEL A (non-regression, reprinted as a stand-alone table)
% =========================================================================
fprintf('=========================================================================\n');
fprintf(' PANEL A  (published Table 8(a), unchanged)\n');
fprintf(' Off diagonal, d_i = d_j = 6 mrad, dtheta = 50 mrad, ASYMPTOTIC kernel.\n');
fprintf(' dtheta / d = %.4f columns -- this separation is not a matrix entry of\n', 50e-3/d_op);
fprintf(' the operational tiling Ms = 1080 (d = %.3f mrad), and a3 = (di+dj)/2 -\n', d_op*1e3);
fprintf(' dtheta = %.4f rad is far from zero, so (29) is well posed here.\n', d6-50e-3);
fprintf('-------------------------------------------------------------------------\n');
fprintf(' %-8s %16s %16s %14s\n','N','summed','(29)','rel. dev.');
for e = 3:7
    N = 10^e;
    s = offdiag_sum(N, d6, d6, 50e-3, false, mu0, L, U_a, U_m, mu_r, GUARD_UA, GUARD_UM, n_guard);
    fprintf(' 1e%-6d %16.4e %16.4e %14.3e\n', e, s, c29_a, abs(s-c29_a)/abs(c29_a));
end
fprintf('=========================================================================\n\n');

% =========================================================================
% PANEL I : the FIRST-NEIGHBOUR case, where (29) fails
% =========================================================================
fprintf('=========================================================================\n');
fprintf(' PANEL I  (new)  Off diagonal at the FIRST NEIGHBOUR of a contiguous\n');
fprintf('                 tiling: dtheta = (d_i + d_j)/2 = d = 6 mrad.\n');
fprintf(' ASYMPTOTIC kernel, N = 1e3..1e7.\n');
fprintf('-------------------------------------------------------------------------\n');
fprintf(' WHY (29) FAILS HERE.  With d_i = d_j = d and dtheta = d the four\n');
fprintf(' arguments are a1 = -d, a2 = +d, a3 = 0, a4 = 2d.  The third series is\n');
fprintf(' SUM cos(0)/n = SUM 1/n, the harmonic series: (29) evaluates\n');
fprintf(' ln|2 sin(a3/2)| = ln 0 = -Inf and the entry has no finite limit.\n');
fprintf(' (29) at a3 = 0 gives : %s\n', num2str(eq29(d6,d6,d6,c)));
fprintf(' This is not a degenerate corner case: on a contiguous tiling it is the\n');
fprintf(' generic situation of the two largest off-diagonal entries of EVERY row.\n');
fprintf(' The correct statement at the first neighbour is (30b), which carries the\n');
fprintf(' ln N explicitly.  That is what is compared below.\n');
fprintf('-------------------------------------------------------------------------\n');
fprintf(' %-8s %16s %16s %14s %16s\n','N','summed','(30b)','rel. dev.','incr./decade');
prev = NaN; inc_last = NaN;
for e = 3:7
    N = 10^e;
    s = offdiag_sum(N, d6, d6, d6, false, mu0, L, U_a, U_m, mu_r, GUARD_UA, GUARD_UM, n_guard);
    b = eq30b(N, d6, c, gamma);
    if isnan(prev)
        fprintf(' 1e%-6d %16.6e %16.6e %14.3e %16s\n', e, s, b, (s-b)/abs(b), '--');
    else
        inc = s - prev; inc_last = inc;
        fprintf(' 1e%-6d %16.6e %16.6e %14.3e %16.4e\n', e, s, b, (s-b)/abs(b), inc);
    end
    prev = s;
end
fprintf('-------------------------------------------------------------------------\n');
fprintf(' predicted increment per decade  -c*ln10 = %.4e\n', -c*log(10));
fprintf(' last measured increment                 = %.4e\n', inc_last);
ok = abs(inc_last - (-c*log(10)))/abs(c*log(10)) < 1e-2;
fprintf(' agreement within 1%%                     : %s  (rel. dev. %.3e)\n', ...
        verdict(ok), abs(inc_last + c*log(10))/abs(c*log(10)));
fprintf(' READING: the first-neighbour entry DOES NOT CONVERGE.  It decreases by\n');
fprintf(' -c ln 10 per decade of N, exactly like the diagonal, and with the same\n');
fprintf(' coefficient c (not 2c).  Panel (a) of Table 8 tests (29) at a separation\n');
fprintf(' that converges and therefore never sees this.\n');
fprintf('=========================================================================\n\n');

% =========================================================================
% PANEL II : the OPERATIONAL truncation, EXACT kernel
% =========================================================================
fprintf('=========================================================================\n');
fprintf(' PANEL II (new)  The OPERATIONAL truncation.  EXACT two-layer kernel,\n');
fprintf('                 uniform contiguous tiling Ms = %d, d = 2*pi/Ms =\n', Ms_op);
fprintf('                 %.9e rad (%.4f mrad).\n', d_op, d_op*1e3);
fprintf(' Entries k = 1,2,3,5,10 at N in {540, 1080, 1e4, 1e5, 1e6}.\n');
fprintf(' The validated chain runs at N = 540 or 630, i.e. the FIRST column below.\n');
fprintf('-------------------------------------------------------------------------\n');
fprintf(' N -> infinity limit for k >= 2 is taken analytically as\n');
fprintf('     lim = (29) + (4 mu0 L/pi) SUM_{n<=%d} ((kappa_n - 1)/n) sin^2(n d/2) cos(n k d)\n', 5*n_guard);
fprintf(' the correction being finite and fully acquired by n ~ 200 (kappa_n = 1\n');
fprintf(' to 1e-6 from n = %d).  A direct summation at N = 1e7 is printed beside\n', nfirst);
fprintf(' it as an independent cross-check; neither is fitted to the other.\n');
fprintf(' k = 1 has NO limit (panel I): it is reported as DIVERGENT.\n');
fprintf('-------------------------------------------------------------------------\n');

kk_list = [1 2 3 5 10];
NN_list = [540 1080 1e4 1e5 1e6];
N_cross = 1e7;
Ncorr   = 5*n_guard;

lim_an  = nan(size(kk_list));
lim_num = nan(size(kk_list));
vals    = nan(numel(kk_list), numel(NN_list));

for a = 1:numel(kk_list)
    k = kk_list(a);
    for b = 1:numel(NN_list)
        vals(a,b) = offdiag_sum(NN_list(b), d_op, d_op, k*d_op, true, ...
                                mu0, L, U_a, U_m, mu_r, GUARD_UA, GUARD_UM, n_guard);
    end
    if k >= 2
        n  = (1:Ncorr).';
        kp = kappa_vec(n, U_a, U_m, mu_r, GUARD_UA, GUARD_UM);
        corr = (4*mu0*L/pi)*sum((kp-1).*sin(n*d_op/2).^2.*cos(n*k*d_op)./n);
        lim_an(a)  = eq29(d_op, d_op, k*d_op, c) + corr;
        lim_num(a) = offdiag_sum(N_cross, d_op, d_op, k*d_op, true, ...
                                 mu0, L, U_a, U_m, mu_r, GUARD_UA, GUARD_UM, n_guard);
    end
end

fprintf('\n (II.1)  Entry values, EXACT kernel, Ms = %d  [H]\n', Ms_op);
fprintf(' %-4s %14s %14s %14s %14s %14s | %16s\n','k','N=540','N=1080','N=1e4','N=1e5','N=1e6','limit N->inf');
for a = 1:numel(kk_list)
    fprintf(' %-4d', kk_list(a));
    fprintf(' %14.6e', vals(a,:));
    if kk_list(a) == 1
        fprintf(' | %16s\n','DIVERGENT');
    else
        fprintf(' | %16.6e\n', lim_an(a));
    end
end

fprintf('\n (II.2)  Deviation from the N -> infinity limit  [%%]\n');
fprintf(' %-4s %14s %14s %14s %14s %14s | %s\n','k','N=540','N=1080','N=1e4','N=1e5','N=1e6','status');
for a = 1:numel(kk_list)
    k = kk_list(a);
    fprintf(' %-4d', k);
    if k == 1
        fprintf(' %14s %14s %14s %14s %14s | %s\n','--','--','--','--','--', ...
                'DIVERGENT - no limit exists, see panel I');
    else
        fprintf(' %14.2f', 100*(vals(a,:)-lim_an(a))/abs(lim_an(a)));
        fprintf(' | converges\n');
    end
end

fprintf('\n (II.3)  k = 1 reported as it must be: a running value, not a number.\n');
fprintf(' %-10s %16s %16s\n','N','-Y(i,i+-1)','incr. vs prev.');
prev = NaN;
for b = 1:numel(NN_list)
    if isnan(prev)
        fprintf(' %-10d %16.6e %16s\n', NN_list(b), vals(1,b), '--');
    else
        fprintf(' %-10d %16.6e %16.4e\n', NN_list(b), vals(1,b), vals(1,b)-prev);
    end
    prev = vals(1,b);
end
fprintf(' The value moves by %.4e between N = 1e5 and N = 1e6, i.e. by\n', vals(1,5)-vals(1,4));
fprintf(' -c ln10 = %.4e per decade.  It has no N -> infinity limit.\n', -c*log(10));

fprintf('\n (II.4)  Cross-check of the analytic limit against a direct sum at N = 1e7\n');
fprintf(' %-4s %18s %18s %14s\n','k','limit (analytic)','sum at N=1e7','rel. diff.');
for a = 1:numel(kk_list)
    if kk_list(a) >= 2
        fprintf(' %-4d %18.8e %18.8e %14.2e\n', kk_list(a), lim_an(a), lim_num(a), ...
                (lim_an(a)-lim_num(a))/abs(lim_num(a)));
    end
end

fprintf('\n (II.5)  CONTROL CHECK on panel II (deviation at N = 540, EXACT kernel)\n');
fprintf(' %-4s %16s %16s %12s  %s\n','k','reference [%]','reproduced [%]','deviation','verdict');
ref540 = [NaN -43.98 39.78 33.33 -4.41];
ok_p2 = true;
for a = 1:numel(kk_list)
    if kk_list(a) >= 2
        got = 100*(vals(a,1)-lim_an(a))/abs(lim_an(a));
        ok  = abs(got - ref540(a)) <= 0.5;
        ok_p2 = ok_p2 && ok;
        fprintf(' %-4d %16.2f %16.3f %12.3f  %s\n', kk_list(a), ref540(a), got, got-ref540(a), verdict(ok));
    end
end
fprintf(' panel II controls: %s\n', verdict(ok_p2));

fprintf('-------------------------------------------------------------------------\n');
fprintf(' READING.  At the truncation the validated chain actually uses (N = 540,\n');
fprintf(' i.e. N = Ms/2), the off-diagonal entries are NOT converged: the second\n');
fprintf(' neighbour is %.1f %% away from its limit, the third %+.1f %%, the fifth\n', ...
        100*(vals(2,1)-lim_an(2))/abs(lim_an(2)), 100*(vals(3,1)-lim_an(3))/abs(lim_an(3)));
fprintf(' %+.1f %%, the tenth %+.1f %%, and the first neighbour has no limit at all.\n', ...
        100*(vals(4,1)-lim_an(4))/abs(lim_an(4)), 100*(vals(5,1)-lim_an(5))/abs(lim_an(5)));
fprintf(' The deviations do not decay monotonically in N either: they oscillate\n');
fprintf(' (see k = 5 and k = 10 at N = 1e4), because the residual tail is the\n');
fprintf(' oscillatory remainder of SUM cos(n k d)/n, not a smooth O(1/N) term.\n');
fprintf('=========================================================================\n\n');

fprintf('END OF RUN  %s\n', RUN_DATE);
diary off;

% =========================================================================
% LOCAL FUNCTIONS
% =========================================================================
function k = kappa_vec(n, U_a, U_m, mu_r, GUARD_UA, GUARD_UM)
% kappa_n with the MANDATORY numerical guard.
    n = double(n(:));
    k = ones(size(n));
    m = (n*U_a) <= GUARD_UA;
    nn = n(m);
    if ~isempty(nn)
        cm = coth(min(nn*U_m, GUARD_UM));
        al = 1.0 ./ ( cosh(nn*U_a) + mu_r.*sinh(nn*U_a).*cm );
        k(m) = ( cosh(nn*U_a) - al ) ./ sinh(nn*U_a);
    end
end

function s = offdiag_sum(N, di, dj, dtheta, exact, mu0, L, U_a, U_m, mu_r, ...
                         GUARD_UA, GUARD_UM, n_guard)
% Direct vectorised chunked summation of
%   -Y(i,j) = (4 mu0 L/pi) SUM_{n=1..N} (kappa_n/n) sin(n di/2) sin(n dj/2) cos(n dtheta)
% Beyond n_guard the kernel factor is EXACTLY 1 by the mandated guard, so the
% exact and asymptotic tails are identical there; the split is not an
% approximation.
    N = round(N);
    pref = 4*mu0*L/pi;
    tot  = 0.0;
    CH   = 2e6;
    n0   = 1;
    while n0 <= N
        n1 = min(N, n0 + CH - 1);
        n  = (n0:n1).';
        t  = sin(n*di/2).*sin(n*dj/2).*cos(n*dtheta)./n;
        if exact && n0 <= n_guard
            w = kappa_vec(n, U_a, U_m, mu_r, GUARD_UA, GUARD_UM);
            tot = tot + sum(w.*t);
        else
            tot = tot + sum(t);
        end
        n0 = n1 + 1;
    end
    s = pref*tot;
end

function v = eq29(di, dj, dtheta, c)
% Closed form (29) in its general four-argument form.
    a1 = (di-dj)/2 - dtheta;  a2 = (di-dj)/2 + dtheta;
    a3 = (di+dj)/2 - dtheta;  a4 = (di+dj)/2 + dtheta;
    f  = @(a) log(abs(2*sin(a/2)));
    v  = c*( f(a3) + f(a4) - f(a1) - f(a2) );
end

function v = eq30b(N, d, c, gamma)
% Closed form (30b), first neighbour of a uniform contiguous tiling.
    v = -c*(log(N) + gamma) + c*( log(abs(2*sin(d))) - 2*log(abs(2*sin(d/2))) );
end

function ok = report_rel(name, ref, got, tol, fmt)
    dev = (got - ref)/abs(ref);
    ok  = abs(dev) <= tol;
    fprintf(['%-58s ' fmt ' ' fmt ' %12.3e  %s\n'], name, ref, got, dev, verdict(ok));
end

function s = verdict(ok)
    if ok; s = 'PASS'; else; s = '*** FAIL ***'; end
end
