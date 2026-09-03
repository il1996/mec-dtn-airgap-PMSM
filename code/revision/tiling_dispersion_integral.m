%% TILING_DISPERSION_INTEGRAL  -  task B1 of the revision specification
%
%  PURPOSE. Table 7(b) of the manuscript publishes a tiling dispersion for two
%  quantities only: the fundamental bore flux density B_g1 and the first slot
%  sideband nu = 8. It publishes none for peak flux linkage, L_a, M or L_d.
%  Sections 7.2 and 7.3 therefore contradict each other on whether +0.24 % on
%  flux linkage and +0.10 % on L_d are above or below the resolution of the
%  chain, and the contradiction is UNDECIDABLE with published data. This script
%  measures the four missing dispersions.
%
%  METHOD. The experiment of Table 7(b) is re-run UNCHANGED - same six tilings,
%  same lock N = floor(Ms/2), same two bases, same iron model - and four output
%  quantities are added. B_g1 and nu = 8 are recomputed by the identical call
%  and serve as a bit-exact non-regression gate: if they do not come back, the
%  chain is not the Table 7(b) chain and nothing here is publishable.
%
%  PROVENANCE OF THE CHAIN. The Table 7(b) chain is code/MEC_BLDC/RUN_V1_PMSM_BASE.m,
%  archived at outputs/MEC_BLDC/V1_pmsm_basis_out.txt. Its call is
%      cogging_mec(M, Ms, 0, 3, M.muI, kfr, 1e-6, basis)
%  i.e. LINEAR iron at relative permeability M.muI, no-load fringing constant
%  kfr, a degenerate rotor span of 1e-6 rad and three positions. cogging_mec
%  forces numax = floor(Nsurf/2): that IS the lock.
%
%  THE ONE THING THAT COULD NOT BE TAKEN OVER UNCHANGED, AND WHY. The published
%  call sweeps a span of 1e-6 rad. Over 1e-6 rad the rotor does not move, so a
%  flux linkage read from it has no peak. Peak flux linkage is therefore taken
%  from a SECOND call at the same tiling, same lock, same basis, same iron model
%  and same fringing constant, differing only in the rotor span and the number
%  of positions, for which the published convention of RUN_A1_TABLE7.m is used
%  (span = 2*pi/p, Np = 61). The first call is left untouched, so the
%  non-regression gate is exact by construction; a guard additionally checks
%  that the second call reproduces B_g1 at its own first position, which proves
%  the span change did not move the operator.
%
%  DECLARED, NOT TUNED. The armature problem is run at kfringe = 0 while the
%  no-load problem is run at kfr. This is not a fitted difference: the published
%  header of inductance_mec.m documents that the slot-opening fringe plays two
%  physically distinct roles (radial entry at no load, tangential tooth-tip
%  bridge under armature current) and that reusing the no-load value in the
%  armature problem double-counts tooth-tip leakage. Both values are printed.
%
%  Nothing under code/MEC_BLDC or outputs/MEC_BLDC is written by this script.
%
%  Output: outputs/revision/B1_tiling_dispersion_out.txt

clear; clc; t0=tic;

here = fileparts(mfilename('fullpath'));
root = fileparts(fileparts(here));
outdir = fullfile(root,'outputs','revision');
if ~exist(outdir,'dir'), mkdir(outdir); end
addpath(here);
cd(fullfile(root,'code','MEC_BLDC'));      % the chain sees its usual directory

%  MATLAB's diary APPENDS to an existing file. A transcript must correspond to
%  exactly one dated run, so any previous one is rotated aside rather than grown.
tfile = fullfile(outdir,'B1_tiling_dispersion_out.txt');
if isfile(tfile)
    movefile(tfile, fullfile(outdir, ...
        sprintf('B1_tiling_dispersion_out.superseded_%s.txt', datestr(now,'yyyymmdd_HHMMSS'))));
end
diary(tfile); diary on;

ISO_DATE = datestr(now,'yyyy-mm-dd');

M=machine_bldc(); p=M.p; Ns=M.Ns; Ntc=M.Ntc;
if isfile('kfringe_ident.mat'), S0=load('kfringe_ident.mat'); kfr=S0.kbest; else, kfr=0.325; end

% --- exactly the reduction used by RUN_V1_PMSM_BASE (do not alter) ---
thu=linspace(0,2*pi,2001); thu(end)=[];
amp=@(y,th,k)abs((2/numel(y))*sum(y(:).'.*exp(-1i*k*th(:).')));

% --- exactly the winding of RUN_A1_TABLE7 (do not alter) ---
PA=[1 -2 -15 3 14]; PB=[6 -7 -5 8 4]; PC=[11 -12 -10 13 9];
SPAN_LAM = 2*pi/p; NP_LAM = 61;

Msl = [540 1080 1260 2160 4320 8640];
BAS = {'p0','p1'}; BASLAB = {'p.c.','hat'};

% --- Table 7(b) as published, for the non-regression gate. These are the only
%     transcribed numbers in this script and they are used ONLY as targets to
%     be reproduced, never as results.
PUB.Bg1  = [1.07390 1.07865 1.07755 1.07793 1.07759 1.07743;   % p.c.
            1.07267 1.07801 1.07698 1.07759 1.07741 1.07734];  % hat
PUB.nu8  = [0.021761 0.016944 0.018062 0.017672 0.018014 0.018180;
            0.023000 0.017596 0.018637 0.018020 0.018195 0.018271];
PUB.disp = [0.440 26.126; 0.496 28.514];   % [Bg1 nu8] per basis

fprintf('=== B1 : tiling dispersion of the integral quantities ===\n');
fprintf('  date (ISO)        : %s\n', ISO_DATE);
fprintf('  MATLAB            : %s\n', version);
fprintf('  machine           : PMSM 15 slots / 14 poles, 750 W (machine_bldc)\n');
fprintf('  chain             : cogging_mec + airgap_magnet  (the Table 7(b) chain,\n');
fprintf('                      RUN_V1_PMSM_BASE.m -> V1_pmsm_basis_out.txt)\n');
fprintf('  IRON MODEL        : LINEAR, relative permeability mu_i = %g\n', M.muI);
fprintf('                      no Newton iteration, no saturation curve\n');
fprintf('  fringing constant : no-load problem   kfringe = %.4f\n', kfr);
fprintf('                      armature problem  kfringe = 0  (declared, not fitted:\n');
fprintf('                      see the header of code/MEC_BLDC/inductance_mec.m)\n');
fprintf('  lock              : N = floor(Ms/2), imposed by cogging_mec\n');
fprintf('  tilings           : %s\n', mat2str(Msl));
fprintf('  bases             : p0 (piecewise constant) and p1 (hat)\n');
fprintf('  B_g1 / nu=8 run   : span = 1e-6 rad, Np = 3   (published call, untouched)\n');
fprintf('  flux-linkage run  : span = 2*pi/p, Np = %d   (RUN_A1_TABLE7 convention)\n', NP_LAM);
fprintf('  inductance run    : inductance_mec_basis, numax = floor(Ms/2), I = 1 A\n');
fprintf('  reduction         : thu = 2000 points, amp = 2|mean(B e^{-i k th})|\n\n');

nM=numel(Msl);
V   = nan(nM,6,2);       % [Bg1 nu8 lambda La M Ld] x basis
GATE = true;

for b=1:2
    bs=BAS{b};
    fprintf('  ---------- basis %s (%s) ----------\n', bs, BASLAB{b});
    for k=1:nM
        Mk=Msl(k); Nk=floor(Mk/2); tk=tic;

        % ---- (1) the published call, verbatim -> B_g1 and nu = 8 ----
        R=cogging_mec(M,Mk,0,3,M.muI,kfr,1e-6,bs);
        [Br,~]=R.AG.field(R.U(1:Mk,1),0,thu); Br=Br(:).';
        Bg1=amp(Br,thu,p); nu8=amp(Br,thu,8);

        % ---- (2) same chain, real rotor span -> peak flux linkage ----
        RL=cogging_mec(M,Mk,0,NP_LAM,M.muI,kfr,SPAN_LAM,bs);
        lamf=@(P) Ntc*sum(sign(P(:)).*RL.PhiT(abs(P(:)),:),1);
        lam=[lamf(PA); lamf(PB); lamf(PC)];
        lam_peak=max(abs(lam(1,:)));

        % guard: the span change must not move the operator at position 1
        [Br2,~]=RL.AG.field(RL.U(1:Mk,1),0,thu); Br2=Br2(:).';
        g_span=abs(amp(Br2,thu,p)-Bg1)/Bg1;

        % ---- (3) same tiling and lock -> L_a, M, L_d ----
        RI=inductance_mec_basis(M,Mk,M.muI,0,bs);
        La=RI.La*1e3; Mu=RI.M*1e3; Ld=RI.Ld*1e3;

        % guard: with p0 the derived copy must equal the published function
        g_copy=NaN;
        if strcmp(bs,'p0')
            RI0=inductance_mec(M,Mk,M.muI,0);
            g_copy=max([abs(RI.La-RI0.La) abs(RI.M-RI0.M) abs(RI.Ld-RI0.Ld)]);
        end

        V(k,:,b)=[Bg1 nu8 lam_peak La Mu Ld];

        % ---- non-regression against Table 7(b), value by value ----
        dB=abs(Bg1-PUB.Bg1(b,k)); dn=abs(nu8-PUB.nu8(b,k));
        okB = dB <= 0.5e-5;      % published to 5 decimals
        okn = dn <= 0.5e-6;      % published to 6 decimals
        GATE = GATE && okB && okn;

        fprintf(['  Ms=%5d N=%5d | B_g1 %.5f (pub %.5f, d=%.1e %s) | nu8 %.6f ' ...
                 '(pub %.6f, d=%.1e %s)\n'], Mk,Nk,Bg1,PUB.Bg1(b,k),dB, tern(okB,'OK','FAIL'), ...
                 nu8,PUB.nu8(b,k),dn, tern(okn,'OK','FAIL'));
        fprintf('                    | lambda %.6f Wb | La %.5f  M %.5f  Ld %.5f mH\n', ...
                 lam_peak,La,Mu,Ld);
        fprintf('                    | guard span-invariance %.2e', g_span);
        if strcmp(bs,'p0'), fprintf(' | guard p0-copy==published %.2e', g_copy); end
        fprintf(' | %.0f s\n', toc(tk));

        clear R RL RI RI0
    end
    fprintf('\n');
end

%% ---------------- dispersion and result blocks -------------------------
%  (max-min)/|mean|. The absolute value matters: the mutual inductance M is a
%  NEGATIVE quantity, and (max-min)/mean would report its dispersion with a
%  spurious minus sign. Taking |mean| makes the measure a magnitude for every
%  quantity while leaving it identical to Table 7(b)'s for the positive ones.
dsp=@(x)100*(max(x)-min(x))/abs(mean(x));
D=nan(2,6);
for b=1:2, for j=1:6, D(b,j)=dsp(V(:,j,b)); end, end

%  Supplementary, in the manuscript's OWN convention: Section 6.2 also reports
%  the dispersion with the coarsest tiling dropped (0.096 / 0.113 % on B_g1,
%  5.7 / 7.0 % on nu = 8). The same statistic is formed here for the four new
%  quantities. It is reported ALONGSIDE the six-tiling measure, never instead
%  of it: the six-tiling figure is the one Table 7(b) publishes.
D5=nan(2,6);
for b=1:2, for j=1:6, D5(b,j)=dsp(V(2:end,j,b)); end, end

for b=1:2
    fprintf('=== B1 : tiling dispersion of the integral quantities, locked chain N = floor(Ms/2) ===\n');
    fprintf('iron model : linear, mu_i = %g, kfringe %.4f no-load / 0 armature          basis : %s          date : %s\n', ...
            M.muI, kfr, BASLAB{b}, ISO_DATE);
    fprintf(' Ms     N     B_g1(T)   nu8(T)    lambda(Wb)  L_a(mH)   M(mH)    L_d(mH)\n');
    for k=1:nM
        fprintf('%4d  %5d   %8.5f  %8.6f  %10.6f  %8.5f  %8.5f  %8.5f\n', ...
            Msl(k), floor(Msl(k)/2), V(k,1,b),V(k,2,b),V(k,3,b),V(k,4,b),V(k,5,b),V(k,6,b));
    end
    fprintf(['dispersion (max-min)/mean :  B_g1 %.3f%%  nu8 %.3f%%  lambda %.3f%%  ' ...
             'L_a %.3f%%  M %.3f%%  L_d %.3f%%\n'], D(b,1),D(b,2),D(b,3),D(b,4),D(b,5),D(b,6));
    fprintf(['  same, coarsest tiling dropped (Section 6.2 convention, five tilings) :\n' ...
             '                             B_g1 %.3f%%  nu8 %.3f%%  lambda %.3f%%  ' ...
             'L_a %.3f%%  M %.3f%%  L_d %.3f%%\n\n'], D5(b,1),D5(b,2),D5(b,3),D5(b,4),D5(b,5),D5(b,6));
end

%% ------- the three contested deviations against the floor measured here -----
%  Purely a comparison of numbers already printed above against the three
%  deviations the abstract and Section 8 put in front. No value is adjusted.
fprintf('---- the contested deviations placed against the tiling floor ----\n');
fprintf('  quantity          published deviation   tiling floor (6 tilings, p.c. / hat)   verdict\n');
CONT = { 'peak flux linkage', 0.24, 3; 'L_d',               0.10, 6; 'B_g1',              0.31, 1 };
for r=1:size(CONT,1)
    j=CONT{r,3}; dev=CONT{r,2};
    below = dev < min(D(1,j),D(2,j));
    fprintf('  %-17s   %+.2f %%              %6.3f %% / %6.3f %%                        %s\n', ...
        CONT{r,1}, dev, D(1,j), D(2,j), tern(below,'BELOW the floor','above the floor'));
end
fprintf(['  Reading. A deviation below the movement the tiling alone produces is a\n' ...
         '  BOUND, not a measurement. This block states which of the three are which;\n' ...
         '  it does not decide how the manuscript should word them.\n\n']);

%% ---------------- the non-regression gate ------------------------------
fprintf('---- NON-REGRESSION GATE against Table 7(b) ----\n');
fprintf('  quantity        recomputed   published   difference\n');
for b=1:2
    fprintf('  dispersion B_g1 %s   %8.3f %%  %8.3f %%   %+.4f point\n', ...
        BASLAB{b}, D(b,1), PUB.disp(b,1), D(b,1)-PUB.disp(b,1));
    fprintf('  dispersion nu8  %s   %8.3f %%  %8.3f %%   %+.4f point\n', ...
        BASLAB{b}, D(b,2), PUB.disp(b,2), D(b,2)-PUB.disp(b,2));
end
okd = max(max(abs(D(:,1:2)-PUB.disp))) <= 0.002;
GATE = GATE && okd;
fprintf('\n  columns B_g1 and nu=8 reproduced at all six tilings, both bases : %s\n', tern(GATE,'YES','NO'));
fprintf('  four published dispersions within 0.002 point                   : %s\n', tern(okd,'YES','NO'));
if GATE
    fprintf('  GATE PASSED. The chain is the Table 7(b) chain; the four new\n');
    fprintf('  dispersions above are produced on that same chain.\n');
else
    fprintf('  *** GATE FAILED ***  The chain is NOT the one that produced Table 7(b).\n');
    fprintf('  Nothing downstream of this transcript is publishable. No value has been\n');
    fprintf('  adjusted to close the gap, and none should be.\n');
end

save(fullfile(outdir,'B1_tiling_dispersion.mat'),'V','D','Msl','BAS','PUB','ISO_DATE');
fprintf('\n  duration %.0f s\n=== B1 complete ===\n',toc(t0));
diary off;

function s=tern(c,a,b), if c, s=a; else, s=b; end, end
