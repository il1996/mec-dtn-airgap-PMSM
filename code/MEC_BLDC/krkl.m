function [mmf_m, mmf_g, kr, kl, disp_m] = krkl(AG, Us, phq, M, L, mu0)
%KRKL  FMM par aimant et par entrefer, facteurs k_r et k_l.
%
%   [mmf_m,mmf_g,kr,kl,disp_m] = krkl(AG,Us,phq,M,L,mu0) evalue, pour un
%   potentiel de surface Us et une position rotor phq :
%     mmf_m  : FMM moyenne par aimant (A)
%     mmf_g  : FMM moyenne par entrefer (A)
%     kr     : facteur de reluctance, mmf_m / mmf_g
%     kl     : facteur de fuite
%     disp_m : dispersion des quatorze FMM d'aimant (%)
%
%   BLOC C1. Cette fonction etait une fonction LOCALE de RUN_A1_TABLE7.m
%   (l. 139), donc non appelable : le bloc A5 a du la dupliquer pour
%   regenerer la colonne Lumped. Promue ici pour qu'une seule definition
%   existe. Les scripts qui la portent encore en local la masquent
%   localement -- leur comportement est inchange -- mais tout nouveau
%   script doit appeler CELLE-CI.
%
%   Voir aussi : cogging_mec, airgap_magnet, RUN_A1_TABLE7, RUN_A5_LUMPED.

nu = AG.nu; Nm = M.Nm;
Usc = AG.Wc*Us; Uss = AG.Ws*Us;
pic = AG.alphau.*Usc + AG.betasrc.*cos(nu*phq);
pis = AG.alphau.*Uss + AG.betasrc.*sin(nu*phq);
thq = linspace(0,2*pi,2001); thq(end) = []; dth_ = thq(2)-thq(1);
phi_i = (pic.'*cos(nu*thq) + pis.'*sin(nu*thq));
phi_b = (Usc.'*cos(nu*thq) + Uss.'*sin(nu*thq));
mm = zeros(1,Nm); mg = zeros(1,Nm);
for k = 1:Nm
    c = (k-1)*2*pi/Nm + phq;
    [~,j] = min(abs(angle(exp(1i*(thq-c)))));
    mm(k) = abs(phi_i(j));
    mg(k) = abs(phi_i(j) - phi_b(j));
end
mmf_m = mean(mm); mmf_g = mean(mg); kr = mmf_m/mmf_g;
disp_m = 100*std(mm)/mmf_m;

%  facteur de fuite : flux utile au bore rapporte au flux d'aimant
Brc_b = AG.g.*Usc + AG.brm.*cos(nu*phq);
Brs_b = AG.g.*Uss + AG.brm.*sin(nu*phq);
Br_bore = (Brc_b.'*cos(nu*thq) + Brs_b.'*sin(nu*thq));
Phi_bore = sum(abs(Br_bore))*dth_*M.Rsi*L/2;
Mn = AG.Mn; Pn = AG.Pn; Um = AG.Um; rmi = M.rmi; mur = M.mu_r;
Mc = Mn.*cos(nu*phq); Msn = Mn.*sin(nu*phq);
Pcv = Pn.*cos(nu*phq); Psv = Pn.*sin(nu*phq);
Sm_ = sinh(nu*Um); Cm_ = cosh(nu*Um);
ac  = (pic - Pcv*M.Rro + Pcv*rmi.*Cm_)./Sm_;
asv = (pis - Psv*M.Rro + Psv*rmi.*Cm_)./Sm_;
r = 0.5*(M.Rro + rmi); v = log(r/rmi);
cmid = -mu0*mur*(Pcv + (ac .*nu.*cosh(nu*v) - Pcv*rmi.*nu.*sinh(nu*v))/r) + mu0*Mc;
smid = -mu0*mur*(Psv + (asv.*nu.*cosh(nu*v) - Psv*rmi.*nu.*sinh(nu*v))/r) + mu0*Msn;
Br_in = (cmid.'*cos(nu*thq) + smid.'*sin(nu*thq));
Phi_in = sum(abs(Br_in))*dth_*r*L/2;
kl = Phi_bore/Phi_in;
end
