%% RUN_DIAG_COGCONV - la detente du maillage est-elle convergee ?
%  RUN_COG_MESH a rapporte 25.0 iterations de Newton PAR POSITION, soit
%  exactement le plafond itmax=25 impose, et cela aux 841 positions. Newton
%  n'a donc atteint la tolerance 1e-10 NULLE PART : le residu final est
%  inconnu. Or la detente est une petite difference de grandes quantites --
%  0.6 mN.m sur un couple electromagnetique de 4.9 N.m, soit 1e-4 en
%  relatif. Une solution non convergee la contamine directement.
%  On mesure donc (a) le plancher de residu reellement atteignable, et
%  (b) la sensibilite du couple a la tolerance.
clear; clc;
M=machine_bldc(); Nm=M.Nm; L=M.ls; mu0=4*pi*1e-7;
if isfile('kfringe_ident.mat'), S0=load('kfringe_ident.mat'); kfr=S0.kbest; else, kfr=0.325; end
Ms=900; nsh=2; phis=linspace(0,2*pi/Nm,841);
qs=[1 100 200 300];
fprintf('=== Convergence de Newton et sensibilite du couple ===\n');
fprintf('  %6s %10s %6s %12s %14s\n','pos','tol','iter','residu','T (mN.m)');
Tt=zeros(numel(qs),4); tols=[1e-6 1e-8 1e-10 1e-12];
for k=1:numel(qs)
    q=qs(k);
    for j=1:4
        ME=mesh_bldc(M,Ms,4,3,[],phis(q),kfr,'nl',nsh);
        S=solve_bldc_mesh(ME,tols(j),80);
        AG=ME.AG; nu=AG.nu;
        Usc=AG.Wc*S.Usurf; Uss=AG.Ws*S.Usurf;
        cq=cos(nu*phis(q)); sq=sin(nu*phis(q));
        Brc=AG.bru.*Usc+AG.brmq.*cq;  Brs=AG.bru.*Uss+AG.brmq.*sq;
        Btc=-AG.btu.*Uss-AG.btmq.*sq; Bts=AG.btu.*Usc+AG.btmq.*cq;
        Tt(k,j)=(L*M.rmid^2/mu0)*pi*sum(Brc.*Btc+Brs.*Bts)*1e3;
        fprintf('  %6d %10.0e %6d %12.2e %14.4f\n',q,tols(j),S.iter,S.err,Tt(k,j));
    end
end
fprintf('\n  --- sensibilite du couple a la tolerance ---\n');
fprintf('  variation max de T entre tol=1e-6 et tol=1e-12 : %.4f mN.m\n', ...
    max(max(Tt,[],2)-min(Tt,[],2)));
fprintf('  a comparer a la detente annoncee : 6.958 mN.m c-c\n');
fprintf('  et a la reference FEA            : 0.604 +/- 1.401 mN.m\n');
