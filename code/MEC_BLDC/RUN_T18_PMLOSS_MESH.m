%% RUN_T18_PMLOSS_MESH - pertes aimant pilotees par le MAILLAGE
%
%  T18, cellule restante du Tableau 7. La formulation des pertes n'est PAS
%  reecrite : on reutilise pm_loss et on lui passe les potentiels de surface
%  du maillage polaire au lieu de ceux du reseau a une dent. Seule la SOURCE
%  du champ change entre les colonnes 'Lumped' et 'Mesh'.
%
%  Une reecriture avait donne 112 W contre 0,334 attendu -- facteur 336 --
%  en omettant la division par nu (le potentiel vecteur est l'INTEGRALE de
%  Br sur theta) et en balayant une periode electrique au lieu d'un PAS
%  D'ENCOCHE. D'ou cette regle : piloter la chaine, ne pas la refaire.
clear; clc; t0=tic;
diary('T18_pmloss_out.txt'); diary on;
M=machine_bldc();
if isfile('kfringe_ident.mat'), S0=load('kfringe_ident.mat'); kfr=S0.kbest; else, kfr=0.325; end
Ms=540; nsh=2; Np=61; n_rpm=M.FEA.n_nl;

fprintf('=== T18 : pertes aimant pilotees par le maillage ===\n');
fprintf('  Ms=%d nsh=%d, %d positions sur UN PAS D''ENCOCHE, %.0f tr/min\n', ...
    Ms,nsh,Np,n_rpm);

%% ---- balayage du maillage sur un pas d'encoche -----------------------
phis=linspace(0,2*pi/M.Ns,Np);
U0=[]; Usurf=[];
for q=1:Np
    ME=mesh_bldc(M,Ms,4,3,[],phis(q),kfr,'nl',nsh);
    S=solve_bldc_mesh(ME,1e-9,30,U0); U0=S.U;
    if isempty(Usurf), Usurf=zeros(numel(S.Usurf),Np); AGm=ME.AG; end
    Usurf(:,q)=S.Usurf; %#ok<SAGROW>
end
fprintf('  balayage maillage termine (%.0f s)\n',toc(t0));

%% ---- pertes : MEME formulation, deux sources -------------------------
ext=struct('AG',AGm,'Usurf',Usurf,'phis',phis);
[Pm,Dm]=pm_loss(M,[],[],kfr,n_rpm,ext);
[Pl,Dl]=pm_loss(M,1260,Np,kfr,n_rpm);

fprintf('\n  %-28s %12s %12s %12s\n','source du champ','P (W)','FEA','ecart');
fprintf('  %-28s %12.4f %12.4f %11.1f %%\n','reseau a une dent',Pl,0.334,100*(Pl-0.334)/0.334);
fprintf('  %-28s %12.4f %12.4f %11.1f %%\n','MAILLAGE polaire',Pm,0.334,100*(Pm-0.334)/0.334);
fprintf('\n  f de denture %.1f Hz | profondeur de peau %.2f mm\n',Dm.fslot,Dm.delta*1e3);
fprintf('  attenuation a travers l''aimant : rang 8 %.3f | rang 22 %.3e\n', ...
    Dm.att8,Dm.att22);
fprintf('\n  duree %.0f s\n=== termine ===\n',toc(t0));
diary off;
