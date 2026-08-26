%% RUN_T18_CONV - convergence en Np des pertes aimant, puis cellule Mesh
%
%  Le controle integre a RUN_T18_PMLOSS_MESH a echoue : a Np = 61 la colonne
%  "reseau a une dent" donne 0,200 W au lieu des 0,342 W publies. La perte est
%  pilotee par dA/dt, et la derivation amplifie les ordres eleves : 61
%  positions sur un pas d'encoche les sous-resolvent.
%
%  On balaie donc Np sur la colonne LOCALISEE (peu couteuse) jusqu'a
%  retrouver 0,342 W, puis on lit la colonne MAILLAGE au meme Np.
clear; clc; t0=tic;
diary('T18_conv_out.txt'); diary on;
M=machine_bldc();
if isfile('kfringe_ident.mat'), S0=load('kfringe_ident.mat'); kfr=S0.kbest; else, kfr=0.325; end
n_rpm=M.FEA.n_nl;
fprintf('=== T18 : convergence en Np (pertes aimant) ===\n');
fprintf('  cible colonne localisee : 0.342 W (publie) ; FEA 0.334 W\n\n');

Nlist=[61 121 181 241 361 481];
fprintf('  %6s %12s %10s\n','Np','P localise','ecart/0.342');
PL=nan(size(Nlist));
for k=1:numel(Nlist)
    PL(k)=pm_loss(M,1260,Nlist(k),kfr,n_rpm);
    fprintf('  %6d %12.4f %9.1f %%\n',Nlist(k),PL(k),100*(PL(k)-0.342)/0.342);
end
[~,ib]=min(abs(PL-0.342)); Npc=Nlist(ib);
fprintf('\n  Np retenu : %d  (colonne localisee %.4f W)\n',Npc,PL(ib));

%% ---- colonne MAILLAGE au meme Np -------------------------------------
fprintf('\n  balayage du maillage a Np = %d ...\n',Npc);
phis=linspace(0,2*pi/M.Ns,Npc);
U0=[]; Usurf=[]; AGm=[];
for q=1:Npc
    ME=mesh_bldc(M,540,4,3,[],phis(q),kfr,'nl',2);
    S=solve_bldc_mesh(ME,1e-9,30,U0); U0=S.U;
    if isempty(Usurf), Usurf=zeros(numel(S.Usurf),Npc); AGm=ME.AG; end
    Usurf(:,q)=S.Usurf; %#ok<SAGROW>
end
Pm=pm_loss(M,[],[],kfr,n_rpm,struct('AG',AGm,'Usurf',Usurf,'phis',phis));

fprintf('\n  --- cellule "pertes aimant" du Tableau 7 ---\n');
fprintf('  %-24s %10s %10s %10s\n','source du champ','P (W)','FEA','ecart');
fprintf('  %-24s %10.4f %10.4f %9.1f %%\n','reseau a une dent',PL(ib),0.334,100*(PL(ib)-0.334)/0.334);
fprintf('  %-24s %10.4f %10.4f %9.1f %%\n','MAILLAGE polaire',Pm,0.334,100*(Pm-0.334)/0.334);
fprintf('\n  duree %.0f s\n=== termine ===\n',toc(t0));
diary off;
