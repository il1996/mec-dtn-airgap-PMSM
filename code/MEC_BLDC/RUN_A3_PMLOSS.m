%% RUN_A3_PMLOSS - pertes aimants : chaine GELEE et regeneree
%
%  A3. Les 0,342 W publies ne sont produits par AUCUNE chaine constructible
%  (six combinaisons testees, cf. RAPPORT_SPEC_CLAUDE_CODE §T18). On ne
%  cherche pas a les retrouver : on GELE une chaine et on regenere.
%
%  CHAINE RETENUE : pm_loss avec Np >= 181. C'est la seule dont la
%  convergence soit etablie (0,2034 a Np=181 -> 0,2038 a Np=481, soit 0,2 %
%  sur un facteur 2,7 en echantillonnage).
%
%  A VIDE : calcule ici, deux sources de champ (reseau a une dent, maillage).
%  EN CHARGE : pm_loss_load exige la chaine d'entrainement complete
%  (instants, courants triphases, position, carte de potentiel d'induit par
%  ampere) -- c'est la section 5b de BLDC_MEC_COMPLET, non un appel isole.
%  La valeur en charge est donc REPRISE de cette source maitresse et
%  DECLAREE comme telle ; seule la normalisation de T11 est calculee ici.
clear; clc; t0=tic;
diary('A3_pmloss_out.txt'); diary on;
M=machine_bldc();
if isfile('kfringe_ident.mat'), S0=load('kfringe_ident.mat'); kfr=S0.kbest; else, kfr=0.325; end
Np=241; Ms=540; nsh=2; n_rpm=M.FEA.n_nl;

fprintf('=== A3 : pertes aimants, chaine gelee ===\n');
fprintf('  chaine : pm_loss | Np = %d | kfringe = %.4f | %.0f tr/min\n', ...
    Np,kfr,n_rpm);
fprintf('  sigma_pm = %.0f S/m (rho = %.3e ohm.m)\n\n',M.sigma_pm,M.rho_pm);

%% ---- A VIDE : deux sources de champ, MEME formulation ----------------
fprintf('  ---- A VIDE ----\n');
Pl=pm_loss(M,1260,Np,kfr,n_rpm);
phis=linspace(0,2*pi/M.Ns,Np); U0=[]; Usurf=[]; AGm=[];
for q=1:Np
    ME=mesh_bldc(M,Ms,4,3,[],phis(q),kfr,'nl',nsh);
    S=solve_bldc_mesh(ME,1e-9,30,U0); U0=S.U;
    if isempty(Usurf), Usurf=zeros(numel(S.Usurf),Np); AGm=ME.AG; end
    Usurf(:,q)=S.Usurf; %#ok<SAGROW>
end
Pm=pm_loss(M,[],[],kfr,n_rpm,struct('AG',AGm,'Usurf',Usurf,'phis',phis));
PF_nl=0.334;
fprintf('  %-26s %10s %10s %10s\n','source du champ','P (W)','FEA','ecart');
fprintf('  %-26s %10.4f %10.4f %9.1f %%\n','reseau a une dent',Pl,PF_nl,100*(Pl-PF_nl)/PF_nl);
fprintf('  %-26s %10.4f %10.4f %9.1f %%\n','MAILLAGE polaire',Pm,PF_nl,100*(Pm-PF_nl)/PF_nl);
fprintf('  (valeur publiee 0.342 W : NON REGENERABLE, cf. T18)\n');

%% ---- EN CHARGE : normalisation du point de fonctionnement (T11) ------
%  Source : BLDC_MEC_COMPLET section 5b. Les deux colonnes de la Table 10
%  ne sont PAS au meme point : 1257 vs 1340 tr/min, 1.419 vs 1.552 A rms.
%  Les pertes par courants de Foucault varient en omega^2 * B_induit^2,
%  et B_induit est proportionnel au courant.
nM=1257; nF=1340; iM=1.419; iF=1.552;
PmL=1.593; PfL=3.183;                       % BLDC_MEC_COMPLET, Table 10
fw=(nM/nF)^2; fi=(iM/iF)^2; fn=fw*fi;
fprintf('\n  ---- EN CHARGE ----\n');
fprintf('  (valeurs de BLDC_MEC_COMPLET section 5b, source declaree)\n');
fprintf('  %-34s %10.4f W\n','MEC   (1257 tr/min, 1.419 A)',PmL);
fprintf('  %-34s %10.4f W\n','FEA   (1340 tr/min, 1.552 A)',PfL);
fprintf('  ecart BRUT                                  %+9.1f %%\n',100*(PmL-PfL)/PfL);
fprintf('\n  normalisation du point de fonctionnement :\n');
fprintf('    facteur vitesse  (1257/1340)^2 = %.4f\n',fw);
fprintf('    facteur courant  (1.419/1.552)^2 = %.4f\n',fi);
fprintf('    produit                          = %.4f\n',fn);
fprintf('  ecart NORMALISE                             %+9.1f %%\n', ...
    100*((PmL/PfL)/fn-1));

%% ---- consequence redactionnelle --------------------------------------
fprintf(['\n  ---- CONSEQUENCE POUR §5.7 ET §8.4 ----\n' ...
  '  L''argument "a vide valide la chaine, en charge est ouvert" NE TIENT\n' ...
  '  PLUS. A vide la chaine gelee donne %+.1f %% (reseau) et %+.1f %%\n' ...
  '  (maillage), non les +2,4 %% publies. En charge, apres normalisation,\n' ...
  '  le residu tombe a %+.1f %%. Les deux points de fonctionnement sont\n' ...
  '  donc du MEME ordre, et la dichotomie disparait.\n'], ...
  100*(Pl-PF_nl)/PF_nl,100*(Pm-PF_nl)/PF_nl,100*((PmL/PfL)/fn-1));
fprintf('\n  duree %.0f s\n=== A3 termine ===\n',toc(t0));
diary off;
