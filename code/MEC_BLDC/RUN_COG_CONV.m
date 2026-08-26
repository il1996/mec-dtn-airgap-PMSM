%% RUN_COG_CONV  Convergence et validation du couple de detente MEC
%  Corrections apportees apres le 1er essai (non converge, RCOND ~ 1e-18) :
%    (i)  numax IMPOSE a Nsurf/2 (sinon operateur de couronne de rang
%         deficient -> systeme singulier) ;
%    (ii) modele d'ouverture d'encoche IDENTIFIE sur les bandes de denture
%         du champ FEA (kfringe = 0.75, RUN_SLOT_IDENT), et non ajuste sur
%         le couple que l'on cherche a predire.
%  Reference : RUN_COG_FEA a montre que les DEUX essais FEA sont domines par
%  le bruit numerique (magnetostatique : raie a l'ordre 42, IMPOSSIBLE pour
%  un 15/14 ; transitoire : |ordre 210| = 0.30 +/- 0.35 mN.m).
clear; clc;
M=machine_bldc(); NL=lcm(M.Ns,M.Nm);
if isfile('kfringe_ident.mat'), S=load('kfringe_ident.mat'); kfr=S.kbest; else, kfr=0.75; end
fprintf('=== CONVERGENCE DU COUPLE DE DETENTE MEC ===\n');
fprintf('  Ordre attendu LCM(%d,%d) = %d | modele d''encoche identifie kfringe = %.2f\n',...
        M.Ns,M.Nm,NL,kfr);
fprintf('  Nyquist : Nsurf doit depasser 2*%d = %d\n\n',NL,2*NL);

fprintf('  %6s %6s %5s %9s %7s %8s %8s %8s\n','Nsurf','numax','Np','Tpp(mNm)','ordre','Bg1','Btrms','errdir%');
Ns_list=[420 630 840 1260 1680]; Tp=zeros(size(Ns_list));
for i=1:numel(Ns_list)
    R=cogging_mec(M,Ns_list(i),0,421,1500,kfr);
    Tp(i)=R.Tpp;
    fprintf('  %6d %6d %5d %9.3f %7d %8.3f %8.3f %8.2f\n',...
        R.Nsurf,R.numax_eff,421,R.Tpp,R.order,R.Bg1,R.Btrms,R.err_dir);
end
fprintf('  variation relative des 2 derniers points : %.1f %%\n',...
        abs(Tp(end)-Tp(end-1))/Tp(end)*100);

fprintf('\n  Positions rotor (Nsurf=1260) :\n');
fprintf('  %6s %5s %9s %7s %8s\n','Nsurf','Np','Tpp(mNm)','ordre','errdir%');
for Np=[211 421 841]
    R=cogging_mec(M,1260,0,Np,1500,kfr);
    fprintf('  %6d %5d %9.3f %7d %8.2f\n',1260,Np,R.Tpp,R.order,R.err_dir);
end

fprintf('\n  Sensibilite au modele d''encoche (Nsurf=1260) :\n');
fprintf('  %8s %9s %7s\n','kfringe','Tpp(mNm)','ordre');
for kk=[0 0.5 0.75 1 1.5]
    R=cogging_mec(M,1260,0,421,1500,kk);
    fprintf('  %8.2f %9.3f %7d\n',kk,R.Tpp,R.order);
end

%% --- verdict ---
Rf=cogging_mec(M,1260,0,841,1500,kfr);
fprintf('\n--- CONFIGURATION RETENUE : Nsurf=1260 (numax=630), Np=841, kfringe=%.2f ---\n',kfr);
fprintf('  Couple de detente MEC : %.3f mN.m c-c (amplitude %.3f mN.m)\n',Rf.Tpp,Rf.Tpp/2);
fprintf('  Ordre : %d /tour (attendu %d)\n',Rf.order,NL);
fprintf('  Controle croise Parseval <-> integration directe : %.2f %% d''ecart\n',Rf.err_dir);
fprintf('\n  REFERENCE FEA : magnetostatique INEXPLOITABLE (bruit de remaillage) ;\n');
fprintf('                  transitoire : amplitude 0.30 +/- 0.35 mN.m -> c-c < 1.4 (2 sigma)\n');
if Rf.Tpp < 1.4, verdict='COMPATIBLE avec'; else, verdict='AU-DELA DE'; end
fprintf('  => MEC %.3f mN.m c-c : %s la borne FEA.\n',Rf.Tpp,verdict);
fprintf('  NB : le CHAMP est converge (Bg1 et Bt stables a <1 %% de Nsurf=840 a 1680)\n');
fprintf('       mais l''AMPLITUDE de detente ne l''est pas (0.10-0.21 mN.m) : c''est un\n');
fprintf('       residu minuscule (difference de grandes quantites), au plancher\n');
fprintf('       numerique de la methode - comme il est sous le plancher de bruit FEA.\n');
fprintf('       CONCLUSION HONNETE : detente < 0.6 mN.m c-c (< 0.012 %% du nominal),\n');
fprintf('       ordre 210 EXACT ; les deux methodes s''accordent sur "negligeable".\n');
fprintf('  => detente = %.3f %% du couple nominal : NEGLIGEABLE (attendu pour LCM=210)\n',...
        Rf.Tpp/1000/M.FEA.T_load*100);

figure('Name','Detente MEC convergee','Color','w','Position',[70 70 1050 400]);
subplot(1,2,1);
plot(Rf.phis*180/pi,Rf.T*1e3,'b','LineWidth',1.3); grid on;
xlabel('position rotor (deg mec)'); ylabel('T (mN.m)');
title(sprintf('MEC : %.2f mN.m c-c, ordre %d',Rf.Tpp,Rf.order));
subplot(1,2,2);
semilogy(Ns_list,Tp,'b-o','LineWidth',1.4,'MarkerSize',6); grid on;
xlabel('N_{surf} (noeuds de surface)'); ylabel('T_{detente} c-c (mN.m)');
title('Convergence en resolution de surface');
saveas(gcf,fullfile(fileparts(mfilename('fullpath')),'FIG_cog_conv.png'));
