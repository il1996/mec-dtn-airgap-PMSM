%% RUN_IND  Inductances par MEC couple + FMM distribuee dans l'encoche (ameli. 3)
%  Cibles FEA (magnetostatique "Armature-Field", i_a = 1 A) :
%      L_a = 50.21 mH, L_ba = L_ca = -2.13 mH, L_d = L_a - M = 52.34 mH.
%  Le modele FEA est 2D -> pas de fuite de tetes de bobines : la confrontation
%  porte donc sur la partie 2D seule (lecon P4 de MEC_IM : ne jamais melanger
%  correction de lambda d'encoche et identification de fuite d'extremite 3D).
%
%  TROIS CORRECTIONS PHYSIQUES par rapport a la version precedente (-11.3 %) :
%   (1) permeance d'encoche integree sur le TRAPEZE REEL (et non a largeur
%       moyenne) : le poids (y/h)^2 est maximal la ou la largeur est MINIMALE ;
%   (2) ampere-tours comptes PAR ENCOCHE : dans un bobinage dentaire, deux
%       cotes de bobine voisins peuvent appartenir a la MEME phase et leurs
%       ampere-tours S'AJOUTENT (encoches 1,2,14,15 ici) ;
%   (3) part d'entrefer issue du RESEAU couple (FMM injectee comme source de
%       branche dans chaque dent), et non d'un entrefer uniforme + Carter.
clear; clc;
M=machine_bldc();
FEA_La=M.FEA.La*1e3; FEA_M=-2.132; FEA_Ld=M.FEA.Ld*1e3;
fprintf('=== INDUCTANCES : MEC couple + FMM distribuee dans l''encoche ===\n');

R=inductance_mec(M,840,5000,0);          % kfringe = 0 pour l'induit (cf. en-tete)

fprintf('\n[1] Permeance d''encoche - integration sur le TRAPEZE REEL\n');
fprintf('    largeur %.3f mm (fond, culasse) -> %.3f mm (sommet, entrefer)\n',...
        M.ws2*1e3,M.ws1*1e3);
fprintf('    lambda corps  EXACT  = %.4f   (largeur moyenne : %.4f -> %+.1f %% d''erreur)\n',...
        R.lam_body,R.lam_body_approx,(R.lam_body_approx-R.lam_body)/R.lam_body*100);
fprintf('    lambda biseau        = %.4f   | ouverture = %.4f | bec = %.4f\n',...
        R.lam_wedge,M.hs0/M.ws0,R.lam_tot-R.lam_slot);
fprintf('    lambda TOTAL         = %.4f\n',R.lam_tot);

fprintf('\n[2] Ampere-tours par encoche (unites Ntc)\n');
fprintf('    encoche : '); fprintf('%5d',1:M.Ns); fprintf('\n');
for k=1:3
    fprintf('    phase %c : ',64+k); fprintf('%5.0f',R.AT(:,k)/M.Ntc); fprintf('\n');
end
fprintf('    Somme des carres phase A = %.0f Ntc^2 ; la formule usuelle\n',sum(R.AT(:,1).^2)/M.Ntc^2);
fprintf('    (4m/Q)*Nph^2 en supposerait %.0f (+%.0f %%) -> comptage exact requis.\n',...
        (4*M.m/M.Ns)*(M.Ntc*5)^2/M.Ntc^2, ...
        ((4*M.m/M.Ns)*25-sum(R.AT(:,1).^2)/M.Ntc^2)/(sum(R.AT(:,1).^2)/M.Ntc^2)*100);

fprintf('\n[3] Decomposition\n');
fprintf('    %-30s %10s %10s\n','','propre','mutuelle');
fprintf('    %-30s %10.3f %10.3f mH\n','part d''entrefer (reseau couple)',R.L_gap*1e3,R.M_gap*1e3);
fprintf('    %-30s %10.3f %10.3f mH\n','fuite d''encoche (distribuee)',R.L_slot*1e3,R.M_slot*1e3);
fprintf('    %-30s %10.3f %10.3f mH\n','TOTAL',R.La*1e3,R.M*1e3);

fprintf('\n[4] CONFRONTATION FEA  (2 cibles INDEPENDANTES : propre ET mutuelle)\n');
row=@(n,a,f) fprintf('    %-30s %10.3f %10.3f %+8.1f %%\n',n,a,f,(a-f)/f*100);
row('Inductance propre La (mH)',R.La*1e3,FEA_La);
row('Mutuelle M (mH)',R.M*1e3,FEA_M);
row('Inductance synchrone Ld (mH)',R.Ld*1e3,FEA_Ld);
fprintf('    (version precedente, lambda approche + entrefer uniforme : 44.5 mH, -11.3 %%)\n');

fprintf('\n[5] Convergence et sensibilites\n');
fprintf('    %-22s %10s %10s %10s\n','parametre','La(mH)','M(mH)','Ld(mH)');
for Nsf=[420 840 1260 1680]
    r=inductance_mec(M,Nsf,5000,0);
    fprintf('    Nsurf = %-14d %10.3f %10.3f %10.3f\n',Nsf,r.La*1e3,r.M*1e3,r.Ld*1e3);
end
for mu=[1000 5000 20000 1e5]
    r=inductance_mec(M,840,mu,0);
    fprintf('    muI   = %-14g %10.3f %10.3f %10.3f\n',mu,r.La*1e3,r.M*1e3,r.Ld*1e3);
end
fprintf('\n    Effet du DOUBLE COMPTAGE si l''on reutilise kfringe a vide :\n');
for kk=[0 0.75]
    r=inductance_mec(M,840,5000,kk);
    fprintf('    kfringe = %-12.2f %10.3f %10.3f %10.3f  (%+.1f %% sur La)\n',...
            kk,r.La*1e3,r.M*1e3,r.Ld*1e3,(r.La*1e3-FEA_La)/FEA_La*100);
end

figure('Name','Inductances MEC','Color','w','Position',[70 70 950 400]);
subplot(1,2,1);
bar([R.L_gap R.L_slot; 0 0]*1e3,'stacked'); hold on;
yline(FEA_La,'r--','LineWidth',1.6);
grid on; xlim([0.4 1.6]); set(gca,'XTick',1,'XTickLabel',{'L_a MEC'});
ylabel('inductance (mH)'); legend('entrefer (reseau)','fuite d''encoche','FEA','Location','best');
title(sprintf('L_a = %.1f mH (FEA %.1f, %+.1f %%)',R.La*1e3,FEA_La,(R.La*1e3-FEA_La)/FEA_La*100));
subplot(1,2,2);
v=[R.La*1e3 FEA_La; -R.M*1e3 -FEA_M; R.Ld*1e3 FEA_Ld];
bar(v); grid on; set(gca,'XTickLabel',{'L_a','-M','L_d'});
ylabel('mH'); legend('MEC','FEA','Location','best'); title('Confrontation');
saveas(gcf,fullfile(fileparts(mfilename('fullpath')),'FIG_ind.png'));
