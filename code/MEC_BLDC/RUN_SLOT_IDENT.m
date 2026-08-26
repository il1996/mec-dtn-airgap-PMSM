%% RUN_SLOT_IDENT  Identification FEM du modele d'ouverture d'encoche (A1)
%  Le COUPLE FEA est domine par le bruit (RUN_COG_FEA) mais le PROFIL DE
%  CHAMP FEA est fiable (1001 pts / 360 deg = 0.36 deg ; creux d'encoche
%  3.53 deg). On identifie le seul parametre libre du modele d'encoche - la
%  frange kfringe - sur la CIBLE PHYSIQUE PERTINENTE : l'amplitude des
%  BANDES LATERALES DE DENTURE du champ d'entrefer, ordres |Ns +/- p| = 8 et
%  22 (et 23, 37 au 2e rang). Une metrique point-a-point serait trompeuse :
%  elle est dominee par le desalignement de phase du fin ripple.
%  Demarche d'identification FEM de Gyselinck [29] (A1 du programme).
clear; clc;
M=machine_bldc(); Ns=M.Ns; p=M.p;
fea=M.FEA.dir;
d4=readmatrix(fullfile(fea,'magnetostique(Magnetic_loading)','Calculator Expressions Plot 4.tab'),'FileType','text','NumHeaderLines',1);
angF=d4(:,2)*pi/180; BrF=d4(:,3);
% re-echantillonnage uniforme sur un tour complet (FFT propre)
thu=linspace(0,2*pi,3601); thu(end)=[];
BrFu=interp1([angF;2*pi],[BrF;BrF(1)],thu,'linear','extrap');

amp=@(B,th,k) 2*abs(mean(B.*exp(-1i*k*th)));
ords=[p, abs(Ns-p), Ns+p, abs(2*Ns-p), 2*Ns+p];    % 7, 8, 22, 23, 37
nam ={'p=7 (fond.)','8 = |Ns-p|','22 = Ns+p','23 = |2Ns-p|','37 = 2Ns+p'};
aF=arrayfun(@(k) amp(BrFu,thu,k), ords);

fprintf('=== IDENTIFICATION DU MODELE D''OUVERTURE D''ENCOCHE (A1) ===\n');
fprintf('  Harmoniques de B_r FEA (mi-entrefer) :\n');
for i=1:numel(ords), fprintf('    ordre %2d  %-14s : %.4f T\n',ords(i),nam{i},aF(i)); end
fprintf('  FEA : Br moy %.4f T, crete %.4f T\n',mean(abs(BrFu)),max(BrFu));

kf=[0 0.25 0.5 0.75 1 1.5 2 3];
fprintf('\n  %6s %8s %8s %8s %8s %8s %9s\n','kfring','a8(T)','a22(T)','a23(T)','Brmoy','Brcrete','err_dent%');
best=inf; kbest=NaN;
for i=1:numel(kf)
    R=cogging_mec(M,840,0,21,M.muI,kf(i));
    th=R.thq(1:end-1); Bm=R.Br(1:end-1);
    aM=arrayfun(@(k) amp(Bm,th,k), ords);
    % erreur sur les bandes de DENTURE seulement (ordres 8, 22, 23)
    ed=norm((aM(2:4)-aF(2:4))./aF(2:4))/sqrt(3)*100;
    fprintf('  %6.2f %8.4f %8.4f %8.4f %8.4f %8.4f %9.2f\n',...
            kf(i),aM(2),aM(3),aM(4),mean(abs(Bm)),max(Bm),ed);
    if ed<best, best=ed; kbest=kf(i); Rb=R; aMb=aM; end
end
fprintf('\n  >> kfringe IDENTIFIE = %.2f  (erreur sur les bandes de denture = %.1f %%)\n',kbest,best);
fprintf('  Harmoniques de denture : MEC %.4f / %.4f / %.4f  vs FEA %.4f / %.4f / %.4f T\n',...
        aMb(2),aMb(3),aMb(4),aF(2),aF(3),aF(4));
fprintf('  Fondamental : MEC %.4f vs FEA %.4f T (%+.1f %%)\n',aMb(1),aF(1),(aMb(1)-aF(1))/aF(1)*100);
fprintf('  Br moyen    : MEC %.4f vs FEA %.4f T (%+.1f %%)\n',...
        mean(abs(Rb.Br)),mean(abs(BrFu)),(mean(abs(Rb.Br))-mean(abs(BrFu)))/mean(abs(BrFu))*100);
save(fullfile(fileparts(mfilename('fullpath')),'kfringe_ident.mat'),'kbest');

figure('Name','Identification denture','Color','w','Position',[60 60 1050 420]);
subplot(1,2,1);
kk=1:40; sF=arrayfun(@(k) amp(BrFu,thu,k),kk);
th=Rb.thq(1:end-1); Bm=Rb.Br(1:end-1); sM=arrayfun(@(k) amp(Bm,th,k),kk);
stem(kk,sF,'r','filled','MarkerSize',4); hold on;
stem(kk+0.3,sM,'b','filled','MarkerSize',4); grid on; set(gca,'YScale','log');
xlabel('ordre spatial'); ylabel('|B_r| (T)'); legend('FEA','MEC identifie');
title(sprintf('Spectre de B_r (k_{fringe}=%.2f)',kbest));
subplot(1,2,2);
plot(thu*180/pi,BrFu,'r--','LineWidth',1.1); hold on;
plot(Rb.thq*180/pi,Rb.Br,'b','LineWidth',1.3); grid on; xlim([0 75]);
xlabel('\theta mec (deg)'); ylabel('B_r (T)'); legend('FEA','MEC','Location','south');
title('Profil (zoom 3 poles)');
saveas(gcf,fullfile(fileparts(mfilename('fullpath')),'FIG_slot_ident.png'));
