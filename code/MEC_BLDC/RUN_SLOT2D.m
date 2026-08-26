%% RUN_SLOT2D  Sous-modele 2D de la bouche d'encoche vs scalaire identifie
%  Compare, sur les BANDES LATERALES DE DENTURE du champ FEA (ordres 8, 22,
%  23, 37 - la cible physique), deux representations de l'ouverture :
%    (a) cogging_mec  : UNE conductance de frange kfringe IDENTIFIEE (0.75) ;
%    (b) cogging_mec2 : la bouche MAILLEE avec sa geometrie reelle
%        (canal ws0 x hs0, biseau ws0->ws1 sur hs1, parois = becs de dents),
%        SANS AUCUN PARAMETRE LIBRE.
clear; clc;
M=machine_bldc(); Ns=M.Ns; p=M.p;
if isfile('kfringe_ident.mat'), S=load('kfringe_ident.mat'); kfr=S.kbest; else, kfr=0.75; end

% --- reference FEA ---
fea=M.FEA.dir;
d4=readmatrix(fullfile(fea,'magnetostique(Magnetic_loading)','Calculator Expressions Plot 4.tab'),'FileType','text','NumHeaderLines',1);
d2=readmatrix(fullfile(fea,'magnetostique(Magnetic_loading)','Calculator Expressions Plot 2.tab'),'FileType','text','NumHeaderLines',1);
angF=d4(:,2)*pi/180; BrF=d4(:,3); BtF=d2(:,2);
thu=linspace(0,2*pi,3601); thu(end)=[];
BrFu=interp1([angF;2*pi],[BrF;BrF(1)],thu,'linear','extrap');
amp=@(B,th,k) 2*abs(mean(B.*exp(-1i*k*th)));
ords=[p 8 22 23 37];
aF=arrayfun(@(k) amp(BrFu,thu,k),ords);
BtF_rms=sqrt(mean(BtF.^2));

fprintf('=== SOUS-MODELE 2D DE LA BOUCHE D''ENCOCHE ===\n');
fprintf('  Geometrie maillee : canal %.1f x %.1f mm, biseau %.1f -> %.2f mm sur %.1f mm\n',...
        M.ws0*1e3,M.hs0*1e3,M.ws0*1e3,M.ws1*1e3,M.hs1*1e3);
fprintf('  FEA : Bg1 %.4f | a8 %.4f | a22 %.4f | a23 %.4f | a37 %.4f | Bt_rms %.4f\n',...
        aF(1),aF(2),aF(3),aF(4),aF(5),BtF_rms);

md=@(R) arrayfun(@(k) amp(R.Br(1:end-1),R.thq(1:end-1),k),ords);
err=@(a) norm((a(2:4)-aF(2:4))./aF(2:4))/sqrt(3)*100;

fprintf('\n  %-26s %8s %8s %8s %8s %8s %9s\n','modele','Bg1','a8','a22','a23','Bt_rms','err_dent%');
R1=cogging_mec(M,840,0,21,1500,kfr);  a1=md(R1);
fprintf('  %-26s %8.4f %8.4f %8.4f %8.4f %8.4f %9.1f\n',...
        sprintf('(a) kfringe=%.2f identifie',kfr),a1(1),a1(2),a1(3),a1(4),R1.Btrms,err(a1));

fprintf('  --- (b) bouche maillee, convergence en nm ---\n');
for nm=[2 4 6 8 12]
    R2=cogging_mec2(M,840,nm,21,1500);  a2=md(R2);
    fprintf('  %-26s %8.4f %8.4f %8.4f %8.4f %8.4f %9.1f\n',...
            sprintf('    nm = %d couches',nm),a2(1),a2(2),a2(3),a2(4),R2.Btrms,err(a2));
end

% --- meilleur nm retenu ---
nmb=8; R2=cogging_mec2(M,840,nmb,21,1500); a2=md(R2);
fprintf('\n  RETENU : bouche maillee nm=%d\n',nmb);
row=@(n,a,f) fprintf('    %-22s %9.4f %9.4f %+8.1f %%\n',n,a,f,(a-f)/f*100);
row('Bg1 fondamental (T)',a2(1),aF(1));
row('B moyen |Br| (T)',mean(abs(R2.Br)),mean(abs(BrFu)));
row('B crete (T)',max(R2.Br),max(BrFu));
row('Bt RMS (T)',R2.Btrms,BtF_rms);
row('harmonique 8 (T)',a2(2),aF(2));
row('harmonique 22 (T)',a2(3),aF(3));
row('harmonique 23 (T)',a2(4),aF(4));
row('harmonique 37 (T)',a2(5),aF(5));
fprintf('    erreur denture : %.1f %% (scalaire identifie : %.1f %%)\n',err(a2),err(a1));

figure('Name','Bouche d''encoche 2D','Color','w','Position',[60 60 1100 420]);
subplot(1,2,1);
kk=1:40;
sF=arrayfun(@(k) amp(BrFu,thu,k),kk);
s1=arrayfun(@(k) amp(R1.Br(1:end-1),R1.thq(1:end-1),k),kk);
s2=arrayfun(@(k) amp(R2.Br(1:end-1),R2.thq(1:end-1),k),kk);
stem(kk-0.25,sF,'r','filled','MarkerSize',3.5); hold on;
stem(kk,s1,'b','filled','MarkerSize',3.5);
stem(kk+0.25,s2,'g','filled','MarkerSize',3.5);
set(gca,'YScale','log'); grid on; xlabel('ordre spatial'); ylabel('|B_r| (T)');
legend('FEA','kfringe identifie','bouche 2D','Location','best');
title('Spectre du champ d''entrefer');
subplot(1,2,2);
plot(thu*180/pi,BrFu,'r--','LineWidth',1.1); hold on;
plot(R1.thq*180/pi,R1.Br,'b','LineWidth',1.0);
plot(R2.thq*180/pi,R2.Br,'g','LineWidth',1.3); grid on; xlim([0 50]);
xlabel('\theta mec (deg)'); ylabel('B_r (T)');
legend('FEA','kfringe','bouche 2D','Location','south'); title('Profil (zoom 2 poles)');
saveas(gcf,fullfile(fileparts(mfilename('fullpath')),'FIG_slot2d.png'));
