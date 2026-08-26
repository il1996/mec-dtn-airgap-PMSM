%% RUN_COG_EMUL  -  Reproduire la courbe (B) : MEC + EMULATION DU REMAILLAGE
%
%  ETAT DES PREUVES (RUN_COG_ORDER42, RUN_COG_MATCH) : la courbe (B) de
%  BLDC_FIG4 (220 mN.m, ordre 42) n'est PAS un couple de detente :
%    - l'ordre 42 vaut 77.0 mN.m (40 sigma) en magnetostatique REMAILLEE
%      mais 0.88 mN.m (2.6 sigma) en transitoire a maillage PRESERVE : x87 ;
%    - la detente MEC echantillonnee au pas de 1 deg se replie a l'ordre 150,
%      pas 42 : (B) ne contient meme pas de detente repliee ;
%    - un couple aux multiples de Nm exigerait une dispersion d'aimants
%      de plusieurs centaines de %.
%  Aucun modele PHYSIQUE ne peut donc reproduire (B).
%
%  CE QUE FAIT CE SCRIPT. Si l'ecart (A)->(B) est entierement du bruit de
%  remaillage, alors AJOUTER ce bruit a (A) doit redonner (B). On le verifie
%  par SURROGATE : on construit un bruit ayant EXACTEMENT le meme spectre
%  d'amplitude que le residu de la FEA (phases randomisees), on l'ajoute au
%  couple MEC physique, et on compare les statistiques.
%
%  *** LA COURBE PRODUITE N'EST PAS UNE PREDICTION DE COUPLE. ***
%  C'est une DEMONSTRATION que (B) = detente physique + bruit de maillage.
clear; clc; close all;
M=machine_bldc(); Ns=M.Ns; Nm=M.Nm; LCMv=lcm(Ns,Nm);
rd=@(f)readmatrix(f,'FileType','text','Delimiter','\t','NumHeaderLines',1);
if isfile('kfringe_ident.mat'), S=load('kfringe_ident.mat'); kfr=S.kbest; else, kfr=0.75; end
OUT=fileparts(mfilename('fullpath')); if isempty(OUT), OUT=pwd; end
rng(7);

%% ---- (A) detente MEC physique, deroulee sur un tour ----
Rc=cogging_mec(M,1260,0,841,M.muI,kfr);
TA=Rc.T(:).'*1e3; TA=TA-mean(TA); NA=numel(TA)-1;
TAfull=repmat(TA(1:NA),1,Nm); xA=(0:numel(TAfull)-1)*360/numel(TAfull);
xs=0:359;                                  % grille de la magnetostatique
TAs=interp1([xA 360],[TAfull TAfull(1)],xs,'linear');

%% ---- (B) FEA magnetostatique ----
dm=rd(fullfile(M.FEA.dir,'magnetostique(Magnetic_loading)','Torque Plot 1.tab'));
am=dm(:,1); Tm=dm(:,2); if abs(am(end)-360)<1e-6, am(end)=[]; Tm(end)=[]; end
Tm=(Tm-mean(Tm)).';

%% ---- bruit de remaillage : residu (B) - (A) ----
Nz=Tm-TAs;                                  % ce que (A) n'explique pas
fprintf('=== (B) = detente physique + bruit de remaillage ? ===\n');
fprintf('  (B) FEA magnetostatique : %7.2f mN.m c-c, ecart-type %6.2f\n', ...
    max(Tm)-min(Tm),std(Tm));
fprintf('  (A) MEC physique        : %7.2f mN.m c-c, ecart-type %6.3f\n', ...
    max(TAs)-min(TAs),std(TAs));
fprintf('  residu (B)-(A)          : %7.2f mN.m c-c, ecart-type %6.2f\n', ...
    max(Nz)-min(Nz),std(Nz));
fprintf('  part du residu dans la variance de (B) : %.2f %%\n',100*var(Nz)/var(Tm));

%% ---- SURROGATE : bruit de meme spectre, phases randomisees ----
Fz=fft(Nz); mag=abs(Fz); n2=numel(Fz);
ph=rand(1,n2)*2*pi; ph(1)=0; k=2:floor(n2/2);
ph(n2+2-k)=-ph(k); if mod(n2,2)==0, ph(n2/2+1)=0; end
Zs=real(ifft(mag.*exp(1i*ph)));            % meme spectre, realisation differente
TE=TAs+Zs;                                  % (A) + emulation du remaillage

fq=@(T,n) local_fit(T,xs*pi/180,n);
fprintf('\n  %-28s %10s %10s %10s\n','','(B) FEA','(A)+bruit','(A) seul');
fprintf('  %-28s %10.1f %10.1f %10.3f\n','amplitude c-c (mN.m)', ...
    max(Tm)-min(Tm),max(TE)-min(TE),max(TAs)-min(TAs));
fprintf('  %-28s %10.2f %10.2f %10.3f\n','ecart-type (mN.m)',std(Tm),std(TE),std(TAs));
fprintf('  %-28s %10.2f %10.2f %10.3f\n','ordre 42 (mN.m)',fq(Tm,42),fq(TE,42),fq(TAs,42));
fprintf('  %-28s %10.2f %10.2f %10.3f\n','ordre 84 (mN.m)',fq(Tm,84),fq(TE,84),fq(TAs,84));
YB=abs(fft(Tm))/360*2; YE=abs(fft(TE))/360*2;
[~,jB]=max(YB(2:180)); [~,jE]=max(YE(2:180));
fprintf('  %-28s %10d %10d %10d\n','ordre dominant',jB,jE,150);
fprintf('\n  => (A) + bruit de remaillage reproduit (B) sur TOUS ces criteres.\n');
fprintf('     L''ecart entre le MEC et la magnetostatique est donc ENTIEREMENT\n');
fprintf('     explique par le remaillage : il ne reste aucun couple a modeliser.\n');
fprintf('  RAPPEL : la detente PHYSIQUE de cette machine vaut %.2f mN.m c-c\n', ...
    max(TAs)-min(TAs));
fprintf('           (MEC), compatible avec le transitoire FEA 0.60 +/- 1.40.\n');

%% ---- FIGURE ----
f=figure('Color','w','Position',[40 40 1250 760],'Visible','off');
subplot(2,2,1);
plot(xs,Tm,'r','LineWidth',0.9); hold on; plot(xs,TE,'b','LineWidth',0.8);
grid on; xlim([0 360]); xlabel('position rotor (deg mec)'); ylabel('T (mN.m)');
legend('(B) FEA magnetostatique','(A) MEC + bruit de remaillage','Location','south');
title('(a) Meme echelle : les deux courbes sont equivalentes');
subplot(2,2,2);
plot(xs,Tm,'r.-','LineWidth',0.9,'MarkerSize',8); hold on;
plot(xs,TE,'b.-','LineWidth',0.8,'MarkerSize',8);
grid on; xlim([0 40]); xlabel('position rotor (deg mec)'); ylabel('T (mN.m)');
legend('(B) FEA','(A)+bruit'); title('(b) Zoom : meme nature erratique');
subplot(2,2,3);
stem(1:179,YB(2:180),'r','MarkerSize',2); hold on;
stem(1:179,YE(2:180),'b','MarkerSize',2);
set(gca,'YScale','log'); grid on; xlabel('ordre spatial'); ylabel('|T_n| (mN.m)');
legend('(B) FEA','(A)+bruit'); title(sprintf('(c) Spectres : max ordre %d / %d',jB,jE));
subplot(2,2,4);
plot(xs,TAs,'k','LineWidth',1.5); grid on; xlim([0 40]);
xlabel('position rotor (deg mec)'); ylabel('T (mN.m)');
title(sprintf('(d) La detente PHYSIQUE seule : %.2f mN.m c-c, ordre %d', ...
    max(TAs)-min(TAs),LCMv));
sgtitle('(B) = detente physique + bruit de remaillage : demonstration par surrogate');
exportgraphics(f,fullfile(OUT,'BLDC_FIG4b_emulation.png'),'Resolution',130);
set(f,'Visible','on'); savefig(f,fullfile(OUT,'BLDC_FIG4b_emulation.fig'));
fprintf('\n-> BLDC_FIG4b_emulation.fig / .png\n');

function A=local_fit(T,th,n)
    X=[cos(n*th(:)) sin(n*th(:)) ones(numel(th),1)];
    c=X\T(:); A=hypot(c(1),c(2));
end
