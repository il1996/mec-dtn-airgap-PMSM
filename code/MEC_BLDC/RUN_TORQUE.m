%% RUN_TORQUE  Couple en charge : courant six-step IDEAL vs courant REEL FEA
%  Le couple MEC etait calcule avec un creneau 120 deg parfait, en phase avec
%  la back-EMF. Le courant reel commute progressivement (inductance 50 mH) :
%  il est arrondi et dephasÃ©. On refait donc le calcul avec le COURANT
%  MESURE par la FEA, pour separer l'erreur de MODELE de l'erreur
%  d'IDEALISATION DU COURANT.
%  Le recalage de phase n'utilise QUE des donnees FEA : la phase du
%  fondamental du flux totalise FL_a fixe la correspondance temps <-> position.
clear; clc;
M=machine_bldc(); Ns=M.Ns; p=M.p; Nm=M.Nm;
if isfile('kfringe_ident.mat'), S=load('kfringe_ident.mat'); kfr=S.kbest; else, kfr=0.5; end
PA=[1 -2 -15 3 14]; PB=[6 -7 -5 8 4]; PC=[11 -12 -10 13 9];
fea=fullfile(M.FEA.dir,'transitoire (en charge)');

%% --- MEC : back-EMF sur une periode electrique --------------------------
Np=721; R=cogging_mec(M,840,0,Np,M.muI,kfr,2*pi/p);
lam=@(P) M.Ntc*sum(sign(P(:)).*R.PhiT(abs(P(:)),:),1);
om=M.speed*2*pi/60; phi=R.phis;
lamA=lam(PA);
eA=-gradient(lamA,phi)*om; eB=-gradient(lam(PB),phi)*om; eC=-gradient(lam(PC),phi)*om;
ang=p*phi;                                  % angle electrique du modele

%% --- FEA : courants et flux en regime etabli ---------------------------
di=readmatrix(fullfile(fea,'BranchCurrent Plot 1.tab'),'FileType','text','NumHeaderLines',1);
dl=readmatrix(fullfile(fea,'Output Variables Plot 1.tab'),'FileType','text','NumHeaderLines',1);
ds=readmatrix(fullfile(fea,'Speed Plot 1.tab'),'FileType','text','NumHeaderLines',1);
n=size(di,1); w=round(n*0.6):n;              % regime etabli
t=di(w,1)*1e-3;                              % s
ia=di(w,2); ib=di(w,3); ic=di(w,4);
nl=size(dl,1); wl=round(nl*0.6):nl;
tl=dl(wl,1)*1e-3; FLa=dl(wl,2);              % col2 = FL_a (col3 = FL_c, col4 = FL_b)
nsp=size(ds,1); nFEA=mean(ds(round(nsp*0.75):nsp,2));
fe=nFEA*Nm/120;
fprintf('=== COUPLE EN CHARGE : courant ideal vs courant reel FEA ===\n');
fprintf('  Essai FEA en charge : %.1f tr/min -> f_elec = %.2f Hz\n',nFEA,fe);
fprintf('  I_phase FEA : %.3f A eff, %.3f A crete\n',sqrt(mean(ia.^2)),max(abs(ia)));

% phase du fondamental (recalage : uniquement des donnees FEA)
zF=mean(FLa.*exp(-1i*2*pi*fe*tl));  psiF=angle(zF);
zM=mean(lamA.*exp(-1i*ang));        psiM=angle(zM);
the=2*pi*fe*t - psiF + psiM;                 % angle electrique de chaque echantillon

% back-EMF du MODELE aux positions correspondantes
% interpolation PERIODIQUE (ang va de 0 a 2*pi inclus -> retirer le doublon)
angp=ang(1:end-1);
ip=@(y,x) interp1([angp 2*pi],[y(1:end-1) y(1)],mod(x,2*pi),'linear');
eAi=ip(eA,the); eBi=ip(eB,the); eCi=ip(eC,the);

T_real = mean(eAi.*ia + eBi.*ib + eCi.*ic)/om;

%% --- couple avec le creneau ideal (pour memoire) ------------------------
Iflat=M.Iph_rms_load*sqrt(3/2);
ss=@(e) sixstep(e,Iflat);
T_ideal=mean((eA.*ss(eA)+eB.*ss(eB)+eC.*ss(eC))/om);

fprintf('\n  %-42s %8s %9s\n','','T (N.m)','ecart');
r=@(nom,v) fprintf('  %-42s %8.3f %+8.1f %%\n',nom,v,(v-M.FEA.T_load)/M.FEA.T_load*100);
r('MEC, creneau 120 deg IDEAL',T_ideal);
r('MEC, COURANT REEL de la FEA',T_real);
fprintf('  %-42s %8.3f\n','FEA (P_em/omega)',M.FEA.T_load);
fprintf('\n  -> l''idealisation du courant explique %.1f points sur %.1f.\n',...
        (T_ideal-T_real)/M.FEA.T_load*100,(T_ideal-M.FEA.T_load)/M.FEA.T_load*100);

figure('Name','Couple : courant reel vs ideal','Color','w','Position',[70 70 1000 400]);
subplot(1,2,1);
[ths,is]=sort(mod(the,2*pi));
plot(ths*180/pi,ia(is),'r','LineWidth',1.2); hold on;
plot(ang*180/pi,ss(eA),'b--','LineWidth',1.2); grid on;
xlabel('angle electrique (deg)'); ylabel('i_A (A)');
legend('FEA (reel)','creneau ideal','Location','best'); title('Courant de phase');
subplot(1,2,2);
bar([T_ideal T_real M.FEA.T_load]); grid on;
set(gca,'XTickLabel',{'ideal','courant reel','FEA'}); ylabel('couple (N.m)');
title(sprintf('%.2f / %.2f / %.2f N.m',T_ideal,T_real,M.FEA.T_load));
saveas(gcf,fullfile(fileparts(mfilename('fullpath')),'FIG_torque.png'));

function i=sixstep(e,Iflat)
    n=numel(e); k=(0:n-1)/n*2*pi;
    psi=atan2(2*mean(e.*sin(k)),2*mean(e.*cos(k)));
    ke=mod(k-psi,2*pi); i=zeros(size(e));
    i(ke<pi/3|ke>5*pi/3)=Iflat; i(ke>2*pi/3&ke<4*pi/3)=-Iflat;
end
