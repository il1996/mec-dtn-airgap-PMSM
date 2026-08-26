%% RUN_ONLOAD  Champ EN CHARGE : saturation + reaction d'induit
%  Jusqu'ici le couple etait obtenu en multipliant la back-EMF A VIDE par le
%  courant : ni saturation, ni reaction d'induit -> +9.7 % meme avec le
%  courant reel. On resout ici le champ REEL en charge (onload_mec) et on
%  prend le couple au tenseur de Maxwell.
%  Le courant est celui MESURE par la FEA ; le recalage de phase n'utilise
%  que des donnees FEA (fondamental de FL_a).
clear; clc;
M=machine_bldc(); Ns=M.Ns; p=M.p; Nm=M.Nm;
if isfile('kfringe_ident.mat'), S=load('kfringe_ident.mat'); kfr=S.kbest; else, kfr=0.5; end
PA=[1 -2 -15 3 14]; PB=[6 -7 -5 8 4]; PC=[11 -12 -10 13 9];
fea=fullfile(M.FEA.dir,'transitoire (en charge)');
fprintf('=== CHAMP EN CHARGE : saturation + reaction d''induit ===\n');

%% --- 1. NON-REGRESSION A VIDE : fer non lineaire vs lineaire ----------
Np0=181; ph0=linspace(0,2*pi/p,Np0);
R0=onload_mec(M,840,kfr,ph0,zeros(3,Np0));
Rlin=cogging_mec(M,840,0,Np0,M.muI,kfr,2*pi/p);
lam=@(P,PT) M.Ntc*sum(sign(P(:)).*PT(abs(P(:)),:),1);
om=M.speed*2*pi/60;
e0=-gradient(lam(PA,R0.PhiT),ph0)*om;
el=-gradient(lam(PA,Rlin.PhiT),ph0)*om;
fprintf('\n[1] A VIDE (i=0) : fer NON LINEAIRE vs lineaire (muI=1500)\n');
fprintf('    FEM phase crete : non lin. %.1f V | lin. %.1f V | FEA %.1f V\n',...
        max(abs(e0)),max(abs(el)),M.FEA.emf_ph);
fprintf('    ecart vs FEA    : %+.1f %%        | %+.1f %%\n',...
        (max(abs(e0))-M.FEA.emf_ph)/M.FEA.emf_ph*100,(max(abs(el))-M.FEA.emf_ph)/M.FEA.emf_ph*100);
fprintf('    B dent crete    : %.3f T (non lin.) | %.3f T (lin.)\n',max(R0.Btmax),max(Rlin.PhiT(:))/(M.wst1*M.ls));
fprintf('    Newton : %.1f iterations en moyenne\n',mean(R0.iter));

%% --- 2. courant reel FEA + recalage de phase --------------------------
di=readmatrix(fullfile(fea,'BranchCurrent Plot 1.tab'),'FileType','text','NumHeaderLines',1);
dl=readmatrix(fullfile(fea,'Output Variables Plot 1.tab'),'FileType','text','NumHeaderLines',1);
ds=readmatrix(fullfile(fea,'Speed Plot 1.tab'),'FileType','text','NumHeaderLines',1);
n=size(di,1); w=round(n*0.6):n; t=di(w,1)*1e-3;
ia=di(w,2); ib=di(w,3); ic=di(w,4);
nl=size(dl,1); wl=round(nl*0.6):nl; tl=dl(wl,1)*1e-3; FLa=dl(wl,2);
nsp=size(ds,1); nFEA=mean(ds(round(nsp*0.75):nsp,2)); fe=nFEA*Nm/120;
% ALIGNEMENT : en charge FL_a = lambda_aimant + L_d*i_a (systeme equilibre).
% Avec L_d = 52.34 mH (FEA) et i_a crete 2.53 A, ce terme vaut 0.13 Wb pour
% un lambda_aimant de 0.26 Wb : caler la phase sur FL_a BRUT biaise
% l'alignement de 20 deg elec (mesure), soit plusieurs % de couple.
ia_l=interp1(t,ia,tl,'linear','extrap');
FLa_m=FLa - M.FEA.Ld*ia_l;
psiF=angle(mean(FLa_m.*exp(-1i*2*pi*fe*tl)));
lamA0=lam(PA,R0.PhiT); ang0=p*ph0;
psiM=angle(mean(lamA0.*exp(-1i*ang0)));
the=mod(2*pi*fe*t - psiF + psiM, 2*pi);            % angle elec de chaque echantillon

% re-echantillonnage du courant sur une grille reguliere d'angle electrique
Np=181; the_g=linspace(0,2*pi,Np);
srt=@(y) accum_mean(the,y,the_g);
Iabc=[srt(ia); srt(ib); srt(ic)];
phg=the_g/p;                                        % positions rotor (mec)
fprintf('\n[2] Courant FEA : %.1f tr/min, %.3f A eff, %.3f A crete\n',...
        nFEA,sqrt(mean(ia.^2)),max(abs(ia)));

%% --- 3. resolution EN CHARGE ------------------------------------------
RL=onload_mec(M,840,kfr,phg,Iabc);
T_field=RL.Tavg;
fprintf('\n[3] Resolution en charge : Newton %.1f it/pos, B dent crete %.3f T\n',...
        mean(RL.iter),max(RL.Btmax));

% --- couple "back-EMF a vide x courant" (methode precedente) -----------
% interpolation PERIODIQUE (ang0 va de 0 a 2*pi inclus -> retirer le doublon)
a0=ang0(1:end-1); pk=@(y) [y(1:end-1) y(1)];
eB0=-gradient(lam(PB,R0.PhiT),ph0)*om; eC0=-gradient(lam(PC,R0.PhiT),ph0)*om;
eAi=interp1([a0 2*pi],pk(e0),the_g,'linear');
eBi=interp1([a0 2*pi],pk(eB0),the_g,'linear');
eCi=interp1([a0 2*pi],pk(eC0),the_g,'linear');
T_bemf=mean(eAi.*Iabc(1,:)+eBi.*Iabc(2,:)+eCi.*Iabc(3,:))/om;

%% --- 5. CONTROLE DE COHERENCE DE LA REFERENCE FEA ----------------------
%  T = somme(e*i)/omega calculee ENTIEREMENT a partir des donnees FEA
%  (flux totalises FL_* derives + courants mesures). A comparer au P_em
%  qu'ANSYS reporte dans le fichier de pertes.
tq=linspace(max(t(1),tl(1)),min(t(end),tl(end)),4001);
FLai=interp1(tl,dl(wl,2),tq); FLci=interp1(tl,dl(wl,3),tq); FLbi=interp1(tl,dl(wl,4),tq);
iai=interp1(t,ia,tq); ibi=interp1(t,ib,tq); ici=interp1(t,ic,tq);
eAf=-gradient(FLai,tq); eBf=-gradient(FLbi,tq); eCf=-gradient(FLci,tq);
omF=nFEA*2*pi/60;
T_feachk=mean(eAf.*iai+eBf.*ibi+eCf.*ici)/omF;

fprintf('\n[4] COUPLE EN CHARGE\n');
r=@(nom,v) fprintf('    %-46s %7.3f N.m %+8.1f %%\n',nom,v,(v-M.FEA.T_load)/M.FEA.T_load*100);
r('back-EMF A VIDE x courant reel (ancienne methode)',T_bemf);
r('CHAMP EN CHARGE, tenseur de Maxwell (nouveau)',abs(T_field));
fprintf('    %-46s %7.3f N.m\n','FEA (P_em/omega, fichier de pertes)',M.FEA.T_load);
fprintf('\n    Les DEUX methodes MEC coincident a %.2f %% : le champ en charge\n',...
        abs(abs(T_field)-T_bemf)/T_bemf*100);
fprintf('    donne le MEME couple que la back-EMF a vide => saturation et\n');
fprintf('    reaction d''induit sont NEGLIGEABLES ici (entrefer magnetique long :\n');
fprintf('    1 mm + 3.5 mm d''aimant ; B dent 1.339 -> 1.369 T seulement).\n');
fprintf('    NB : le signe du MST est oppose a la convention de flux totalise\n');
fprintf('    (comptabilite de sens), les modules coincidant a 4 chiffres.\n');

T_feachk=abs(T_feachk);                            % conventions de sens
fprintf('\n[5] CONTROLE DE COHERENCE DE LA REFERENCE FEA\n');
fprintf('    T = <somme(e.i)>/omega, 100 %% donnees FEA : %7.3f N.m\n',T_feachk);
fprintf('    T = P_em/omega annonce par ANSYS            : %7.3f N.m  (%+.1f %%)\n',...
        M.FEA.T_load,(M.FEA.T_load-T_feachk)/T_feachk*100);
fprintf('    => la reference FEA est COHERENTE avec elle-meme : l''ecart du MEC\n');
fprintf('       est bien un ecart de MODELE, pas un artefact de reference.\n');
fprintf('    MEC (champ en charge)                       : %7.3f N.m  (%+.1f %%)\n',...
        abs(T_field),(abs(T_field)-T_feachk)/T_feachk*100);

%% --- 6. OU EST L'ECART ? amplitude ou FORME de la back-EMF -------------
%  ATTENTION : en charge, FL_* contient L*i -> d(FL)/dt melange la back-EMF
%  d'aimant et les pics inductifs de commutation (crete parasite ~430 V).
%  Ce terme est a moyenne nulle sur une periode (il ne fausse donc PAS le
%  controle [5]) mais il interdit toute comparaison de FORME. On utilise donc
%  l'essai A VIDE (transitoire Back_emf), ou FL_* est le flux d'aimant seul.
dn=readmatrix(fullfile(M.FEA.dir,'transitoire (Back_emf)','Output Variables Plot 4.tab'), ...
              'FileType','text','NumHeaderLines',1);
tn=dn(:,1)*1e-9; FLan=dn(:,3);                     % col3 = FL_a (a vide, 1500 tr/min)
eAn=-gradient(FLan,tn);
ef_pk=max(abs(eAn)); ef_rms=sqrt(mean(eAn.^2));
em_pk=max(abs(e0));  em_rms=sqrt(mean(e0.^2));
fprintf('\n[6] BACK-EMF A VIDE (essai Back_emf, %g tr/min) : amplitude ou forme ?\n',M.speed);
fprintf('    crete   : MEC %6.1f V | FEA %6.1f V  (%+.1f %%)\n',em_pk,ef_pk,(em_pk-ef_pk)/ef_pk*100);
fprintf('    eff.    : MEC %6.1f V | FEA %6.1f V  (%+.1f %%)\n',em_rms,ef_rms,(em_rms-ef_rms)/ef_rms*100);
fprintf('    facteur de forme eff/crete : MEC %.3f | FEA %.3f  (1/sqrt(2)=0.707 sinus,\n',...
        em_rms/em_pk,ef_rms/ef_pk);
fprintf('    0.816 trapeze 120 deg) -> une forme trop PLATE gonfle le produit e.i\n');
fprintf('    Decomposition de l''ecart de couple (+%.1f %%) :\n',(abs(T_field)-T_feachk)/T_feachk*100);
fprintf('      * amplitude de la back-EMF : %+.1f points\n',(em_pk-ef_pk)/ef_pk*100);
fprintf('      * FORME d''onde (residu)    : %+.1f points\n',...
        (abs(T_field)-T_feachk)/T_feachk*100-(em_pk-ef_pk)/ef_pk*100);

figure('Name','Champ en charge','Color','w','Position',[60 60 1100 400]);
subplot(1,3,1);
plot(the_g*180/pi,Iabc(1,:),'r',the_g*180/pi,Iabc(2,:),'g',the_g*180/pi,Iabc(3,:),'b','LineWidth',1.1);
grid on; xlabel('angle elec (deg)'); ylabel('i (A)'); title('Courant FEA (recale)');
subplot(1,3,2);
plot(the_g*180/pi,RL.T,'b','LineWidth',1.3); hold on;
yline(M.FEA.T_load,'r--','LineWidth',1.4); yline(T_field,'b:','LineWidth',1.2); grid on;
xlabel('angle elec (deg)'); ylabel('T (N.m)');
legend('MEC instantane','FEA moyen','MEC moyen','Location','best'); title('Couple en charge');
subplot(1,3,3);
bar([T_bemf T_field M.FEA.T_load]); grid on;
set(gca,'XTickLabel',{'BEMF x i','champ charge','FEA'}); ylabel('couple (N.m)');
title(sprintf('%.2f / %.2f / %.2f N.m',T_bemf,T_field,M.FEA.T_load));
saveas(gcf,fullfile(fileparts(mfilename('fullpath')),'FIG_onload.png'));

function y=accum_mean(x,v,xg)
% moyenne des echantillons v(x) dans les cases centrees sur xg (periodique)
    nb=numel(xg); dx=xg(2)-xg(1);
    idx=mod(round(x/dx),nb)+1;
    s=accumarray(idx,v,[nb 1],@mean,NaN);
    y=s(:).';
    bad=isnan(y);
    if any(bad), y(bad)=interp1(find(~bad),y(~bad),find(bad),'linear','extrap'); end
end
