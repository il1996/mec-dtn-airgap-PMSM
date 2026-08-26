%% RUN_R6B_FEA - R6, volet reference : l'inductance incrementale mesuree sur l'EF
%
%  CE QUE R6 DEMANDE. Comparer L_ll(i) du reseau a l'inductance INCREMENTALE
%  dLambda/di extraite de la reference elements finis, a 2, 5, 10, 15 et 20 A.
%  Le piege nomme par la specification -- confondre l'incrementale et la
%  secante lambda/i -- est traite frontalement : les deux sont calculees et
%  affichees cote a cote, sur les memes donnees.
%
%  OU EST LA MESURE. Le projet EF ne contient AUCUN balayage en courant :
%  magnetostique(Armature-Field) ne donne qu'un point, i_a = 1 A. La seule
%  source de lambda(i) est la MONTEE DE COURANT du transitoire en charge,
%  ou le rotor est encore quasi immobile.
%
%  DEUX ROUTES INDEPENDANTES, et c'est la garde.
%    route 1 : L = (dLambda - (dLambda/dtheta_e)*dtheta_e) / di
%              lue sur les FLUX TOTALISES exportes ;
%    route 2 : L = (V_applique - 2*R*i - e_ligne) / (di/dt)
%              lue sur les TENSIONS DE NOEUD exportees.
%  Elles ne partagent aucune colonne de donnees.
%
%  DEUX PIEGES DE MISE EN OEUVRE, CORRIGES ICI APRES LES AVOIR SUBIS.
%   (a) di/dt doit etre forme A L'INTERIEUR de la fenetre. Le forme sur le
%       vecteur complet fait entrer l'echantillon de COMMUTATION dans le
%       dernier point, ou i_c retombe de 12,27 a 8,95 A : la route 2 y
%       sortait -60 mH.
%   (b) local_Lline MOYENNE sur les 61 positions rotor, alors que la fenetre
%       EF est a UNE position (excursion 0,34 deg elec). Les deux nombres ne
%       sont pas la meme grandeur. On evalue donc le reseau DES DEUX FACONS.
%
%  LE TERME DE ROTATION EST SOUSTRAIT, PAS SEULEMENT MAJORE, avec un
%  dLambda/dtheta relu de l'essai a vide a la position reconstruite.
%
%  DIAGNOSTIC DECISIF (section 7). La correction de double comptage de
%  mec_map.m:203-231 cale L(0) sur inductance_mec(...,0), reseau LINEAIRE et
%  AIMANTS ETEINTS. Si les aimants pre-saturent la denture, cette calibration
%  detruit l'effet -- a courant nul ET, la soustraction etant lineaire en i,
%  a TOUS les courants. On le teste en rejouant la carte a kfringe = 0, ou il
%  n'y a plus de double comptage a corriger : ce qui reste de l'ecart entre
%  la valeur BRUTE et la cible lineaire est l'effet des aimants.
%
%  AUCUN CHIFFRE TRANSCRIT.
clear; clc; t0=tic;
diary('R6_satmap_out.txt'); diary on;
M=machine_bldc(); p=M.p;
KFR_A=0.75; NSURF=720; NTH=61; NA=43;
FEA=M.FEA.dir;
rd=@(f)readmatrix(f,'FileType','text','Delimiter','\t','NumHeaderLines',1);

fprintf('=== R6 : carte de saturation, inductance incrementale contre EF ===\n\n');
fprintf('  CONFIGURATION\n');
fprintf('    machine        : PMSM 15/14, 750 W (machine_bldc)\n');
fprintf('    chaine MEC     : mec_map + local_Lline (recopiee verbatim)\n');
fprintf('    Nsurf / Nth    : %d / %d\n',NSURF,NTH);
fprintf('    grille courant : na = %d (production = 19 ; convergence en R6a)\n',NA);
fprintf('    k_fringe       : %.4f (role induit) ET 0 (diagnostic §7)\n',KFR_A);
fprintf('    solveur        : Newton, B(H) M350-50A\n');
fprintf('    reference EF   : %s\n',FEA);
fprintf('    cache          : AUCUN\n');

%% ================= 1. REFERENCE EF, DONNEES BRUTES ====================
LD=fullfile(FEA,'transitoire (en charge)');
FL=rd(fullfile(LD,'Output Variables Plot 1.tab'));   % t, FL_a, FL_c, FL_b
IB=rd(fullfile(LD,'BranchCurrent Plot 1.tab'));      % t, i_a, i_b, i_c
SP=rd(fullfile(LD,'Speed Plot 1.tab'));              % t, n
NV=rd(fullfile(LD,'NodeVoltage Plot 2.tab'));        % t[ns], Ie_abc, Iv_abc [mV]
NL=rd(fullfile(FEA,'transitoire (Back_emf)','Output Variables Plot 4.tab'));
AR=rd(fullfile(FEA,'magnetostique(Armature-Field)','Output Variables Table 1.tab'));
tFL=FL(:,1); tIB=IB(:,1); tSP=SP(:,1); tNV=NV(:,1)*1e-6;

fprintf('\n  ---- 1. DONNEES EF ----\n');
fprintf('    flux %d pts | courants %d pts | tensions %d pts\n', ...
    numel(tFL),numel(tIB),numel(tNV));
fprintf('    le fichier de courant n''a pas la ligne t = 0 : appariement PAR TEMPS.\n');
iNV=interp1(tIB,IB(:,2:4),tNV,'linear','extrap');
Ie=NV(:,2:4)*1e-3; Iv=NV(:,5:7)*1e-3; msk=abs(iNV)>1;
Rmes=mean((Iv(msk)-Ie(msk))./iNV(msk));
fprintf('    R_ph MESUREE sur la reference : %.6f ohm (code : %.6f, ecart %+.3f %%)\n', ...
    Rmes,M.Rph,100*(M.Rph-Rmes)/Rmes);

%% ---- fenetre W1 : conduction c->b, rotor quasi immobile ---------------
ia=IB(:,2); ib=IB(:,3); ic=IB(:,4);
w=find(abs(ia)<1e-3 & ic>0 & tIB<1.3);
q0=ic(w)>0.1;
fprintf('\n  ---- 2. FENETRE EXPLOITABLE ----\n');
fprintf('    conduction c -> b : %d points, t = %.4f -> %.4f ms\n', ...
    numel(w),tIB(w(1)),tIB(w(end)));
fprintf('    ecart max sur i_c = -i_b (hors point a courant nul) : %.1f ppm\n', ...
    1e6*max(abs(ic(w(q0))+ib(w(q0)))./ic(w(q0))));
fprintf('    courant maximal dans la fenetre : %.4f A\n',max(ic(w)));
fprintf('    echantillon suivant (commutation) : i_c = %.4f A -> EXCLU\n',ic(w(end)+1));

om_m=SP(:,2)*2*pi/60; th_m=cumtrapz(tSP*1e-3,om_m);
the_ib=p*interp1(tSP,th_m,tIB,'linear','extrap');
om_ib=p*interp1(tSP,om_m,tIB,'linear','extrap');
fprintf('    excursion electrique sur la fenetre : %.4f deg\n', ...
    (max(the_ib(w))-min(the_ib(w)))*180/pi);

the_nl=NL(:,2)*p*pi/180; lcb_nl=NL(:,5)-NL(:,7);
[the_s,is]=sort(mod(the_nl,2*pi)); lcb_s=lcb_nl(is);
dl_dthe=gradient(lcb_s,the_s);
dldt_ib=interp1(the_s,dl_dthe,mod(the_ib,2*pi),'linear','extrap');
fprintf('    |dLambda_cb/dtheta_e| max a vide : %.4f Wb/rad elec\n',max(abs(dl_dthe)));

%% ---- les deux routes, formees DANS la fenetre -------------------------
lcb=interp1(tFL,FL(:,3)-FL(:,4),tIB,'linear','extrap');
iW=ic(w); lW=lcb(w); thW=the_ib(w); dlW=dldt_ib(w);
Vw=interp1(tNV,Iv(:,3)-Iv(:,2),tIB(w),'linear','extrap');
omW=om_ib(w);
dt_s=(tIB(2)-tIB(1))*1e-3;
n=numel(w); idx=(1:n).';
dlam=gradient(lW,idx); dii=gradient(iW,idx); dthe=gradient(thW,idx);
L1=(dlam-dlW.*dthe)./dii;
L2=(Vw-2*Rmes*iW-omW.*dlW)./(dii/dt_s);
lam_arm=lW-lW(1);
L1sec=lam_arm./max(iW,eps); L1sec(iW<1e-6)=NaN;
rot=100*abs(dlW.*dthe)./max(abs(dlam),eps);

fprintf('\n  ---- 3. REFERENCE EF : DEUX ROUTES INDEPENDANTES ----\n');
fprintf('  %9s %10s %12s %12s %12s %9s\n', ...
    'i (A)','lam_arm','L route1','L route2','L secante','rot. %');
for k=1:n
    fprintf('  %9.4f %10.5f %10.2f mH %9.2f mH %9.2f mH %8.2f\n', ...
        iW(k),lam_arm(k),L1(k)*1e3,L2(k)*1e3,L1sec(k)*1e3,rot(k));
end
%  GARDE sur les points INTERIEURS : aux deux bords les differences sont
%  decentrees et le premier point est a courant nul.
in=2:n-1;
ec=100*abs(L1(in)-L2(in))./abs(L1(in));
fprintf('\n  GARDE 1 : accord des deux routes, points interieurs (%d a %d)\n',in(1),in(end));
fprintf('    ecart moyen %.2f %% | maximal %.2f %% (a i = %.2f A)\n', ...
    mean(ec),max(ec),iW(in(find(ec==max(ec),1))));
G1=max(ec)<15;
if G1, fprintf('    GARDE PASSEE : les deux routes mesurent la meme grandeur.\n');
else,  fprintf('    GARDE ECHOUEE : extraction non fiable, ne rien conclure.\n'); end

%% ================= 4. ANCRAGE A COURANT NUL ===========================
La_F=AR(1,2); Mf=(AR(1,3)+AR(1,4))/2; Lms=2*(La_F-Mf);
Lchord0=(lW(2)-lW(1))/(iW(2)-iW(1));
fprintf('\n  ---- 4. ANCRAGE A COURANT NUL : L''EF SE CONTREDIT ----\n');
fprintf('    EF magnetostatique, i_a = %.0f A, AIMANTS ETEINTS :\n',AR(1,1));
fprintf('      L_a = %.6f H, M = %.6f H -> L_ligne = 2(L_a-M) = %.4f mH\n',La_F,Mf,Lms*1e3);
fprintf('    EF transitoire, AIMANTS PRESENTS, corde [0 ; %.3f A] : %.4f mH\n', ...
    iW(2),Lchord0*1e3);
fprintf('    ecart entre les deux etudes EF : %+.2f %%\n',100*(Lchord0-Lms)/Lms);
fprintf('    A ce courant aucune saturation d''induit n''est possible :\n');
fprintf('      FMM de bobine = Ntc*i = %.1f At. L''ecart ne vient donc pas\n',M.Ntc*iW(2));
fprintf('      du courant, mais de ce qui differe entre les deux etudes.\n');

%% ================= 5. LE RESEAU ========================================
fprintf('\n  ---- 5. RESEAU : reconstruction de la carte ----\n');
uu=linspace(-1,1,NA); iax=28*sign(uu).*abs(uu).^2;
tm=tic; S=mec_map(M,NSURF,KFR_A,NTH,iax);
fprintf('    kfringe = %.3f : %.0f s\n',KFR_A,toc(tm));
%  MOYENNE CIRCULAIRE, et non mean(mod(.,2pi)). La fenetre est a cheval sur
%  zero (de -0,266 a +0,075 deg) : la moyenne arithmetique des angles replies
%  melange des valeurs proches de 0 et de 2*pi et sort 327 deg, une position
%  qui n'a jamais ete occupee. Faute commise, puis corrigee.
theW=mod(angle(mean(exp(1i*thW))),2*pi);
fprintf('    position de la fenetre EF (moyenne circulaire) : theta_e = %.4f rad\n',theW);
fprintf('    (la moyenne arithmetique des angles replies donnerait %.4f rad, faux)\n', ...
    mean(mod(thW,2*pi)));
IQ=[2 5 10];
Lm =arrayfun(@(x)local_Lline(S,x),IQ);          % moyenne sur theta
LmT=arrayfun(@(x)local_Lline_at(S,x,theW),IQ);  % a la position de la fenetre
Lm0=local_Lline(S,0); Lm0T=local_Lline_at(S,0,theW);
dcorr=S.Leq0-2*S.Ldref;
%  dispersion de L sur le tour, pour dire ce que vaut la moyenne
fprintf('\n    sensibilite a la position rotor du nombre publie par le reseau :\n');
fprintf('    %8s %12s %12s %12s %12s\n','i (A)','moyenne','a theta_W','min sur tour','max sur tour');
for k=[0 IQ]
    Lk=arrayfun(@(t)local_Lline_at(S,k,t),S.the);
    if k==0, mk=Lm0; tk=Lm0T; else, j=find(IQ==k,1); mk=Lm(j); tk=LmT(j); end
    fprintf('    %8.0f %10.2f mH %9.2f mH %9.2f mH %9.2f mH\n', ...
        k,mk*1e3,tk*1e3,min(Lk)*1e3,max(Lk)*1e3);
end
fprintf('\n  GARDE 2 : la correction cale-t-elle L(0) sur la cible LINEAIRE ?\n');
fprintf('    L(0) apres correction %.4f mH | cible 2*Ld_lin %.4f mH | ecart %.2e\n', ...
    Lm0*1e3,2*S.Ldref*1e3,abs(Lm0-2*S.Ldref));
G2=abs(Lm0-2*S.Ldref)<1e-9;
if G2, fprintf('    GARDE PASSEE : le calage est exact, par construction.\n');
else,  fprintf('    GARDE ECHOUEE.\n'); end
fprintf('    et cette cible contre l''EF magnetostatique : %+.2f %%\n',100*(Lm0-Lms)/Lms);

%% ================= 6. CONFRONTATION ====================================
fprintf('\n  ---- 6. CONFRONTATION, aux courants ATTEIGNABLES ----\n');
L1q=interp1(iW,L1,IQ); L2q=interp1(iW,L2,IQ); Lsq=interp1(iW(q0),L1sec(q0),IQ);
LFq=0.5*(L1q+L2q);
fprintf('  %6s %12s %12s %12s %12s %11s %11s\n', ...
    'i (A)','EF r1','EF r2','EF moy','MEC <th>','MEC a th_W','ecart');
for k=1:numel(IQ)
    fprintf('  %6.0f %10.2f mH %9.2f mH %9.2f mH %9.2f mH %8.2f mH %9.0f %%\n', ...
        IQ(k),L1q(k)*1e3,L2q(k)*1e3,LFq(k)*1e3,Lm(k)*1e3,LmT(k)*1e3, ...
        100*(Lm(k)-LFq(k))/LFq(k));
end
fprintf('\n  Le piege nomme par la specification, chiffre : comparer la SECANTE\n');
fprintf('  de l''EF a l''INCREMENTALE du reseau donnerait\n');
for k=1:numel(IQ)
    fprintf('    %2d A : %+.0f %% au lieu de %+.0f %%\n', ...
        IQ(k),100*(Lm(k)-Lsq(k))/Lsq(k),100*(Lm(k)-LFq(k))/LFq(k));
end

%% ================= 7. DIAGNOSTIC : QUE DETRUIT LA CALIBRATION ? =======
fprintf('\n  ---- 7. DIAGNOSTIC : l''effet des aimants et son effacement ----\n');
tm=tic; S0=mec_map(M,NSURF,0,NTH,iax);
fprintf('    carte rejouee a kfringe = 0 (aucun double comptage) : %.0f s\n',toc(tm));
fprintf('  %-46s %12s\n','grandeur','mH');
fprintf('  %-46s %12.4f\n','L_ligne BRUTE a i=0, aimants ON, non lineaire',S0.Leq0*1e3);
fprintf('  %-46s %12.4f\n','cible du code : 2*Ld LINEAIRE, aimants OFF',2*S0.Ldref*1e3);
fprintf('  %-46s %12.4f\n','EF magnetostatique, aimants OFF',Lms*1e3);
fprintf('  %-46s %12.4f\n','EF transitoire, aimants ON, i -> 0',Lchord0*1e3);
demag=S0.Leq0-2*S0.Ldref;
fprintf('\n    effet des aimants VU PAR LE RESEAU  : %+.4f mH (%+.1f %%)\n', ...
    demag*1e3,100*demag/(2*S0.Ldref));
fprintf('    effet des aimants VU PAR L''EF       : %+.4f mH (%+.1f %%)\n', ...
    (Lchord0-Lms)*1e3,100*(Lchord0-Lms)/Lms);
fprintf('\n    La correction retranche %+.4f mH a kfringe = %.3f et %+.4f mH a\n', ...
    dcorr*1e3,KFR_A,demag*1e3);
fprintf('    kfringe = 0, de facon LINEAIRE en courant, donc a TOUS les courants.\n');

save('R6_satmap.mat','iW','lW','L1','L2','L1sec','IQ','Lm','LmT','Lm0','LFq', ...
     'Lms','Lchord0','S','S0','Rmes','G1','G2','NA','dcorr','demag','theW');
fprintf('\n  duree totale %.0f s\n=== R6 termine ===\n',toc(t0));
diary off;

% ======================================================================
function L = local_Lline(S,ii)
%LOCAL_LLINE  RECOPIE VERBATIM de BLDC_MEC_COMPLET.m:1277-1283.
    ia=ii; ib=-ii/sqrt(3); s3=sqrt(3);
    q=@(X,k) interpn(1:3,S.the,S.iax,S.iax,X,k,S.the,ia,ib,'linear');
    L=mean((q(S.La,1)-q(S.La,2))-(q(S.Lb,1)-q(S.Lb,2))/s3);
end

function L = local_Lline_at(S,ii,th)
%LOCAL_LLINE_AT  Meme combinaison, a UNE position rotor au lieu de la moyenne.
%   C'est la grandeur homologue de celle que la fenetre EF fournit.
%
%   PERIODISATION OBLIGATOIRE. S.the = linspace(0,2*pi,Nth+1) prive de son
%   dernier point s'arrete a 2*pi - dth = 6.1795 rad. La fenetre EF est a
%   theta_e = 6.2808 rad (soit -0.14 deg), DEHORS : interpn y renvoie NaN.
%   On rajoute donc la tranche theta = 2*pi, egale a celle de theta = 0,
%   exactement comme drive_mec.m:78-88 le fait pour la meme carte.
    ia=ii; ib=-ii/sqrt(3); s3=sqrt(3);
    nt=numel(S.the); na=numel(S.iax); th4=[S.the(:).' 2*pi];
    wrap=@(X,k) cat(1,reshape(X(k,:,:,:),[nt na na]),reshape(X(k,1,:,:),[1 na na]));
    q=@(X,k) interpn(th4,S.iax,S.iax,wrap(X,k),mod(th,2*pi),ia,ib,'linear');
    L=(q(S.La,1)-q(S.La,2))-(q(S.Lb,1)-q(S.Lb,2))/s3;
end
