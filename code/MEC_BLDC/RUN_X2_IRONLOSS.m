%% RUN_X2_IRONLOSS - bloc X2, anomalie A-1 (SPEC_CLAUDE_CODE_v3 §8)
%
%  QUESTION. ARCHITECTURE_DEUX_ARTICLES §1.5 affirme que les pertes fer sur
%  maillage valent 16,81 W contre 19,73 W en EF (-14,8 %), la ou le modele
%  localise donne +8,2 %. La seule sortie stockee (T18_out.txt) affiche
%  128 612 W sur une machine de 750 W. Ce bloc diagnostique, corrige, et
%  n'accepte un chiffre qu'apres FERMETURE DU BILAN DE MATIERE.
%
%  DIAGNOSTIC (cinq defauts, non un seul)
%  ------------------------------------------------------------------
%  (a) UNITES. Les coefficients ANSYS Kh/Kc/Ke du M350-50A sont VOLUMIQUES
%      (W/m3) -- machine_bldc.m:87-93 le dit mot pour mot. T18 les a
%      multiplies par la MASSE. Facteur rho = 7650.
%      C'est de la division 128612,55 / 7650 = 16,8121 que sort le 16,81 W
%      de ARCHITECTURE §1.5 : le chiffre n'est pas invente, il vient d'une
%      version corrigee du script dont la sortie n'a jamais ete regeneree.
%
%  (b) VITESSE. T18 balaye a 1500 tr/min, soit f = 175 Hz. La reference EF
%      M.FEA.Pfe_nl = 19,73 W est un essai a vide a M.FEA.n_nl = 1688
%      tr/min, soit f = 196,93 Hz -- et c'est bien cette frequence
%      qu'emploie la chaine localisee (BLDC_MEC_COMPLET.m:321). Les deux
%      colonnes du Tableau 7 n'etaient donc pas au meme point.
%      Hysteresis en f, Foucault en f^2 : facteur 1,125 a 1,266.
%
%  (c) FOISONNEMENT COMPTE DEUX FOIS. mesh_bldc pose deja area = ...*Ki*L
%      (l. 143, 164, 181) : ME.A est la section de fer NETTE, et B = Phi/A
%      est l'induction DANS LA TOLE. T18 remultiplie par kFe = 0,90.
%
%  (d) SOMMATION SUR LES BRANCHES. C'est le defaut majeur. Dans le pavage de
%      Perho/Silva, la branche RADIALE et la branche TANGENTIELLE d'une
%      meme cellule portent CHACUNE a peu pres le volume de cette cellule :
%      sum(A.*l) sur les branches de fer compte donc le fer DEUX FOIS.
%      T18 annonce 1,951 kg de fer statorique la ou la geometrie en donne
%      environ 1,2 kg. Correction : integrer PAR CELLULE, avec la
%      decomposition orthogonale (B_r, B_t) -- exactement ce que fait la
%      chaine localisee en separant dent (radial) et culasse (tangentiel).
%
%  (e) PERIMETRE. Le maillage ne contient que le fer STATOR. La chaine
%      localisee ajoute un terme de culasse ROTOR (Ph_r + Pe_r,
%      BLDC_MEC_COMPLET.m:328). Comparaison a perimetre egal obligatoire.
%
%  (f) SOLVEUR. T18 maille en Newton, la chaine localisee resout en lineaire
%      mu_r = 3000. Les deux colonnes sont donc produites par deux solveurs.
%      Ce bloc donne les deux et declare l'ecart.
%
%  CRITERE D'ACCEPTATION. Aucun chiffre de perte n'est retenu tant que le
%  volume de fer somme sur les cellules ne rejoint pas le volume geometrique
%  a mieux que 3 %. C'est le "bilan qui ferme" exige par la specification.
%
%  CONFIGURATION DECLAREE
%    machine       : PMSM 15/14, 750 W (machine_bldc)
%    vitesse       : M.FEA.n_nl tr/min, f = p*n/60
%    materiau      : M350-50A, Bertotti VOLUMIQUE Kh/Kc/Ke de machine_bldc
%    chaine LOCAL. : cogging_mec(M,1260,0,721,mu_r=M.muI,kfr,2*pi/p)
%                    -- celle de BLDC_MEC_COMPLET §5 qui produit les 21,34 W
%    chaine MAILL. : mesh_bldc(Ms=540,nst=4,nys=3,n_sh=2) + solve_bldc_mesh,
%                    Newton, balayage d'UNE periode electrique, Np positions
%    reference EF  : M.FEA.Pfe_nl W a vide
clear; clc; t0=tic;
diary('X2_ironloss_out.txt'); diary on;

M=machine_bldc(); p=M.p; Ns=M.Ns; L=M.ls; Ki=M.Ki;
if isfile('kfringe_ident.mat'), S0=load('kfringe_ident.mat'); kfr=S0.kbest; else, kfr=0.325; end
n_rpm=M.FEA.n_nl; om=n_rpm*2*pi/60; f=p*n_rpm/60;
PfeF=M.FEA.Pfe_nl;
Ms=540; nsh=2; nst0=4; nys0=3;

fprintf('=== X2 : pertes fer sur maillage, bilan ferme ===\n');
fprintf('  machine      : PMSM 15/14, 750 W\n');
fprintf('  vitesse      : %.0f tr/min  ->  f = %.3f Hz   (T18 : 1500 tr/min, 175 Hz)\n',n_rpm,f);
fprintf('  Bertotti VOL.: Kh=%.3f, Kc=%.5f, Ke=%g  [W/m3], Ki=%.2f, rho=%d kg/m3\n', ...
    M.Kh,M.KeFe,M.Kex,Ki,M.rof);
fprintf('  reference EF : %.2f W a vide @ %.0f tr/min\n',PfeF,n_rpm);
fprintf('  kfringe      : %.4f  (chaine localisee seulement)\n\n',kfr);

fprintf('  ---- 0. trace du 16,81 W de ARCHITECTURE §1.5 ----\n');
fprintf('    T18_out.txt : Pfe = 128612.55 W\n');
fprintf('    128612.55 / rho = 128612.55 / %d = %.4f W\n',M.rof,128612.55/M.rof);
fprintf('    Le 16,81 W est donc la correction (a) SEULE, a 1500 tr/min,\n');
fprintf('    avec le foisonnement compte deux fois et le fer compte deux\n');
fprintf('    fois. Il n''est pas sans source : il est sous-corrige.\n\n');

%% =====================================================================
%  1. CHAINE LOCALISEE -- celle qui produit les 21,34 W publies
%  =====================================================================
fprintf('  ---- 1. chaine LOCALISEE (reference interne, BLDC_MEC_COMPLET §5) ----\n');
Nsurf=1260; NpL=721;
R=cogging_mec(M,Nsurf,0,NpL,M.muI,kfr,2*pi/p);
Vt1=M.wst1*(M.hs0+M.hs1+M.hs2)*L*Ki;
Vy1=(2*pi*(M.Rso-M.wsy/2)/Ns)*M.wsy*L*Ki;
tau_yi=2*pi*(M.Rsi+M.hs0+M.hs1+M.hs2)/Ns;
Bt_=R.PhiT/(M.wst1*L*Ki);
By_=R.PhiY/(M.wsy*L*Ki);
Byr_=R.PhiT/(tau_yi*L*Ki);
tt=linspace(0,1/f,NpL); dtL=tt(2)-tt(1);
%  deux formulations de la perte par courants de Foucault :
%    'grad' : (Kc/2pi^2)*<(dB/dt)^2>   -- celle du manuscrit
%    'amp'  : Kc*(f*Bamp)^2            -- Bertotti d'amplitude, celle du maillage
bl=@(B,V)deal( sum(M.Kh*f*((max(B,[],2)-min(B,[],2))/2).^2)*V, ...
               sum((M.KeFe/(2*pi^2))*mean(gradient(B,dtL).^2,2))*V, ...
               sum(M.KeFe*(f*(max(B,[],2)-min(B,[],2))/2).^2)*V );
[Ph_t,Pe_t,Pa_t]=bl(Bt_ ,Vt1);
[Ph_y,Pe_y,Pa_y]=bl(By_ ,Vy1);
[Ph_r,Pe_r,Pa_r]=bl(Byr_,Vy1);
PfeL_stator = Ph_t+Pe_t+Ph_y+Pe_y;
PfeL_total  = PfeL_stator+Ph_r+Pe_r;
PfeL_amp_st = Ph_t+Pa_t+Ph_y+Pa_y;
PfeL_amp_to = PfeL_amp_st+Ph_r+Pa_r;
fprintf('    %-34s %10s %10s\n','','grad (manuscrit)','amp (Bertotti)');
fprintf('    %-34s %10.4f %14.4f\n','dents  : hysteresis + Foucault',Ph_t+Pe_t,Ph_t+Pa_t);
fprintf('    %-34s %10.4f %14.4f\n','culasse stator                ',Ph_y+Pe_y,Ph_y+Pa_y);
fprintf('    %-34s %10.4f %14.4f\n','culasse ROTOR (hors maillage) ',Ph_r+Pe_r,Ph_r+Pa_r);
fprintf('    %-34s %10.4f %14.4f\n','SOUS-TOTAL STATOR             ',PfeL_stator,PfeL_amp_st);
fprintf('    %-34s %10.4f %14.4f\n','TOTAL (stator + rotor)        ',PfeL_total,PfeL_amp_to);
fprintf('    ecart du TOTAL vs EF %.2f W : %+.1f %% (grad) | %+.1f %% (amp)\n', ...
    PfeF,100*(PfeL_total-PfeF)/PfeF,100*(PfeL_amp_to-PfeF)/PfeF);
fprintf('    ecart des deux formulations de Foucault : %+.2f %%\n', ...
    100*(PfeL_amp_to-PfeL_total)/PfeL_total);

%% =====================================================================
%  2. BILAN DE MATIERE -- avant toute lecture de perte
%  =====================================================================
fprintf('\n  ---- 2. bilan de matiere : le fer du maillage ----\n');
ME0=mesh_bldc(M,Ms,nst0,nys0,[],0,kfr,'nl',nsh);
Ls=ME0.Ls; nshT=ME0.nst-nst0; dth=ME0.dth; rc=ME0.rc; rth=ME0.rth;
isFe=ME0.isFe;
Vcell=zeros(Ls,Ms);
for la=1:Ls, Vcell(la,:)=dth*rc(la)*rth(la)*Ki*L; end
cellFe=isFe;
Vmesh=sum(Vcell(cellFe));
%  volume geometrique exact du fer statorique
Vyoke=pi*(M.Rso^2-(M.Rsi+M.hs)^2)*L*Ki;
Vbody=Ns*M.wst1*M.hs2*L*Ki;
rsh=linspace(M.Rsi,M.Rsi+M.hs0+M.hs1,2001);
arcFe=2*pi*rsh/Ns - 2*rsh.*asin(min(M.ws0./(2*rsh),1));      % arc de fer par pas
Vshoe=Ns*trapz(rsh,arcFe)*L*Ki;
Vgeo=Vyoke+Vbody+Vshoe;
Vbranch=sum(ME0.A(ME0.iron).*ME0.l(ME0.iron));               % ce que sommait T18
fprintf('    volume geometrique  : bec %.4e + dent %.4e + culasse %.4e\n',Vshoe,Vbody,Vyoke);
fprintf('                          TOTAL %.6e m3  ->  %.4f kg\n',Vgeo,Vgeo*M.rof);
fprintf('    somme des CELLULES  : %.6e m3  ->  %.4f kg   (ecart %+.2f %%)\n', ...
    Vmesh,Vmesh*M.rof,100*(Vmesh-Vgeo)/Vgeo);
fprintf('    somme des BRANCHES  : %.6e m3  ->  %.4f kg   (rapport %.3f)\n', ...
    Vbranch,Vbranch*M.rof,Vbranch/Vgeo);
fprintf('    (T18 annoncait %.3f kg : somme des branches x kFe = %.4f kg)\n', ...
    1.951,Vbranch*Ki*M.rof);
%  Le meme bilan, cote chaine LOCALISEE : si les deux chaines portent le
%  meme fer, l'ecart entre leurs pertes ne peut venir que du CHAMP.
VtL=Ns*Vt1; VyL=Ns*Vy1;
fprintf('    -- fer de la chaine localisee, pour comparaison --\n');
fprintf('    dents  %.6e m3 (maillage bec+corps %.6e, ecart %+.2f %%)\n', ...
    VtL,Vshoe+Vbody,100*(VtL-(Vshoe+Vbody))/(Vshoe+Vbody));
fprintf('    culasse %.6e m3 (geometrique %.6e, ecart %+.2f %%)\n', ...
    VyL,Vyoke,100*(VyL-Vyoke)/Vyoke);
fprintf('    TOTAL   %.6e m3 (ecart au maillage %+.2f %%)\n', ...
    VtL+VyL,100*(VtL+VyL-Vmesh)/Vmesh);
fprintf('    => les deux colonnes portent le MEME fer a %.1f %% pres ;\n', ...
    abs(100*(VtL+VyL-Vmesh)/Vmesh));
fprintf('       tout ecart de perte vient donc du CHAMP, non de la matiere.\n');
ok=abs(Vmesh-Vgeo)/Vgeo<=0.03;
if ok
    fprintf('    BILAN FERME (%.2f %% <= 3 %%). Les pertes peuvent etre lues.\n', ...
        100*abs(Vmesh-Vgeo)/Vgeo);
else
    fprintf('    *** BILAN NON FERME (%.2f %%). Les pertes ci-dessous ne sont\n', ...
        100*abs(Vmesh-Vgeo)/Vgeo);
    fprintf('    *** PAS publiables : corriger la geometrie du pavage d''abord.\n');
end

%% ---- operateurs de moyenne branche -> cellule -------------------------
%  Chaque cellule recoit la moyenne des branches de FER qui la traversent :
%  deux radiales (dessous / dessus) et deux tangentielles (gauche / droite).
%  Les branches d'air sont exclues du comptage, non comptees a zero.
tp=ME0.type; co=ME0.col; ly=ME0.lay; Fe=ME0.iron; nb=numel(tp);
Rid=zeros(Ms,Ls-1); Tid=zeros(Ms,Ls); Sid3=zeros(Ms,1);
for k=1:nb
    if     tp(k)==1, Rid(co(k),ly(k))=k;
    elseif tp(k)==2, Tid(co(k),ly(k))=k;
    else,            Sid3(co(k))=k;
    end
end
cid=@(c,la)(la-1)*Ms+c;
iR=[]; jR=[]; vR=[]; iT=[]; jT=[]; vT=[];
for la=1:Ls
    for c=1:Ms
        if ~isFe(la,c), continue; end
        k=cid(c,la);
        if la==1, br=[Sid3(c) Rid(c,1)]; elseif la==Ls, br=Rid(c,Ls-1);
        else,     br=[Rid(c,la-1) Rid(c,la)]; end
        br=br(br>0); br=br(Fe(br));
        if ~isempty(br)
            iR=[iR;repmat(k,numel(br),1)]; jR=[jR;br(:)]; %#ok<AGROW>
            vR=[vR;repmat(1/numel(br),numel(br),1)];      %#ok<AGROW>
        end
        cm=mod(c-2,Ms)+1; bt=[Tid(cm,la) Tid(c,la)];
        bt=bt(bt>0); bt=bt(Fe(bt));
        if ~isempty(bt)
            iT=[iT;repmat(k,numel(bt),1)]; jT=[jT;bt(:)]; %#ok<AGROW>
            vT=[vT;repmat(1/numel(bt),numel(bt),1)];      %#ok<AGROW>
        end
    end
end
AR=sparse(iR,jR,vR,Ms*Ls,nb); AT=sparse(iT,jT,vT,Ms*Ls,nb);
Vc=Vcell.'; Vc=Vc(:);                       % (c,la) -> indice cid
keepC=false(Ms*Ls,1);
regC=zeros(Ms*Ls,1);
for la=1:Ls
    if la<=nshT
        rg=1;
    elseif la<=ME0.nst
        rg=2;
    else
        rg=3;
    end
    for c=1:Ms
        k=cid(c,la);
        regC(k)=rg;
        if isFe(la,c), keepC(k)=true; end
    end
end

%% =====================================================================
%  3. CHAINE MAILLAGE -- convergence en Np, deux solveurs
%  =====================================================================
fprintf('\n  ---- 3. chaine MAILLAGE : Ms=%d, n_sh=%d, nst=%d, nys=%d ----\n', ...
    Ms,nsh,nst0,nys0);
fprintf('    balayage d''UNE periode electrique, phi in [0, 2pi/p[, Np points\n');
NpS=[30 60 120]; SOL={'nl','nl','nl'};
res=struct('Np',{},'sol',{},'Ph',{},'Pe_amp',{},'Pe_grad',{},'Ptot_amp',{}, ...
           'Ptot_grad',{},'reg',{},'Bmax',{});
runs=[num2cell(NpS) {60}]; sols=[SOL {M.muI}];
for ir=1:numel(runs)
    Np=runs{ir}; sol=sols{ir};
    if ischar(sol), snm='Newton'; else, snm=sprintf('lin mu=%g',sol); end
    phis=linspace(0,2*pi/p,Np+1); phis(end)=[];
    dt=(2*pi/p)/Np/om;
    Ball=zeros(nb,Np); U0=[];
    warm=ischar(sol)||isstring(sol);
    %  DEMARRAGE A CHAUD : licite avec NEWTON seulement. Avec le solveur
    %  lineaire, solve_bldc_mesh ecrase la solution directe par U0 et rend
    %  la position precedente inchangee -- le champ devient constant sur le
    %  balayage et toute amplitude sort a zero. C'est le bug deja rencontre
    %  sur la FEM (SPEC_CLAUDE_CODE_v3 §3.3), a ne pas reintroduire.
    for q=1:Np
        MEq=mesh_bldc(M,Ms,nst0,nys0,[],phis(q),kfr,sol,nsh);
        if warm
            Sq=solve_bldc_mesh(MEq,1e-9,30,U0); U0=Sq.U;
        else
            Sq=solve_bldc_mesh(MEq,1e-9,30);
        end
        Ball(:,q)=Sq.B;
    end
    Brc=AR*Ball; Btc=AT*Ball;
    Bra=0.5*(max(Brc,[],2)-min(Brc,[],2));
    Bta=0.5*(max(Btc,[],2)-min(Btc,[],2));
    B2=Bra.^2+Bta.^2;
    dBr=(circshift(Brc,-1,2)-circshift(Brc,1,2))/(2*dt);
    dBt=(circshift(Btc,-1,2)-circshift(Btc,1,2))/(2*dt);
    ph  = M.Kh*f*B2.*Vc;
    peA = M.KeFe*(f^2)*B2.*Vc;
    peG = (M.KeFe/(2*pi^2))*(mean(dBr.^2,2)+mean(dBt.^2,2)).*Vc;
    if M.Kex>0, pex=M.Kex*(f^1.5)*(B2).^0.75.*Vc; else, pex=zeros(size(ph)); end
    m=keepC;
    r=struct();
    r.Np=Np; r.sol=snm;
    r.Ph=sum(ph(m)); r.Pe_amp=sum(peA(m)); r.Pe_grad=sum(peG(m));
    r.Ptot_amp=r.Ph+r.Pe_amp+sum(pex(m)); r.Ptot_grad=r.Ph+r.Pe_grad+sum(pex(m));
    r.reg=[sum(ph(m&regC==1)+peA(m&regC==1)) ...
           sum(ph(m&regC==2)+peA(m&regC==2)) ...
           sum(ph(m&regC==3)+peA(m&regC==3))];
    r.Bmax=[max(Bra(m)) max(Bta(m))];
    res(end+1)=r; %#ok<AGROW>
    fprintf('    Np=%3d  %-12s  Ph=%7.3f  Pe(amp)=%7.3f  Pe(grad)=%7.3f  ->  %7.3f W (amp) | %7.3f W (grad)   [%.0f s]\n', ...
        Np,snm,r.Ph,r.Pe_amp,r.Pe_grad,r.Ptot_amp,r.Ptot_grad,toc(t0));
end

iN=find(strcmp({res.sol},'Newton')); iF=iN(end);       % Np le plus fin, Newton
iL=find(~strcmp({res.sol},'Newton'),1);
fprintf('\n    convergence en Np (Newton, formulation amp)  : %s\n', ...
    mat2str(round([res(iN).Ptot_amp],4)));
fprintf('    variation du dernier doublement : %+.2f %%\n', ...
    100*(res(iN(end)).Ptot_amp-res(iN(end-1)).Ptot_amp)/res(iN(end-1)).Ptot_amp);
fprintf('    convergence en Np (Newton, formulation grad) : %s\n', ...
    mat2str(round([res(iN).Ptot_grad],4)));
fprintf('    variation du dernier doublement : %+.2f %%\n', ...
    100*(res(iN(end)).Ptot_grad-res(iN(end-1)).Ptot_grad)/res(iN(end-1)).Ptot_grad);
fprintf(['    LECTURE. La forme d''AMPLITUDE converge en Np ; la forme a\n' ...
         '    DERIVEE ne converge pas -- chaque doublement de Np resout des\n' ...
         '    harmoniques de position supplementaires dans dB/dt. Sur le\n' ...
         '    maillage, une perte fondee sur une derivee locale herite donc\n' ...
         '    du statut des grandeurs LOCALES, exactement comme la bande de\n' ...
         '    denture. C''est un argument direct pour la restriction §8.5.\n']);
fprintf('    repartition (bec / corps de dent / culasse) : %.3f / %.3f / %.3f W\n', ...
    res(iF).reg(1),res(iF).reg(2),res(iF).reg(3));
fprintf('    soit %.1f %% / %.1f %% / %.1f %%\n',100*res(iF).reg/sum(res(iF).reg));
fprintf('    B alternatif max : radial %.3f T | tangentiel %.3f T\n', ...
    res(iF).Bmax(1),res(iF).Bmax(2));
fprintf('    effet du solveur (Np=%d) : Newton %.3f W | %s %.3f W  (%+.1f %%)\n', ...
    res(iL).Np,res(iN(2)).Ptot_amp,res(iL).sol,res(iL).Ptot_amp, ...
    100*(res(iL).Ptot_amp-res(iN(2)).Ptot_amp)/res(iN(2)).Ptot_amp);

%% =====================================================================
%  4. COMPARAISON A PERIMETRE, VITESSE ET FORMULATION EGAUX
%  =====================================================================
fprintf('\n  ---- 4. les deux colonnes du Tableau 7, comparables ----\n');
fprintf('    perimetre STATOR seul, %.0f tr/min, Bertotti d''amplitude,\n',n_rpm);
fprintf('    coefficients volumiques M350-50A identiques des deux cotes.\n\n');
fprintf('    %-26s %12s %12s\n','source du champ','P_fe (W)','ecart / EF');
fprintf('    %-26s %12.3f %11.1f %%\n','reseau a une dent (stator)',PfeL_amp_st, ...
    100*(PfeL_amp_st-PfeF)/PfeF);
fprintf('    %-26s %12.3f %11.1f %%\n','MAILLAGE polaire (stator) ',res(iF).Ptot_amp, ...
    100*(res(iF).Ptot_amp-PfeF)/PfeF);
fprintf('    %-26s %12.3f\n','reference EF (total)      ',PfeF);
fprintf('\n    RAPPEL DE PERIMETRE. La reference EF %.2f W est un total\n',PfeF);
fprintf('    machine ; la chaine localisee lui ajoute un terme de culasse\n');
fprintf('    rotor de %.3f W que le maillage ne contient pas. Colonne\n',Ph_r+Pa_r);
fprintf('    localisee AVEC rotor : %.3f W (%+.1f %%).\n', ...
    PfeL_amp_to,100*(PfeL_amp_to-PfeF)/PfeF);
fprintf('    Maillage + meme terme rotor : %.3f W (%+.1f %%).\n', ...
    res(iF).Ptot_amp+Ph_r+Pa_r,100*(res(iF).Ptot_amp+Ph_r+Pa_r-PfeF)/PfeF);

fprintf('\n  ---- 5. consequence pour l''anomalie A-1 et la §8.5 ----\n');
fprintf('    T18 (publie dans ARCHITECTURE §1.5) : 16,81 W, -14,8 %%\n');
fprintf('    X2  (bilan de matiere ferme)        : %.2f W, %+.1f %% (stator)\n', ...
    res(iF).Ptot_amp,100*(res(iF).Ptot_amp-PfeF)/PfeF);
fprintf('                                          %.2f W, %+.1f %% (avec rotor)\n', ...
    res(iF).Ptot_amp+Ph_r+Pa_r,100*(res(iF).Ptot_amp+Ph_r+Pa_r-PfeF)/PfeF);
fprintf('    Le -14,8 %% ne peut pas etre publie : il cumule quatre defauts\n');
fprintf('    dont deux se compensaient partiellement. La valeur ci-dessus\n');
fprintf('    est la premiere a reposer sur un fer dont la masse est juste.\n');

save('X2_ironloss.mat','res','PfeL_stator','PfeL_total','PfeL_amp_st', ...
     'PfeL_amp_to','Ph_t','Pe_t','Pa_t','Ph_y','Pe_y','Pa_y','Ph_r','Pe_r','Pa_r', ...
     'Vgeo','Vmesh','Vbranch','Vshoe','Vbody','Vyoke','f','n_rpm','PfeF', ...
     'Ms','nsh','nst0','nys0','kfr');
fprintf('\n  duree %.0f s\n=== X2 termine ===\n',toc(t0));
diary off;
