%% RUN_R6A_GRID - R6, volet MEC : l'inductance incrementale et son pas de grille
%
%  CE QUE CE BLOC TESTE. Le transitoire de drive_mec n'est PAS integre sur la
%  surface lambda(theta,i) : il l'est sur les DERIVEES La, Lb lues dans la
%  carte (drive_mec.m:147-148, Leq = sum(d.*(La*da+Lb*db)) puis di/dt = .../Leq).
%  Or mec_map.m:196-201 et 227-230 forment ces derivees par la fonction
%  gradient() sur la grille de courant iax. Le PAS de cette grille fixe donc
%  directement la constante de temps electrique simulee.
%
%  En production (BLDC_MEC_COMPLET.m:388-396) cette grille compte 19 points sur
%  [-28,28] en loi quadratique : le pas local vaut plusieurs amperes entre 5 et
%  20 A, precisement la ou l'inductance s'effondre de 104 a 12 mH. On balaye
%  donc le NOMBRE DE POINTS a loi et a chaine inchangees, et on regarde si L
%  converge. C'est le protocole de R4 transpose.
%
%  CACHE PERIME, ET C'EST UN RESULTAT. mec_map.mat date du 28/07 alors que
%  machine_bldc.m et airgap_magnet.m ont ete modifies le 04/08 : le cache de
%  BLDC_MEC_COMPLET n'est indexe que sur le NOM de fichier, sans invalidation
%  sur la source. On ne lit donc aucun cache ici : tout est reconstruit.
%
%  GRANDEUR ET CHAINE. L'inductance de ligne est formee par la MEME expression
%  que celle du programme maitre (local_Lline, recopiee VERBATIM en fin de
%  fichier), evaluee au meme motif de conduction : phase a en serie avec phase
%  b, soit (i_alpha, i_beta) = (i, -i/sqrt(3)).
%
%  GARDE. A la grille de production, L(0) et L(25) doivent reproduire les deux
%  valeurs publiees en Fig. 12, RELUES du manuscrit et non transcrites.
clear; clc; t0=tic;
diary('R6a_grid_out.txt'); diary on;
M=machine_bldc();
KFR_A=0.75;              % role induit/entrainement (BLDC_MEC_COMPLET.m:62)
NSURF=720; NTH=61;       % valeurs de production (BLDC_MEC_COMPLET.m:396)
NA=[19 27 35 43];        % grille de courant balayee ; 19 = production
IQ=[0 2 5 10 15 20 25];  % courants d'evaluation (R6 en demande cinq)
TEX=fullfile('..','article','ArticleI_DtN_PMSM.tex');

fprintf('=== R6a : inductance incrementale et pas de grille en courant ===\n\n');
fprintf('  CONFIGURATION\n');
fprintf('    machine      : PMSM 15/14, 750 W (machine_bldc)\n');
fprintf('    chaine       : mec_map + local_Lline (recopiee verbatim)\n');
fprintf('    Nsurf        : %d noeuds de surface\n',NSURF);
fprintf('    Nth          : %d positions rotor par point de courant\n',NTH);
fprintf('    k_fringe     : %.4f (role induit)\n',KFR_A);
fprintf('    solveur      : Newton, courbe B(H) M350-50A (bh_curve)\n');
fprintf('    loi de grille: 28*sign(u)*|u|^2, u = linspace(-1,1,na)\n');
fprintf('    na balaye    : %s   (production = %d)\n',mat2str(NA),NA(1));
fprintf('    conduction   : phase a en serie avec phase b, (ia,ib)=(i,-i/sqrt(3))\n');
fprintf('    cache        : AUCUN, tout reconstruit\n');

%% ---- valeurs publiees, RELUES du manuscrit ---------------------------
if ~isfile(TEX), error('R6a:tex','manuscrit introuvable : %s',TEX); end
tx=fileread(TEX);
tk=regexp(tx,['collapsing from\s*\\SI\{([\d.]+)\}\{\\milli\\henry\}\s*to\s*' ...
              '\\SI\{([\d.]+)\}\{\\milli\\henry\}'],'tokens','once');
if isempty(tk), error('R6a:parse','valeurs de la Fig. 12 non reconnues dans le .tex'); end
Lpub=[str2double(tk{1}) str2double(tk{2})];
fprintf('\n  valeurs publiees relues du manuscrit (Fig. 12) :\n');
fprintf('    L_ligne(0 A) = %.1f mH   |   L_ligne(25 A) = %.1f mH\n',Lpub);

%% ---- balayage du pas de grille ---------------------------------------
K=numel(NA); nq=numel(IQ);
Lg=nan(K,nq); Lsec=nan(K,nq); lam_arm=nan(K,nq); pas=nan(K,nq); tps=nan(1,K);
for k=1:K
    na=NA(k); uu=linspace(-1,1,na); iax=28*sign(uu).*abs(uu).^2;
    tk0=tic;
    S=mec_map(M,NSURF,KFR_A,NTH,iax);
    tps(k)=toc(tk0);
    for q=1:nq
        ii=IQ(q);
        Lg(k,q)=local_Lline(S,ii);
        [lam_arm(k,q),Lsec(k,q)]=ray_secant(S,ii);
        pas(k,q)=local_step(iax,ii);
    end
    fprintf('\n  ---- na = %d (%.0f s) ----\n',na,tps(k));
    fprintf('  %8s %10s %14s %14s %14s\n', ...
        'i (A)','pas (A)','L increm (mH)','L secante (mH)','lam_arm (Wb)');
    for q=1:nq
        fprintf('  %8.1f %10.3f %14.4f %14.4f %14.5f\n', ...
            IQ(q),pas(k,q),Lg(k,q)*1e3,Lsec(k,q)*1e3,lam_arm(k,q));
    end
end

%% ---- GARDE : la grille de production reproduit-elle la Fig. 12 ? ------
fprintf('\n  ---- GARDE : grille de production contre Fig. 12 ----\n');
i0=find(IQ==0,1); i25=find(IQ==25,1);
Lg0=Lg(1,i0)*1e3; Lg25=Lg(1,i25)*1e3;
fprintf('    L(0 A)  : chaine %.2f mH | publie %.1f mH | ecart %+.2f %%\n', ...
    Lg0,Lpub(1),100*(Lg0-Lpub(1))/Lpub(1));
fprintf('    L(25 A) : chaine %.2f mH | publie %.1f mH | ecart %+.2f %%\n', ...
    Lg25,Lpub(2),100*(Lg25-Lpub(2))/Lpub(2));
tolg=0.5;                                  % demi-unite du dernier chiffre publie
G1=abs(Lg0-Lpub(1))<=tolg && abs(Lg25-Lpub(2))<=tolg;
if G1
    fprintf('    GARDE PASSEE : la chaine reconstruite reproduit la figure publiee.\n');
else
    fprintf(['    GARDE ECHOUEE : la chaine reconstruite NE reproduit PAS la figure.\n' ...
             '    Cause a etablir avant toute conclusion -- le cache mec_map.mat est\n' ...
             '    anterieur a la derniere modification de machine_bldc.m.\n']);
end

%% ---- convergence en pas de grille -------------------------------------
fprintf('\n  ---- CONVERGENCE DE L''INDUCTANCE INCREMENTALE ----\n');
fprintf('  %8s',' i (A)'); fprintf(' %12s',compose("na=%d",NA)); fprintf(' %12s %12s\n','deriv./prod.','dispersion');
for q=1:nq
    fprintf('  %8.1f',IQ(q));
    fprintf(' %12.4f',Lg(:,q)*1e3);
    d=100*(Lg(end,q)-Lg(1,q))/Lg(1,q);
    dp=100*(max(Lg(:,q))-min(Lg(:,q)))/mean(Lg(:,q));
    fprintf(' %11.2f %% %11.2f %%\n',d,dp);
end

fprintf('\n  ---- SECANTE CONTRE INCREMENTALE (grille la plus fine) ----\n');
fprintf('  Les confondre fabriquerait exactement l''ecart cherche : on les affiche.\n');
fprintf('  %8s %16s %16s %12s\n','i (A)','increm. (mH)','secante (mH)','rapport');
for q=1:nq
    fprintf('  %8.1f %16.4f %16.4f %12.4f\n', ...
        IQ(q),Lg(end,q)*1e3,Lsec(end,q)*1e3,Lsec(end,q)/Lg(end,q));
end

save('R6a_grid.mat','NA','IQ','Lg','Lsec','lam_arm','pas','tps','Lpub','G1', ...
     'NSURF','NTH','KFR_A');
fprintf('\n  duree totale %.0f s\n=== R6a termine ===\n',toc(t0));
diary off;

% ======================================================================
function L = local_Lline(S,ii)
%LOCAL_LLINE  RECOPIE VERBATIM de BLDC_MEC_COMPLET.m:1277-1283.
%   Inductance de ligne incrementale lue dans la cartographie, pour le motif
%   de conduction reel (phase a en serie avec phase b).
    ia=ii; ib=-ii/sqrt(3); s3=sqrt(3);
    q=@(X,k) interpn(1:3,S.the,S.iax,S.iax,X,k,S.the,ia,ib,'linear');
    L=mean((q(S.La,1)-q(S.La,2))-(q(S.Lb,1)-q(S.Lb,2))/s3);
end

function [lam_arm,Lsec] = ray_secant(S,ii)
%RAY_SECANT  Flux d'INDUIT de ligne et inductance SECANTE au meme point.
%   lam_arm(i) = <lambda_a - lambda_b>_theta a (i,-i/sqrt3) MOINS la meme
%   moyenne a courant nul : on retranche ainsi la contribution des aimants,
%   qui n'est pas une inductance. Lsec = lam_arm/i.
    ia=ii; ib=-ii/sqrt(3);
    q=@(k,a,b) interpn(1:3,S.the,S.iax,S.iax,S.lam,k,S.the,a,b,'linear');
    d1=mean(q(1,ia,ib)-q(2,ia,ib));
    d0=mean(q(1,0,0)  -q(2,0,0));
    lam_arm=d1-d0;
    if abs(ii)<eps, Lsec=NaN; else, Lsec=lam_arm/ii; end
end

function h = local_step(iax,ii)
%LOCAL_STEP  Demi-largeur effective de la difference centree de gradient()
%   au voisinage du courant ii : distance entre les deux noeuds encadrants.
    v=sort(iax(:).');
    j=find(v<=ii,1,'last'); k=find(v>=ii,1,'first');
    if isempty(j), j=1; end
    if isempty(k), k=numel(v); end
    if j==k
        jl=max(1,j-1); jr=min(numel(v),j+1); h=(v(jr)-v(jl))/2;
    else
        h=v(k)-v(j);
    end
end
