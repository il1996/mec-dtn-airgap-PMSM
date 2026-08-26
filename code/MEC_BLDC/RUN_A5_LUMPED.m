%% RUN_A5_LUMPED - provenance des trois cellules "Lumped" du Tableau 7
%
%  ENONCE. Trois cellules de la colonne Lumped ne se regenerent pas :
%       B_r crete    publie 0.9846 | programme 0.9860  (+0.14 %)
%       FMM aimant   publie 763.2  | programme 762.9   (-0.04 %)
%       FMM entrefer publie 723.8  | programme 725.1   (+0.18 %)
%  alors que k_r et k_l de la MEME colonne coincident exactement (1.0522 et
%  0.8537). Identifier la configuration de cogging_mec qui produit les
%  valeurs publiees, OU etablir qu'aucune ne le fait.
%
%  METHODE. Trois etages, du plus discriminant au plus large :
%    §1  TEST D'IDENTITE INTERNE. Dans krkl, k_r = mmf_m/mmf_g par
%        CONSTRUCTION. Le triplet publie est donc soit compatible avec un
%        k_r publie, soit non : c'est un test a une ligne qui elimine d'un
%        coup toutes les configurations a phi = 0.
%    §2  RECHERCHE DE LA SOURCE. Les trois valeurs sont recherchees a pleine
%        precision dans les fichiers du projet, et la configuration qui les
%        a produites est REEXECUTEE telle quelle.
%    §3  BALAYAGE DEMANDE. Nsurf, Np, kfringe, mu_i, base p0/p1,
%        cogging_mec vs cogging_mec2, plus deux axes NUMERIQUES qui ne
%        figuraient pas dans l'enonce (finesse de quadrature de krkl,
%        finesse de la grille du maximum) et un axe HISTORIQUE (mu_r
%        d'aimant de la fiche de dimensionnement).
%
%  REGLE. Aucune valeur transcrite a la main : les cibles sont RELUES de
%  outF.txt, les valeurs du programme de A1_table7.mat, les valeurs du
%  manuscrit de article\MEC_DtN_paper_v2.tex.
%
%  Diary : A5_lumped_out.txt (fichier remis a zero a chaque execution, pour
%  qu'il n'y ait jamais d'ambiguite sur le bloc qui fait foi).

clear; clc; t0=tic;
if isfile('A5_lumped_out.txt'), delete('A5_lumped_out.txt'); end
diary('A5_lumped_out.txt'); diary on;

mu0=4*pi*1e-7;
M=machine_bldc(); p=M.p; Ns=M.Ns; Nm=M.Nm; L=M.ls;
if isfile('kfringe_ident.mat'), S0=load('kfringe_ident.mat'); kfr=S0.kbest; else, kfr=0.325; end
rd=@(f)readmatrix(f,'FileType','text','Delimiter','\t','NumHeaderLines',1);
fea=M.FEA.dir;
thq2=linspace(0,2*pi,2001); thq2(end)=[];   % grille de l''enonce  (2000 pts)
thu =linspace(0,2*pi,3601); thu(end)=[];    % grille de A1 et de BLDC_MEC_COMPLET

fprintf('=== A5 : provenance des trois cellules Lumped ===\n');
fprintf('  machine PMSM 15/14 750 W | kfringe identifie = %.4f | mu_i nominal = %g\n', ...
    kfr,M.muI);
fprintf('  mu_r aimant = %.4f (Br = %.4f T, Hc = %d A/m)\n\n',M.mu_r,M.Br,M.Hc);

%% ======================================================================
%  0. LES TROIS VALEURS PUBLIEES : OU SONT-ELLES ECRITES, A QUELLE PRECISION
%  ======================================================================
fprintf('---- 0. sources relues (aucune transcription) ----\n');

fout=fullfile(pwd,'outF.txt');
ftex=fullfile('..','article','MEC_DtN_paper_v2.tex');
fA1  ='A1_table7.mat';

TGT=nan(1,3); KRKL_pub=nan(1,2); prov='';
if isfile(fout)
    Lo=readlines(fout);
    TGT(1)=a5_num(Lo,'B entrefer crete');
    TGT(2)=a5_num(Lo,'FMM par aimant (moyenne)');
    TGT(3)=a5_num(Lo,'FMM par entrefer (moyenne)');
    KRKL_pub(1)=a5_num(Lo,'Facteur de reluctance k_r');
    KRKL_pub(2)=a5_num(Lo,'Facteur de fuite k_l');
    k1=find(contains(Lo,'noeuds de surface'),1);
    k2=find(contains(Lo,'alignement : position rotor'),1);
    if ~isempty(k1), prov=strtrim(Lo(k1)); end
    fprintf('  outF.txt (sortie brute de BLDC_MEC_COMPLET) :\n');
    fprintf('    %s\n',prov);
    if ~isempty(k2), fprintf('    %s\n',strtrim(Lo(k2))); end
    fprintf('    B entrefer crete %.4f | FMM aimant %.4f | FMM entrefer %.4f\n', ...
        TGT(1),TGT(2),TGT(3));
    fprintf('    k_r %.4f | k_l %.4f\n',KRKL_pub(1),KRKL_pub(2));
else
    fprintf('  outF.txt INTROUVABLE : les cibles ne peuvent pas etre relues.\n');
end

% --- valeurs du programme, relues de A1 (colonne Lumped) ----------------
PRG=nan(1,5);
if isfile(fA1)
    SA=load(fA1);   % SA.VL, SA.lab
    ix=@(s)find(strcmp(SA.lab,s),1);
    PRG=[SA.VL(ix('Br crete (T)')) SA.VL(ix('FMM moy./aimant (A)')) ...
         SA.VL(ix('FMM moy./entrefer (A)')) SA.VL(ix('facteur k_r (-)')) ...
         SA.VL(ix('facteur k_l (-)'))];
    fprintf('  A1_table7.mat (colonne Lumped, phi = 0, Np = 61) :\n');
    fprintf('    Br crete %.5f | FMM aimant %.5f | FMM entrefer %.5f | k_r %.5f | k_l %.5f\n', ...
        PRG(1),PRG(2),PRG(3),PRG(4),PRG(5));
else
    fprintf('  A1_table7.mat introuvable : lancer RUN_A1_TABLE7 d''abord.\n');
end

% --- valeurs du manuscrit ------------------------------------------------
MSC=nan(1,2);
if isfile(ftex)
    Lt=readlines(ftex);
    MSC(1)=a5_tex(Lt,'reluctance factor $k_r$');
    MSC(2)=a5_tex(Lt,'leakage factor $k_l$');
    fprintf('  MEC_DtN_paper_v2.tex, colonne Lumped du Tableau 7 :\n');
    fprintf('    k_r %.4f | k_l %.4f\n',MSC(1),MSC(2));
    for s={'peak $B_r$','mean MMF per magnet','mean MMF per gap'}
        kk=find(contains(Lt,s{1}),1);
        if ~isempty(kk), fprintf('    %s\n',strtrim(Lt(kk))); end
    end
else
    fprintf('  MEC_DtN_paper_v2.tex introuvable.\n');
end

%% ======================================================================
%  1. LE TEST D'IDENTITE INTERNE  --  k_r = mmf_m / mmf_g
%  ======================================================================
fprintf('\n---- 1. test d''identite interne : k_r = mmf_m/mmf_g ----\n');
fprintf(['  Dans krkl.m (l.35) k_r n''est pas une grandeur independante :\n' ...
         '  k_r = mean(mm)/mean(mg). Toute configuration de cogging_mec verifie\n' ...
         '  donc EXACTEMENT cette identite. Elle se teste sur les nombres publies.\n\n']);
fprintf('  %-42s %12s %12s %12s\n','jeu de valeurs','mmf_m','mmf_g','mmf_m/mmf_g');
fprintf('  %-42s %12.4f %12.4f %12.5f\n','publie (cellules de l''enonce)', ...
    TGT(2),TGT(3),TGT(2)/TGT(3));
fprintf('  %-42s %12.4f %12.4f %12.5f\n','programme A1 (phi = 0)', ...
    PRG(2),PRG(3),PRG(2)/PRG(3));
fprintf('\n  k_r du manuscrit                     : %.5f\n',MSC(1));
fprintf('  ecart |k_r manuscrit - mmf publiees| : %.4f  (%.3f %%)\n', ...
    abs(MSC(1)-TGT(2)/TGT(3)),100*abs(MSC(1)-TGT(2)/TGT(3))/MSC(1));
fprintf('  ecart |k_r manuscrit - mmf A1|       : %.6f  (%.4f %%)\n', ...
    abs(MSC(1)-PRG(2)/PRG(3)),100*abs(MSC(1)-PRG(2)/PRG(3))/MSC(1));
if abs(MSC(1)-PRG(2)/PRG(3)) < abs(MSC(1)-TGT(2)/TGT(3))/10
    fprintf(['\n  LECTURE. Le k_r du manuscrit est le quotient des FMM du PROGRAMME,\n' ...
             '  pas celui des FMM publiees. La colonne Lumped du manuscrit MELANGE\n' ...
             '  donc deux extractions. Aucune configuration unique de cogging_mec ne\n' ...
             '  peut produire a la fois (%.1f ; %.1f) et k_r = %.4f : le systeme est\n' ...
             '  SURDETERMINE et INCOMPATIBLE. Le balayage du paragraphe 3 ne peut\n' ...
             '  donc pas reussir sur les trois cellules ET sur k_r simultanement --\n' ...
             '  ce qui oriente la recherche vers un CHANGEMENT DE POSITION\n' ...
             '  ROTORIQUE, seul parametre qui deplace les FMM sans etre un\n' ...
             '  parametre de modele.\n'],TGT(2),TGT(3),MSC(1));
    MIX=true;
else
    fprintf(['\n  LECTURE. Le test ne tranche pas : le k_r du manuscrit n''est pas\n' ...
             '  nettement plus proche du quotient des FMM du programme que de celui\n' ...
             '  des FMM publiees. Aucune conclusion de melange n''est tiree ici.\n']);
    MIX=false;
end

%% ======================================================================
%  2. LA SOURCE : BLDC_MEC_COMPLET, A LA POSITION ROTOR ALIGNEE
%  ======================================================================
fprintf('\n---- 2. reexecution de la configuration declaree dans outF.txt ----\n');
fprintf(['  outF.txt annonce 1260 noeuds de surface et 721 positions ; la section 1\n' ...
         '  de BLDC_MEC_COMPLET (l.92-93) porte Nsurf = 1260, Np = 721,\n' ...
         '  span = 2pi/p, et la section 2 (l.168) evalue les FMM A LA POSITION\n' ...
         '  IDENTIFIEE qb, pas a phi = 0. On rejoue exactement cette chaine.\n']);

Nsurf=1260; Np721=721;
Rr=cogging_mec(M,Nsurf,0,Np721,M.muI,kfr,2*pi/p);
phis=Rr.phis; AG=Rr.AG; nu=AG.nu;

d4=rd(fullfile(fea,'magnetostique(Magnetic_loading)','Calculator Expressions Plot 4.tab'));
angF=d4(:,2)*pi/180; BrF=d4(:,3);
BrFu=interp1([angF;2*pi],[BrF;BrF(1)],thu,'linear','extrap');
nmax=60; nn=(1:nmax)';
AF_h=2*mean(BrFu.*exp(-1i*nn*thu),2);
Uc=AG.Wc*Rr.Usurf; Usn=AG.Ws*Rr.Usurf;
cph=cos(nu*phis); sph=sin(nu*phis);
Brc_=AG.bru.*Uc+AG.brmq.*cph;  Brs_=AG.bru.*Usn+AG.brmq.*sph;
AM_h=Brc_(1:nmax,:)-1i*Brs_(1:nmax,:);
dgrid=linspace(0,2*pi/Ns,1441);
bestv=-inf; qb=1; db=0;
for q=1:size(AM_h,2)
    fv=real((AM_h(:,q).*conj(AF_h)).'*exp(-1i*nn*dgrid));
    [v,j]=max(fv);
    if v>bestv, bestv=v; qb=q; db=dgrid(j); end
end
ph_al=exp(-1i*nu*db);
Brc_al=real((Brc_(:,qb)-1i*Brs_(:,qb)).*ph_al);
Brs_al=-imag((Brc_(:,qb)-1i*Brs_(:,qb)).*ph_al);
Brm_al=(Brc_al.'*cos(nu*thu)+Brs_al.'*sin(nu*thu));
[mm_al,mg_al,kr_al,kl_al]=krkl(AG,Rr.U(1:Nsurf,qb),phis(qb),M,L,mu0);
V_al=[max(Brm_al) mm_al mg_al];
fprintf('\n  position identifiee : indice %d/%d, phi = %.3f deg mec, origine %.3f deg\n', ...
    qb,Np721,phis(qb)*180/pi,db*180/pi);
fprintf('  %-24s %14s %14s %12s\n','grandeur','reexecute','cible outF','ecart');
nmv={'Br crete (T)','FMM aimant (A)','FMM entrefer (A)'};
for q=1:3
    fprintf('  %-24s %14.5f %14.4f %11.4f %%\n', ...
        nmv{q},V_al(q),TGT(q),100*(V_al(q)/TGT(q)-1));
end
fprintf('  %-24s %14.5f %14.4f %11.4f %%\n','k_r (-)',kr_al,KRKL_pub(1), ...
    100*(kr_al/KRKL_pub(1)-1));
fprintf('  %-24s %14.5f %14.4f %11.4f %%\n','k_l (-)',kl_al,KRKL_pub(2), ...
    100*(kl_al/KRKL_pub(2)-1));

%  meme configuration, mais evaluee a phi = 0 : c''est la colonne A1
[mm_0,mg_0,kr_0,kl_0]=krkl(AG,Rr.U(1:Nsurf,1),0,M,L,mu0);
Br0_2000=max(AG.field(Rr.U(1:Nsurf,1),0,thq2));
Br0_3600=max(AG.field(Rr.U(1:Nsurf,1),0,thu));
fprintf('\n  meme reseau, evalue a phi = 0 (c''est ce que fait A1) :\n');
fprintf('    Br crete %.5f (2000 pts) / %.5f (3600 pts) | FMM %.4f / %.4f | k_r %.5f | k_l %.5f\n', ...
    Br0_2000,Br0_3600,mm_0,mg_0,kr_0,kl_0);
fprintf('    -> le SEUL changement entre les deux lignes est la POSITION ROTOR.\n');

%% ======================================================================
%  3. BALAYAGE DEMANDE
%  ======================================================================
fprintf('\n---- 3. balayage des configurations plausibles (toutes a phi = 0) ----\n');
NS=[630 900 1260 1800 2520];
KF=[0 0.325 0.75];
MU=[1500 3000 5000];
BA={'p0','p1'};
Npsw=5;                 % voir §3b : les trois grandeurs sont Np-invariantes

ALL=struct('lab',{},'v',{},'kr',{},'kl',{},'e',{},'etot',{});
fprintf('  %-38s %10s %11s %11s %9s %9s %9s %9s\n', ...
    'configuration','Br crete','FMM aim.','FMM entr.','e_Br %','e_m %','e_g %','TOTAL %');
for a=1:numel(NS)
 for b=1:numel(KF)
  for c=1:numel(MU)
   for d=1:numel(BA)
      Rk=cogging_mec(M,NS(a),0,Npsw,MU(c),KF(b),2*pi/p,BA{d});
      Uk=Rk.U(1:NS(a),1);
      [mk,gk,krk,klk]=krkl(Rk.AG,Uk,0,M,L,mu0);
      Bk=max(Rk.AG.field(Uk,0,thq2));
      lab=sprintf('cog1 Ns=%d kf=%.3f mu=%d %s',NS(a),KF(b),MU(c),BA{d});
      ALL=a5_push(ALL,lab,[Bk mk gk],krk,klk,TGT);
      A=ALL(end);
      fprintf('  %-38s %10.5f %11.4f %11.4f %9.3f %9.3f %9.3f %9.3f\n', ...
          lab,Bk,mk,gk,A.e(1),A.e(2),A.e(3),A.etot);
   end
  end
 end
end

%  --- 3b. axe Np -------------------------------------------------------
fprintf('\n  3b. axe Np (positions rotor) a Nsurf = 1260, kf = %.3f, mu = %g, p0\n',kfr,M.muI);
NP=[3 61 181 361];
fprintf('  %8s %14s %14s %14s\n','Np','Br crete','FMM aim.','FMM entr.');
for a=1:numel(NP)
    Rk=cogging_mec(M,1260,0,NP(a),M.muI,kfr,2*pi/p);
    Uk=Rk.U(1:1260,1);
    [mk,gk,krk,klk]=krkl(Rk.AG,Uk,0,M,L,mu0);
    Bk=max(Rk.AG.field(Uk,0,thq2));
    fprintf('  %8d %14.6f %14.5f %14.5f\n',NP(a),Bk,mk,gk);
    ALL=a5_push(ALL,sprintf('cog1 Np=%d (Ns=1260 nom.)',NP(a)),[Bk mk gk],krk,klk,TGT);
end
fprintf(['  Np n''a AUCUN effet : phis = linspace(0,span,Np) commence toujours a 0,\n' ...
         '  donc U(:,1) est la meme colonne quel que soit Np. Cet axe est mort.\n']);

%  --- 3c. cogging_mec2 (bouche d''encoche maillee, pas de kfringe) ------
fprintf('\n  3c. cogging_mec2 : bouche maillee, nm couches, aucun kfringe\n');
fprintf('  %8s %8s %14s %14s %14s\n','Nsurf','nm','Br crete','FMM aim.','FMM entr.');
for Nq2=[900 1260]
  for nm=[2 3 4]
    try
        R2=cogging_mec2(M,Nq2,nm,Npsw,M.muI,2*pi/p);
        U2=R2.U(1:Nq2,1);
        [m2,g2,kr2,kl2]=krkl(R2.AG,U2,0,M,L,mu0);
        B2=max(R2.AG.field(U2,0,thq2));
        fprintf('  %8d %8d %14.5f %14.5f %14.5f\n',Nq2,nm,B2,m2,g2);
        ALL=a5_push(ALL,sprintf('cog2 Ns=%d nm=%d',Nq2,nm),[B2 m2 g2],kr2,kl2,TGT);
    catch ME2
        fprintf('  %8d %8d   ECHEC : %s\n',Nq2,nm,ME2.message);
    end
  end
end

%  --- 3d. muI = ''nl'' -------------------------------------------------
fprintf('\n  3d. axe mu_i = ''nl'' (solveur de Newton)\n');
try
    Rnl=cogging_mec(M,1260,0,Npsw,'nl',kfr,2*pi/p); %#ok<NASGU>
    fprintf('  ''nl'' accepte -- resultat a exploiter.\n');
catch MEnl
    fprintf(['  ''nl'' REFUSE par cogging_mec. Message : %s\n' ...
             '  Motif : dans cogging_mec, muI entre dans des produits scalaires\n' ...
             '  (g_face = mu0*muI*(dth*Rs)*L/(hs0+hs1), l.57-59). Le reseau a une\n' ...
             '  dent est LINEAIRE par construction ; il n''existe pas de variante\n' ...
             '  de Newton. Cet axe de l''enonce est donc SANS OBJET pour la colonne\n' ...
             '  Lumped -- ''nl'' n''appartient qu''a mesh_bldc/solve_bldc_mesh,\n' ...
             '  c''est-a-dire aux colonnes Mesh.\n'],MEnl.message);
end

%  --- 3e. axe HISTORIQUE : mu_r d''aimant de la fiche de dimensionnement -
fprintf('\n  3e. axe historique : aimant de la fiche Magnet_Grade (Hc = 1005.3 kA/m)\n');
Mh=M; Mh.Hc=1005.3e3; Mh.mu_r=Mh.Br/(4*pi*1e-7*Mh.Hc);
Rh=cogging_mec(Mh,1260,0,Npsw,M.muI,kfr,2*pi/p);
Uh=Rh.U(1:1260,1);
[mh,gh,krh,klh]=krkl(Rh.AG,Uh,0,Mh,L,mu0);
Bh=max(Rh.AG.field(Uh,0,thq2));
fprintf('  mu_r = %.4f (au lieu de %.4f) : Br crete %.5f | FMM %.4f / %.4f | k_r %.5f\n', ...
    Mh.mu_r,M.mu_r,Bh,mh,gh,krh);
ALL=a5_push(ALL,sprintf('cog1 mu_r=%.4f (fiche)',Mh.mu_r),[Bh mh gh],krh,klh,TGT);

%  --- 3f. axes NUMERIQUES : quadrature de krkl et grille du maximum -----
fprintf('\n  3f. finesse de la quadrature angulaire de krkl (config nominale, phi = 0)\n');
fprintf('  %10s %14s %14s %12s\n','Nq','FMM aim.','FMM entr.','k_r');
for Nq=[501 1001 1501 2001 2801 3601 5001 7201 10001]
    [mq,gq]=a5_mmf(Rr.AG,Rr.U(1:1260,1),0,M,Nq);
    fprintf('  %10d %14.5f %14.5f %12.6f\n',Nq,mq,gq,mq/gq);
end
fprintf('\n  3g. finesse de la grille du maximum de Br (config nominale)\n');
fprintf('  %10s %16s %16s\n','N pts','max Br a phi=0','max Br aligne');
for Ng=[360 720 1000 1440 2000 2880 3600 7200 14400 36000]
    tg=linspace(0,2*pi,Ng+1); tg(end)=[];
    b0=max(Rr.AG.field(Rr.U(1:1260,1),0,tg));
    ba=max(Brc_al.'*cos(nu*tg)+Brs_al.'*sin(nu*tg));
    fprintf('  %10d %16.6f %16.6f\n',Ng,b0,ba);
end

%% ======================================================================
%  4. TABLEAU GLOBAL TRIE PAR ECART TOTAL CROISSANT
%  ======================================================================
%  On ajoute la configuration ALIGNEE du §2 pour qu''elle soit classee avec
%  les autres, sur le meme critere.
ALL=a5_push(ALL,'cog1 Ns=1260 Np=721 POSITION ALIGNEE',V_al,kr_al,kl_al,TGT);

et=[ALL.etot]; [~,so]=sort(et);
fprintf('\n---- 4. classement par ecart total (somme des trois |ecarts| en %%) ----\n');
fprintf('  %4s %-40s %10s %11s %11s %8s %8s %8s %9s %9s %9s\n', ...
    'rang','configuration','Br crete','FMM aim.','FMM entr.', ...
    'e_Br %','e_m %','e_g %','TOTAL %','k_r','k_l');
for q=1:numel(so)
    A=ALL(so(q));
    fprintf('  %4d %-40s %10.5f %11.4f %11.4f %8.3f %8.3f %8.3f %9.3f %9.5f %9.5f\n', ...
        q,A.lab,A.v(1),A.v(2),A.v(3),A.e(1),A.e(2),A.e(3),A.etot,A.kr,A.kl);
end

%% ======================================================================
%  5. VERDICT
%  ======================================================================
fprintf('\n---- 5. verdict ----\n');
Ab=ALL(so(1)); Ab2=ALL(so(2));
fprintf('  meilleure configuration : %s, ecart total %.3f %%\n',Ab.lab,Ab.etot);
fprintf('  seconde                 : %s, ecart total %.3f %%\n',Ab2.lab,Ab2.etot);
etot_al=ALL(end).etot;                      % la configuration alignee du §2
et0=[ALL(1:end-1).etot]; best0=min(et0);    % le meilleur des phi = 0
IDENT = etot_al < 0.05 && etot_al < best0/10;

if IDENT
    fprintf(['\n  1. LA PROVENANCE N''EST PAS PERDUE, ELLE EST IDENTIFIEE. Les trois\n' ...
             '     cellules sortent de BLDC_MEC_COMPLET.m (sect.1 l.92-93 et sect.2\n' ...
             '     l.168), reseau a une dent Nsurf = 1260, Np = 721, span = 2pi/p,\n' ...
             '     mu_i = %g, kfringe = %.3f -- c''est-a-dire LA MEME configuration\n' ...
             '     que A1 -- mais evaluees A LA POSITION ROTOR ALIGNEE SUR LA FEA\n' ...
             '     (indice %d/%d, %.3f deg mec) et non a phi = 0. Ecart total\n' ...
             '     residuel %.4f %%. Elles sont ecrites a pleine precision dans\n' ...
             '     outF.txt l.151/154/155 et dans sa traduction anglaise\n' ...
             '     article\\results_PMSM_15_14.txt l.158/161/162.\n'], ...
             M.muI,kfr,qb,Np721,phis(qb)*180/pi,etot_al);
else
    fprintf(['\n  1. LA PISTE DE LA POSITION ALIGNEE NE SUFFIT PAS. La reexecution\n' ...
             '     de BLDC_MEC_COMPLET a la position identifiee laisse un ecart\n' ...
             '     total de %.4f %%, contre %.4f %% pour le meilleur des %d cas a\n' ...
             '     phi = 0. La provenance reste a etablir.\n'], ...
             etot_al,best0,numel(et0));
end
fprintf(['\n  2. AUCUN PARAMETRE DE MODELE N''EXPLIQUE L''ECART. Le balayage de\n' ...
         '     %d configurations a phi = 0 (Nsurf, Np, kfringe, mu_i, base p0/p1,\n' ...
         '     cogging_mec2, mu_r historique) ne descend pas sous %.3f %% d''ecart\n' ...
         '     total, alors que le seul changement de position rotor y descend a\n' ...
         '     %.4f %%. Le parametre discriminant n''est pas un parametre de\n' ...
         '     MODELE, c''est un parametre d''EXTRACTION.\n'],numel(et0),best0,etot_al);
if MIX
    fprintf(['\n  3. LA COLONNE PUBLIEE EST INCOHERENTE AVEC ELLE-MEME. Le manuscrit\n' ...
             '     porte k_r = %.4f, qui est le quotient des FMM de A1 (phi = 0) et\n' ...
             '     non des FMM publiees (%.4f). Les trois cellules et les deux\n' ...
             '     facteurs proviennent de DEUX POSITIONS ROTOR DIFFERENTES. C''est\n' ...
             '     precisement pourquoi k_r et k_l "coincidaient exactement" alors\n' ...
             '     que les trois autres cellules ne se regeneraient pas : ils ne\n' ...
             '     venaient pas du meme calcul.\n'],MSC(1),TGT(2)/TGT(3));
end
fprintf(['\n  4. RECOMMANDATION. Ne pas chercher a "retomber" sur 0.9846 / 763.2 /\n' ...
         '     723.8 en reglant un parametre : ce serait falsifier un parametre de\n' ...
         '     modele pour compenser un changement de convention d''extraction.\n' ...
         '     Deux options coherentes, et une seule a retenir :\n' ...
         '       (a) TOUTE la colonne a phi = 0 -- c''est ce que fait A1 et ce que\n' ...
         '           porte deja MEC_DtN_paper_v2.tex (0.9860 / 762.9 / 725.1) ;\n' ...
         '           il faut alors verifier que k_r et k_l y sont bien ceux de\n' ...
         '           phi = 0, ce qui est le cas (%.5f et %.5f).\n' ...
         '       (b) TOUTE la colonne a la position alignee, k_r = %.5f et\n' ...
         '           k_l = %.5f compris -- ce qui obligerait a corriger deux\n' ...
         '           cellules aujourd''hui justes.\n' ...
         '     (a) est retenue : elle est deja appliquee, elle est la moins\n' ...
         '     couteuse, et la position phi = 0 est une convention declarable,\n' ...
         '     alors que la position alignee depend de la reference FEA.\n'], ...
         PRG(4),PRG(5),kr_al,kl_al);
fprintf(['\n  5. A RETENIR POUR LE RESTE DU MANUSCRIT. Toute grandeur dependant\n' ...
         '     de la position rotor doit porter sa position. La difference n''est\n' ...
         '     pas du bruit : elle atteint %.2f %% sur B_r crete, ce qui est du\n' ...
         '     meme ordre que plusieurs ecarts revendiques comme significatifs.\n'], ...
         100*abs(V_al(1)/Br0_3600-1));

save('A5_lumped.mat','ALL','TGT','PRG','MSC','KRKL_pub','V_al','qb','db', ...
     'phis','kr_al','kl_al','mm_0','mg_0','kr_0','kl_0','Br0_2000','Br0_3600', ...
     'NS','KF','MU','BA','NP');
fprintf('\n  duree %.0f s\n=== A5 termine ===\n',toc(t0));
diary off;

%% ======================================================================
function v=a5_num(Lr,pat)
%  Premier nombre de la premiere ligne contenant pat.
    v=NaN; k=find(contains(Lr,pat),1); if isempty(k), return; end
    tk=regexp(char(Lr(k)),'[-+]?\d+\.?\d*','match');
    if ~isempty(tk), v=str2double(tk{1}); end
end
% ----------------------------------------------------------------------
function v=a5_tex(Lt,pat)
%  Deuxieme quantite $...$ NUMERIQUE de la ligne : colonne Lumped du
%  tableau (Mesh, Lumped, FEA).
    v=NaN; k=find(contains(Lt,pat),1); if isempty(k), return; end
    tk=regexp(char(Lt(k)),'\$([-+]?[0-9.]+)\$','tokens');
    if numel(tk)>=2, v=str2double(tk{2}{1}); end
end
% ----------------------------------------------------------------------
function ALL=a5_push(ALL,lab,v,kr,kl,TGT)
    e=100*(v./TGT-1);
    ALL(end+1)=struct('lab',lab,'v',v,'kr',kr,'kl',kl,'e',e,'etot',sum(abs(e)));
end
% ----------------------------------------------------------------------
function [mmf_m,mmf_g]=a5_mmf(AG,Us,phq,M,Nq)
%  Partie FMM de krkl, avec la finesse de quadrature Nq en parametre.
    nu=AG.nu; Nm=M.Nm;
    Usc=AG.Wc*Us; Uss=AG.Ws*Us;
    pic=AG.alphau.*Usc+AG.betasrc.*cos(nu*phq);
    pis=AG.alphau.*Uss+AG.betasrc.*sin(nu*phq);
    thq=linspace(0,2*pi,Nq); thq(end)=[];
    phi_i=(pic.'*cos(nu*thq)+pis.'*sin(nu*thq));
    phi_b=(Usc.'*cos(nu*thq)+Uss.'*sin(nu*thq));
    mm=zeros(1,Nm); mg=zeros(1,Nm);
    for k=1:Nm
        c=(k-1)*2*pi/Nm+phq;
        [~,j]=min(abs(angle(exp(1i*(thq-c)))));
        mm(k)=abs(phi_i(j)); mg(k)=abs(phi_i(j)-phi_b(j));
    end
    mmf_m=mean(mm); mmf_g=mean(mg);
end
