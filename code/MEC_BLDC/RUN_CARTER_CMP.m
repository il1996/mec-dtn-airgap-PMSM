%% RUN_CARTER_CMP - ce que l'operateur DtN apporte face au MEC classique
%
%  L'article doit montrer non seulement que le DtN est proche de l'EF,
%  mais ce qu'il APPORTE face a la fermeture d'entrefer usuelle. On
%  compare donc TROIS niveaux de modele sur la MEME machine et la MEME
%  reference, du plus classique au plus complet :
%
%   (N1) ENTREFER LISSE + CARTER. Le stator est remplace par une surface
%        lisse infiniment permeable (potentiel de surface nul) et
%        l'entrefer geometrique g par g_eff = kC*g. C'est la fermeture
%        d'Ostovic et de la quasi-totalite des MEC : elle restitue le flux
%        moyen et rien d'autre.
%   (N2) LISSE + CARTER + FONCTION DE PERMEANCE. On module le champ lisse
%        par la permeance relative d'encoche lambda(theta) de forme
%        classique (creux sur l'ouverture, profondeur fixee par Carter).
%        C'est le niveau de Zarko-Ban-Lipo en version reelle.
%   (N3) DtN COMPLET. Le potentiel de surface est porte par le reseau de
%        reluctances : la denture EMERGE.
%
%  Grandeurs jugees : fondamental, les deux bandes laterales de denture
%  nu = |Ns -+ p| = 8 et 22, et la composante TANGENTIELLE, qu'aucun
%  reseau a branches radiales ne produit.
%  BLOC A6 (SPEC v3 §8) : mise en forme finale. Les colonnes gamma_m,
%  gamma_z et nombre d'elements exigees par A6 sont ajoutees ICI, dans la
%  chaine qui produit deja les six lignes -- et non dans un script
%  parallele, qui creerait une seconde chaine pour la meme grandeur.
clear; clc;
diary('A6_n2b_out.txt'); diary on;
M=machine_bldc(); Ns=M.Ns; p=M.p; g=M.lag; taus=2*pi/Ns;
if isfile('kfringe_ident.mat'), S0=load('kfringe_ident.mat'); kfr=S0.kbest; else, kfr=0.325; end
fprintf('=== A6 : comparaison de fermetures, mise en forme finale ===\n');
fprintf('  machine  : PMSM 15/14, 750 W | Ns = %d, p = %d\n',Ns,p);
fprintf('  reference: %s\n',M.FEA.dir);
fprintf('  chaine   : RUN_CARTER_CMP (N1, N2, trois N2b, N3)\n');
fprintf('  kfringe  : %.4f | solveur lineaire mu_r = %g\n\n',kfr,M.muI);
rd=@(f)readmatrix(f,'FileType','text','Delimiter','\t','NumHeaderLines',1);
fea=M.FEA.dir;
d4=rd(fullfile(fea,'magnetostique(Magnetic_loading)','Calculator Expressions Plot 4.tab'));
d2=rd(fullfile(fea,'magnetostique(Magnetic_loading)','Calculator Expressions Plot 2.tab'));
thu=linspace(0,2*pi,3601); thu(end)=[];
BrF=interp1([d4(:,2)*pi/180;2*pi],[d4(:,3);d4(1,3)],thu,'linear','extrap');
BtF=interp1(linspace(0,2*pi,numel(d2(:,2))).',d2(:,2),thu,'linear','extrap');
amp=@(B,th,k) 2*abs(mean(B(:).'.*exp(-1i*k*th(:).')));

%% ---------- coefficient de Carter -------------------------------------
%  gamma = (4/pi)[ (b/2g) atan(b/2g) - ln(sqrt(1+(b/2g)^2)) ]
%  L'entrefer magnetique inclut l'aimant : g_m = g + hm/mu_rec.
%  ATTENTION AU NOM : g_mag est l'ENTREFER MAGNETIQUE (metres). Il n'a
%  RIEN A VOIR avec le gamma_m de la fonction de Silva plus bas, qui est un
%  ARC (radians). Les deux sont notes "gm" dans la litterature ; ici ils
%  portent des noms distincts pour qu'aucune ligne ne puisse les confondre.
g_mag=g+M.hm/M.mu_r; x=M.ws0/(2*g_mag);
gam=(4/pi)*( x*atan(x) - log(sqrt(1+x^2)) );
kC=taus*M.Rsi/(taus*M.Rsi - gam*g_mag);
fprintf('=== Ce que le DtN apporte face au MEC classique ===\n');
fprintf('  ouverture %.2f mm, entrefer magnetique %.3f mm -> kC = %.4f\n', ...
    M.ws0*1e3,g_mag*1e3,kC);

%% ---------- N1 : entrefer lisse + Carter -------------------------------
%  Surface lisse infiniment permeable : Us = 0, donc Br = brm(nu).
%  L'effet d'encoche est porte par le seul allongement g -> kC*g.
M1=M; M1.lag=kC*g; M1.Rsi=M.Rro+M1.lag; M1.rmid=M1.Rsi-M1.lag/2;
ths=(0:1259)*2*pi/1260; AG1=airgap_magnet(M1,ths,(2*pi/1260)*ones(1,1260),630);
[Br1,Bt1]=AG1.field(zeros(1260,1),0,thu);

%% ---------- N2 : + fonction de permeance d'encoche ---------------------
%  Permeance relative classique : creux en cos^2 sur l'ouverture, de
%  profondeur telle que la moyenne redonne 1/kC.
lam=ones(size(thu));
for k=0:Ns-1
    c=(k+0.5)*taus; d=abs(angle(exp(1i*(thu-c))));
    m=d<=(M.ws0/M.Rsi)/2;
    lam(m)=1-0.5*(1+cos(pi*d(m)/((M.ws0/M.Rsi)/2)));
end
%  NORMALISATION. lam/mean(lam) impose <lam> = 1 (les deux kC s'annulent).
%  CONSEQUENCE A DECLARER : la modulation a la periode du pas dentaire, donc
%  elle ne cree d'harmoniques qu'aux rangs p +- k*Ns, JAMAIS au rang p. Avec
%  une moyenne unitaire, N2 rend donc EXACTEMENT le meme Bg1 que N1 -- ce
%  n'est pas une coincidence numerique mais une propriete de construction.
lam=lam/mean(lam);                       % normalisation : <lam> = 1
Br2=Br1.*lam; Bt2=Bt1;                   % la permeance reelle ne cree pas de Bt

%% ---------- N2b : la VRAIE fonction de Silva et al. --------------------
%  T6. La ligne N2 ci-dessus utilise un creux en cos^2 GENERIQUE, non
%  identifie : ce n'est pas un test loyal de la famille (iii). La fonction
%  reellement publiee par Silva et al. (eq. 14) est un PLATEAU suivi d'une
%  DECROISSANCE GAUSSIENNE, reprise de Lannoo [34] sur la base d'Ostovic [3] :
%
%      Lambda(sigma) = Lambda_max                                si sigma <= gm
%                    = Lambda_max * exp(-[(sigma-gm)/(gz-gm)]^2) sinon
%
%  C1 au raccord (derivee nulle des deux cotes de sigma = gm).
%
%  DEFAUT M1, ET C'EST LE POINT. gm et gz ne sont PAS des grandeurs
%  physiques : ce sont les demi-sommes et demi-differences des ARCS
%  D'ELEMENTS des deux surfaces. La fonction depend donc du maillage, et
%  tout raffinement tangentiel impose une recalibration. On le montre en
%  balayant l'arc d'element rotorique a_r, seul parametre libre.
a_t = (taus - M.ws0/M.Rsi);              % arc de FACE de dent [rad]
silva=@(d,gmm,gzz) (abs(d)<=gmm) + ...
      (abs(d)>gmm).*exp(-((abs(d)-gmm)./max(gzz-gmm,eps)).^2);
arList = a_t*[0.25 0.50 1.00];           % arcs d'element rotorique testes
Br2b=cell(1,numel(arList)); lamS=cell(1,numel(arList));
gmv=zeros(1,numel(arList)); gzv=gmv;     % gamma_m, gamma_z de chaque variante
for ia=1:numel(arList)
    a_r=arList(ia);
    gmv(ia)=abs(a_t-a_r)/2;  gzv(ia)=(a_t+a_r)/2;   % construction d'Ostovic
    L=zeros(size(thu));
    for k=0:Ns-1
        c=k*taus;                        % centre de FACE de dent
        d=angle(exp(1i*(thu-c)));
        L=L+silva(d,gmv(ia),gzv(ia));
    end
    L=L/mean(L);                         % meme normalisation que N2
    lamS{ia}=L; Br2b{ia}=Br1.*L;
end

%% ---------- A6 : les parametres de Silva, explicites -------------------
%  C'est ce qui rend le defaut de DETERMINATION lisible : gamma_m et
%  gamma_z ne sont pas des grandeurs physiques mais des demi-differences
%  et demi-sommes d'ARCS D'ELEMENTS. Publier le tableau sans ces colonnes
%  laisserait croire a un modele parametre par la geometrie.
%  A6.1 -- les deux parametres, en radians ET en degres, a cote du nombre
%  d'elements 2*pi/a_r dont ils sortent. Mises l'une a cote de l'autre, ces
%  colonnes disent tout : gamma_m et gamma_z SUIVENT le compte d'elements et
%  rien d'autre. Le nombre d'elements n'est meme pas entier (69.58, 34.79,
%  17.40 par tour) : c'est le signe qu'il ne s'agit pas d'un decoupage de la
%  machine mais d'un choix de discretisation.
fprintf('\n  --- A6.1 : parametres de la fonction de Silva, variante par variante ---\n');
fprintf('  reference stator : arc de FACE DE DENT a_t = %.6f rad = %.4f deg\n', ...
    a_t,a_t*180/pi);
fprintf('                     soit 2*pi/a_t = %.2f faces de dent par tour\n',2*pi/a_t);
fprintf('  construction d''Ostovic : gamma_m = |a_t - a_r|/2 , gamma_z = (a_t + a_r)/2\n\n');
%  Les valeurs imprimees ici sont gmv/gzv, c'est-a-dire EXACTEMENT les
%  arguments passes a silva() dans la boucle N2b ci-dessus -- pas une
%  seconde evaluation de la formule. Le tableau ne peut donc pas decrire
%  une fonction autre que celle qui a produit les lignes N2b.
fprintf('  %-16s %10s %10s %11s %11s %11s %11s %11s\n', ...
    'variante','a_r (rad)','a_r (deg)','2pi/a_r', ...
    'gam_m (rad)','gam_m (deg)','gam_z (rad)','gam_z (deg)');
for ia=1:numel(arList)
    fprintf('  %-16s %10.6f %10.4f %11.2f %11.6f %11.4f %11.6f %11.4f\n', ...
        sprintf('a_r = %.2f a_t',arList(ia)/a_t),arList(ia),arList(ia)*180/pi, ...
        2*pi/arList(ia),gmv(ia),gmv(ia)*180/pi,gzv(ia),gzv(ia)*180/pi);
end
%  Largeur relative du PLATEAU : gamma_m/gamma_z. Elle dit, sans commentaire,
%  quelle part de la fonction (14) est encore un plateau et quelle part est
%  deja la gaussienne. Elle passe de 60 % a 0 % sans que la machine change.
fprintf('\n  %-16s %14s %16s\n','variante','gam_m / gam_z','forme obtenue');
for ia=1:numel(arList)
    if gmv(ia)>0, frm='plateau + gaussienne'; else, frm='gaussienne seule'; end
    fprintf('  %-16s %13.1f %% %2s%s\n', ...
        sprintf('a_r = %.2f a_t',arList(ia)/a_t),100*gmv(ia)/gzv(ia),'',frm);
end
fprintf(['\n  LECTURE DE CE TABLEAU. gamma_m va de %.4f deg a 0 et gamma_z de\n' ...
         '  %.4f a %.4f deg alors que la MACHINE n''a pas bouge d''un micron :\n' ...
         '  seul le nombre d''elements rotoriques par tour a change (%.2f -> %.2f).\n' ...
         '  Ni gamma_m ni gamma_z ne se deduisent donc de la geometrie ; ils se\n' ...
         '  deduisent du MAILLAGE des deux surfaces. Le cas a_r = a_t donne meme\n' ...
         '  gamma_m = 0 exactement : le plateau disparait et la fonction (14) se\n' ...
         '  reduit a la seule gaussienne. Tout raffinement tangentiel change donc\n' ...
         '  la fermeture d''entrefer, et impose une recalibration.\n'], ...
    gmv(1)*180/pi,gzv(1)*180/pi,gzv(3)*180/pi,2*pi/arList(1),2*pi/arList(3));

%  A6.2 -- <lambda> = 1 pour les QUATRE modulations. C'est la condition qui
%  rend le rang p intouchable (voir le bloc final) ; elle doit etre VERIFIEE
%  et non postulee, donc imprimee a 6 decimales avec son ecart a l'unite.
fprintf('\n  --- A6.2 : moyenne des QUATRE modulations d''entrefer ---\n');
fprintf('  %-30s %12s %14s\n','modulation','<lambda>','|<lambda> - 1|');
lamAll=[{lam} lamS]; lamNom={'N2   creux cos^2 generique', ...
    sprintf('N2b  Silva a_r = %.2f a_t',arList(1)/a_t), ...
    sprintf('N2b  Silva a_r = %.2f a_t',arList(2)/a_t), ...
    sprintf('N2b  Silva a_r = %.2f a_t',arList(3)/a_t)};
mL=cellfun(@mean,lamAll);
for ia=1:4
    fprintf('  %-30s %12.6f %14.2e\n',lamNom{ia},mL(ia),abs(mL(ia)-1));
end
fprintf('  ecart maximal a l''unite sur les quatre : %.2e (arrondi machine)\n', ...
    max(abs(mL-1)));

%% ---------- N3 : DtN complet (reseau) ---------------------------------
R3=cogging_mec(M,1260,0,3,M.muI,kfr,1e-6);
Br3=R3.Br(1:end-1); Bt3=R3.Bt(1:end-1); th3=R3.thq(1:end-1);

%% ---------- comparaison ------------------------------------------------
q=@(B,th,k)amp(B,th,k);
tF=[q(BrF,thu,p) q(BrF,thu,8) q(BrF,thu,22) sqrt(mean(BtF.^2))];
%  Les six lignes sont EVALUEES UNE SEULE FOIS, dans V (6 x 4). L'impression,
%  la dispersion A6.3 et le fichier .mat lisent tous cette meme matrice : il
%  n'existe aucun chemin par lequel le tableau publie et le chiffre de
%  dispersion pourraient diverger.
mdl={'(N1) lisse + Carter','(N2) + permeance cos^2 generique'};
V=zeros(6,4);
V(1,:)=[q(Br1,thu,p) q(Br1,thu,8) q(Br1,thu,22) sqrt(mean(Bt1.^2))];
V(2,:)=[q(Br2,thu,p) q(Br2,thu,8) q(Br2,thu,22) sqrt(mean(Bt2.^2))];
for ia=1:numel(arList)
    B=Br2b{ia};
    V(2+ia,:)=[q(B,thu,p) q(B,thu,8) q(B,thu,22) sqrt(mean(Bt1.^2))];
    mdl{2+ia}=sprintf('(N2b) Silva, a_r = %.2f a_t',arList(ia)/a_t);
end
V(6,:)=[q(Br3,th3,p) q(Br3,th3,8) q(Br3,th3,22) sqrt(mean(Bt3.^2))];
mdl{6}='(N3) DtN complet (reseau)';
fprintf('\n  %-34s %9s %10s %9s %10s %9s %10s %8s %10s\n', ...
    'modele','Bg1 (T)','ecart','rang 8','ecart','rang 22','ecart','Bt rms','ecart');
for ir=1:6
    v=V(ir,:);
    fprintf('  %-34s %9.4f %8.1f %% %9.5f %8.1f %% %9.5f %8.1f %% %8.4f %8.1f %%\n', ...
        mdl{ir},v(1),100*(v(1)-tF(1))/tF(1),v(2),100*(v(2)-tF(2))/tF(2), ...
        v(3),100*(v(3)-tF(3))/tF(3),v(4),100*(v(4)-tF(4))/tF(4));
end
fprintf('  %-34s %9.4f %10s %9.5f %10s %9.5f %10s %8.4f %10s\n', ...
    'reference FEA',tF(1),'--',tF(2),'--',tF(3),'--',tF(4),'--');

fprintf('\n  --- lecture ---\n');
fprintf('  N1 ne produit AUCUNE bande de denture : par construction, un entrefer\n');
fprintf('     lisse n''a pas d''harmonique d''encoche. Le Bt qui subsiste (%.3f T)\n',sqrt(mean(Bt1.^2)));
fprintf('     vient des transitions entre aimants, PAS de la denture : d''ou -20 %%.\n');
fprintf('  N2 restitue des bandes laterales, mais la PROFONDEUR du creux de\n');
fprintf('     permeance n''est pas identifiee ici (forme cos^2 generique) et elle\n');
fprintf('     depasse la cible d''un facteur 3 a 4. C''est precisement le point :\n');
fprintf('     ce niveau de modele EXIGE un calage, que N3 n''exige pas. Et une\n');
fprintf('     permeance REELLE module l''amplitude sans rien ajouter a Bt.\n');
fprintf('  N2b est la fonction REELLE de Silva (plateau + gaussienne, eq. 14).\n');
fprintf('     Ses parametres gm et gz sont des ARCS D''ELEMENTS, pas des grandeurs\n');
fprintf('     physiques : le balayage de a_r ci-dessus mesure cette dependance.\n');
fprintf('  N3 produit les bandes ET la denture dans Bt, sans coefficient de Carter,\n');
fprintf('     sans fonction de permeance et sans aucune constante identifiee.\n');

%  A6.3 -- la dispersion, chiffree. C'est la mesure du defaut : la meme
%  machine, la meme fonction publiee, trois maillages -> trois reponses.
a8=V(3:5,2).';                           % relu dans V, pas re-evalue
fprintf('\n  --- A6.3 : dispersion du rang 8 sur les trois variantes de Silva ---\n');
fprintf('  %-22s %11s %12s %12s\n','variante','2pi/a_r','rang 8 (T)','ecart / FEA');
for ia=1:numel(arList)
    fprintf('  %-22s %11.2f %12.5f %11.1f %%\n', ...
        sprintf('a_r = %.2f a_t',arList(ia)/a_t),2*pi/arList(ia), ...
        a8(ia),100*(a8(ia)-tF(2))/tF(2));
end
fprintf('  min %.5f T | max %.5f T | moyenne %.5f T | reference FEA %.5f T\n', ...
    min(a8),max(a8),mean(a8),tF(2));
fprintf('  DISPERSION (max - min)/moyenne = %.0f %%\n',100*(max(a8)-min(a8))/mean(a8));
fprintf(['  Un seul parametre de MAILLAGE, non contraint par la machine, fait\n' ...
         '  varier la bande de denture de %.0f %% : la famille (iii) n''a pas de\n' ...
         '  determination propre. Le DtN (N3), lui, donne %.1f %% sans reglage.\n'], ...
    100*(max(a8)-min(a8))/mean(a8),100*(q(Br3,th3,8)-tF(2))/tF(2));
fprintf('\n  --- pourquoi N1 et N2/N2b partagent le meme Bg1 ---\n');
fprintf('  Toute modulation de periode le pas dentaire ne cree d''harmoniques\n');
fprintf('  qu''aux rangs p +- k*Ns. Normalisee a <lam> = 1, elle laisse le rang p\n');
fprintf('  INCHANGE. L''egalite est structurelle, pas fortuite. Verification :\n');
fprintf('  <lam> N2 = %.6f | N2b = %.6f %.6f %.6f  (A6.2)\n',mL(1),mL(2),mL(3),mL(4));
fprintf('  ecart maximal des SIX Bg1 modulees a Bg1(N1) : %.2e T\n', ...
    max(abs(V(2:5,1)-V(1,1))));

%  SAUVEGARDE. Tout ce que l'article a besoin de citer part d'ici : la
%  matrice V (six lignes x quatre grandeurs), la reference tF, les
%  parametres de Silva gmv/gzv et les moyennes mL. Aucun de ces chiffres
%  n'aura donc a etre retape.
res=struct('modeles',{mdl},'V',V,'colonnes',{{'Bg1','rang8','rang22','Btrms'}}, ...
    'tF',tF,'arList',arList,'a_t',a_t,'gamma_m',gmv,'gamma_z',gzv, ...
    'nelem',2*pi./arList,'lam_moy',mL,'kC',kC, ...
    'disp8_pct',100*(max(a8)-min(a8))/mean(a8));
save('carter_cmp.mat','kC','tF','arList','a_t','gmv','gzv','mL','V','mdl','res','-v7.3');
fprintf('\n=== A6 termine ===\n');
diary off;
