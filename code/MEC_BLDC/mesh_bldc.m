function ME = mesh_bldc(M, Ms, nst, nys, iabc, phi_rot, kfringe, muI, nsh)
%MESH_BLDC  Reseau de reluctances 2D RAFFINE du stator BLDC + DtN etendu.
%
%   POURQUOI. Le reseau BLDC utilise jusqu'ici UNE branche par dent et UNE
%   par segment de culasse : 30 degres de liberte dans le fer. Les grandeurs
%   globales (flux, FEM, couple, inductances) en sortent justes -- elles ne
%   dependent que des flux de branche -- mais aucune carte de champ ne peut
%   etre plus fine que cela. Le programme du moteur asynchrone (MEC_IM,
%   mec.mesh_refined) resout le meme probleme par un MAILLAGE POLAIRE : Ms
%   colonnes angulaires x Ls couches radiales, chaque cellule portant quatre
%   branches (radiales haut/bas, tangentielles gauche/droite) -- l'element de
%   Perho/Silva. C'est cette finesse, et non le trace, qui rend ses cartes
%   comparables a l'EF. On la porte ici.
%
%   CE QUI CHANGE PAR RAPPORT A L'IM. L'IM couple stator et rotor par des
%   permeances dent-a-dent (noyau gaussien identifie par EF). Le BLDC dispose
%   de mieux : l'operateur DtN ETENDU (airgap_magnet) donne la solution
%   EXACTE de Laplace sur les deux couches entrefer + aimant, avec la
%   magnetisation en source. On garde donc le DtN sur les noeuds de SURFACE
%   et on ne raffine que le fer statorique -- la ou le modele manquait de
%   degres de liberte.
%
%   ME = mesh_bldc(M,Ms,nst,nys,iabc,phi_rot,kfringe,muI)
%     Ms      : colonnes angulaires (grille UNIFORME, requise par le DtN)
%     nst,nys : couches radiales dans les dents / dans la culasse
%     iabc    : courants de phase [A] (3x1), [] a vide
%     phi_rot : position rotor [rad]
%     muI     : permeabilite relative du fer, ou 'nl' pour la courbe B(H)
%
%   Sortie ME : listes de branches (a,b,iron,l,A,P,E), la matrice DtN sur
%   les noeuds de surface, et toutes les metadonnees geometriques
%   necessaires au trace des cartes.

mu0=4*pi*1e-7;
if nargin<2||isempty(Ms),   Ms=360;  end
if nargin<3||isempty(nst),  nst=4;   end
if nargin<4||isempty(nys),  nys=3;   end
if nargin<5,                iabc=[]; end
if nargin<6||isempty(phi_rot), phi_rot=0; end
if nargin<7||isempty(kfringe), kfringe=0.325; end
if nargin<8||isempty(muI),  muI=M.muI; end
Ns=M.Ns; L=M.ls; Ki=M.Ki; taus=2*pi/Ns; Ntc=M.Ntc;
PH={[1 -2 -15 3 14],[6 -7 -5 8 4],[11 -12 -10 13 9]};
Rs=M.Rsi; Rso=M.Rso; hs=M.hs0+M.hs1+M.hs2;
Ls=nst+nys;

%% ---------- grille angulaire UNIFORME (imposee par le DtN) -------------
dth=2*pi/Ms; ths=(0:Ms-1)*dth;
tooth=zeros(1,Ms); dd=zeros(1,Ms);
for c=1:Ms
    i0=mod(round(ths(c)/taus),Ns); tooth(c)=i0+1;
    dd(c)=angle(exp(1i*(ths(c)-i0*taus)));    % ecart au centre de dent
end

%% ---------- couches radiales CONFORMES A LA GEOMETRIE -------------------
%  Les frontieres de couches tombent sur les changements de profil : isthme
%  (hs0), biseau (hs1), puis nst couches dans le corps d'encoche (hs2) et
%  nys dans la culasse. Un decoupage uniforme sur hs melangerait l'isthme
%  (1 mm) et le corps (19.6 mm) dans une meme couche, ce qui remplacait la
%  hauteur d'ouverture reelle par 5 mm d'air.
%  nsh = nombre de couches PAR sous-region de bec (isthme hs0, biseau hs1).
%  C'est le parametre critique : la saturation des cornes de bec se
%  concentre sur quelques dixiemes de millimetre et c'est elle qui creuse
%  la permÃ©ance d'entrefer, donc la bande de denture.
if nargin<9||isempty(nsh), nsh=1; end
re=[Rs+(0:nsh)/nsh*M.hs0, ...
    Rs+M.hs0+(1:nsh)/nsh*M.hs1, ...
    Rs+M.hs0+M.hs1+(1:nst)/nst*M.hs2, ...
    Rs+hs+(1:nys)/nys*(Rso-Rs-hs)];
rc=0.5*(re(1:end-1)+re(2:end)); rth=diff(re);
nshT=2*nsh; Ls=nshT+nst+nys;                 % bec + corps + culasse

%% ---------- profil reel : DENT A FLANCS PARALLELES ----------------------
%  La dent de cette machine est a flancs PARALLELES (wst1 constant) : son
%  demi-arc a un rayon r vaut asin(wst1/2r) et NON (taus - b)/2 -- une
%  soustraction directe de la largeur d'encoche a l'arc du pas sous-estime
%  fortement le fer en profondeur (elle donnait 4.0 mm de dent au lieu de
%  8.49, d'ou Bg1 a -19 % et une modulation d'encoche multipliee par 10).
%  Dans le BEC (hs0 + hs1) c'est l'ouverture qui decoupe : demi-arc de fer
%  = taus/2 - asin(ws0/2r).
isFe=false(Ls,Ms);
for la=1:Ls
    if la>nshT+nst, isFe(la,:)=true; continue; end    % culasse : tout fer
    r=rc(la);
    if la<=nshT, half=taus/2-asin(min(M.ws0/(2*r),1));% bec : isthme/biseau
    else,        half=asin(min(M.wst1/(2*r),1));      % corps : flancs paralleles
    end
    isFe(la,:)=abs(dd)<half;
end
nst=nshT+nst;                                % couches "region dents" au total

%% ---------- FMM de bobinage par colonne ---------------------------------
AT=zeros(Ns,3);
for k=1:3
    for c=PH{k}
        i=abs(c); s=sign(c);
        AT(mod(i-2,Ns)+1,k)=AT(mod(i-2,Ns)+1,k)+s*Ntc;
        AT(i,k)            =AT(i,k)            -s*Ntc;
    end
end
Fcol=[];
if ~isempty(iabc)
    ATs=AT*iabc(:);                                     % At par encoche
    %  BOBINAGE DENTAIRE : la bobine entoure SA dent. Sa FMM Ntc*i
    %  s'applique donc LE LONG DE LA DENT, et non comme un cumul
    %  d'ampere-tours au niveau de l'alesage. C'est la convention de
    %  inductance_mec (Fs = sign*Ntc*I0 dans la branche de dent), la seule
    %  qui soit non ambigue ici : injecter le cumul Fcol donnait une source
    %  de plusieurs centaines d'ampere-tours par dent au lieu de 152, d'ou
    %  L_a a 1140 mH. La FMM est repartie sur la HAUTEUR DE BOBINE, donc
    %  sur les couches du corps d'encoche.
    Fc=zeros(Ns,1);
    for k=1:3
        for c=PH{k}, Fc(abs(c))=Fc(abs(c))+sign(c)*Ntc*iabc(k); end
    end
    Fcol=Fc.';                                  % FMM par DENT [1 x Ns]
end
%  couches du corps d'encoche : celles ou se trouve la bobine
labod=find(rc>Rs+M.hs0+M.hs1 & rc<Rs+hs);
if isempty(labod), labod=nsh+1; end

%% ---------- assemblage des branches -------------------------------------
Sid=@(c,la)(la-1)*Ms+c; Ng0=Ms*Ls; Su=@(c)Ng0+c;
N=Ng0+Ms;
maxB=3*Ms*Ls+2*Ms;
a=zeros(maxB,1); b=a; irn=a; ll=a; AA=a; PP=a; EE=a; nb=0;
%  metadonnees de branche : type (1 radial, 2 tangentiel, 3 surface),
%  colonne et couche. Elles servent a extraire les flux de DENT et de
%  CULASSE, donc le flux totalise, la FEM et les pertes fer.
tp=zeros(maxB,1); cl=tp; ly=tp; ctyp=0; ccol=0; clay=0;
    function push(n1,n2,ir,len,area,perm,emf)
        nb=nb+1; a(nb)=n1; b(nb)=n2; irn(nb)=ir;
        ll(nb)=len; AA(nb)=area; PP(nb)=perm; EE(nb)=emf;
        tp(nb)=ctyp; cl(nb)=ccol; ly(nb)=clay;
    end
for c=1:Ms
    ccol=c;
    % --- branches RADIALES (couche la -> la+1) ---
    ctyp=1;
    for la=1:Ls-1
        clay=la;
        len=rc(la+1)-rc(la); rm=0.5*(rc(la)+rc(la+1));
        area=dth*rm*Ki*L;
        if isFe(la,c)&&isFe(la+1,c)
            %  INJECTION DE LA FMM DE BOBINAGE. Fcol(c) est le potentiel
            %  magnetique scalaire du bobinage a la colonne c (cumul des
            %  ampere-tours, a moyenne nulle). On l'injecte dans la branche
            %  RADIALE de la couche d'entrefer, comme le fait mec.mesh_refined
            %  de MEC_IM. La graduation en profondeur de la fuite d'encoche
            %  n'a pas a etre imposee : elle EMERGE du reseau, le flux devant
            %  venir de la dent pour traverser l'encoche a chaque profondeur.
            emf=0;
            if ~isempty(Fcol)&&any(la==labod)
                emf=Fcol(tooth(c))/numel(labod);   % FMM repartie sur la bobine
            end
            push(Sid(c,la),Sid(c,la+1),1,len,area,0,emf);
        else
            push(Sid(c,la),Sid(c,la+1),0,0,0,mu0*(dth*rm*L)/len,0);
        end
    end
    % --- branches TANGENTIELLES (colonne c -> c+1) par couche ---
    c2=mod(c,Ms)+1; ctyp=2;
    for la=1:Ls
        clay=la; arc=dth*rc(la); area=rth(la)*Ki*L;
        %  AUCUNE FMM sur les branches tangentielles. En y mettant
        %  Fcol(c2)-Fcol(c) a TOUTES les couches, la source devenait un
        %  pur gradient : sa circulation sur toute boucle fermee etait
        %  nulle (les contributions tangentielles de deux couches
        %  differentes s'annulent), donc elle n'excitait rien. C'est ce qui
        %  ramenait L_a a 12.2 mH au lieu de 50.2.
        emf=0;
        if isFe(la,c)&&isFe(la,c2)
            push(Sid(c,la),Sid(c2,la),1,arc,area,0,emf);
        else
            push(Sid(c,la),Sid(c2,la),0,0,0,mu0*(rth(la)*L)/arc,emf);
        end
    end
    % --- demi-branche vers le noeud de SURFACE (r = Rs) ---
    ctyp=3; clay=0; len=rth(1)/2; rm=Rs+rth(1)/4;
    if isFe(1,c)
        push(Sid(c,1),Su(c),1,len,dth*rm*Ki*L,0,0);
    else
        %  Sous l'ouverture d'encoche : AIR pur. Pas de kfringe ici -- dans
        %  le reseau a une dent, kfringe representait forfaitairement toute
        %  l'entree radiale par la bouche d'encoche, faute de maillage. Le
        %  maillage la resout maintenant explicitement (couches d'isthme et
        %  de biseau, branches tangentielles de fuite) : le remettre
        %  reviendrait a compter deux fois la meme frange.
        push(Sid(c,1),Su(c),0,0,0,mu0*(dth*rm*L)/len,0);
    end
end
a=a(1:nb); b=b(1:nb); irn=logical(irn(1:nb));
ll=ll(1:nb); AA=AA(1:nb); PP=PP(1:nb); EE=EE(1:nb);
tp=tp(1:nb); cl=cl(1:nb); ly=ly(1:nb);
%  couche de reference pour le flux de DENT : mi-corps d'encoche (la ou la
%  dent est a flancs paralleles, donc de section exacte)
[~,lamT]=min(abs(rc-(Rs+M.hs0+M.hs1+M.hs2/2)));
%  couches de CULASSE : au-dela de la region dents
layY=(nst+1):Ls;

%% ---------- operateur DtN etendu sur les noeuds de surface --------------
AG=airgap_magnet(M,ths,dth*ones(1,Ms),floor(Ms/2));
Isrc=zeros(N,1);
Isrc(Ng0+(1:Ms))=L*Rs*pi*( AG.Wc.'*(AG.brm.*cos(AG.nu*phi_rot)) ...
                         + AG.Ws.'*(AG.brm.*sin(AG.nu*phi_rot)) );

ME=struct('a',a,'b',b,'iron',irn,'l',ll,'A',AA,'P',PP,'E',EE, ...
    'N',N,'Ms',Ms,'Ls',Ls,'nst',nst,'nys',nys,'ths',ths,'dth',dth, ...
    're',re,'rc',rc,'rth',rth,'isFe',isFe,'tooth',tooth,'Sid',Sid, ...
    'Su',Su,'Ng0',Ng0,'AG',AG,'Isrc',Isrc,'muI',muI,'phi',phi_rot, ...
    'Fcol',Fcol,'iabc',iabc,'kfringe',kfringe, ...
    'type',tp,'col',cl,'lay',ly,'lamT',lamT,'layY',layY,'taus',taus);
end
