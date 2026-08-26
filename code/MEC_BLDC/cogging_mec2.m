function R = cogging_mec2(M, Nsurf, nm, Np, muI, span)
%COGGING_MEC2  Idem cogging_mec, mais avec un SOUS-MODELE 2D RESOLU de la
%              BOUCHE D'ENCOCHE (a la place du scalaire identifie kfringe).
%
%   MOTIVATION. Dans cogging_mec, tout le comportement de l'ouverture
%   d'encoche est porte par UNE conductance de frange kfringe, identifiee sur
%   les harmoniques FEA. Cela laisse deux defauts :
%     * les bandes de denture restent a +/-25-32 % (ordres 8 et 23) ;
%     * un meme scalaire doit servir a deux geometries differentes (entree
%       RADIALE a vide / pont TANGENTIEL dent-a-dent en induit).
%   Ici la bouche est MAILLEE, avec sa geometrie reelle et AUCUN parametre
%   libre :
%     - canal d'ouverture : largeur ws0 constante sur la profondeur hs0 ;
%     - biseau : largeur ws0 -> ws1 sur la profondeur hs1 ;
%     - parois laterales = becs des deux dents adjacentes ;
%     - bec de dent (fer) represente par un noeud propre TT, relie au corps.
%   Le pont dent-a-dent EMERGE alors de la meme geometrie que l'entree
%   radiale : plus d'ambiguite, et plus de lambda_tip analytique a ajouter.
%
%   nm = nombre de couches radiales dans la bouche (0 -> voir cogging_mec).

mu0=4*pi*1e-7;
Rs=M.Rsi; L=M.ls; Ns=M.Ns; Nm=M.Nm; taus=2*pi/Ns;
if nargin<6||isempty(span), span=2*pi/Nm; end
numax=floor(Nsurf/2);

% ---- grille de surface : fer (face de dent) / air (ouverture) ----------
dth=2*pi/Nsurf; ths=(0:Nsurf-1)*dth; dths=dth*ones(1,Nsurf);
wtf=taus-M.ws0/Rs;
tooth=zeros(1,Nsurf); isFe=false(1,Nsurf); slot=zeros(1,Nsurf);
for j=1:Nsurf
    i0=mod(round(ths(j)/taus),Ns); dd=angle(exp(1i*(ths(j)-i0*taus)));
    tooth(j)=i0+1;
    if abs(dd)<wtf/2, isFe(j)=true;
    elseif dd>0,      slot(j)=i0+1;                 % entre dent i0+1 et i0+2
    else,             slot(j)=mod(i0-1+Ns,Ns)+1;    % entre dent i0 et i0+1
    end
end
AG=airgap_magnet(M,ths,dths,numax);

% ---- geometrie de la bouche ------------------------------------------
D=M.hs0+M.hs1; dd_=D/nm;                     % profondeur totale, pas radial
dc=((1:nm)-0.5)*dd_;                          % profondeur des centres
wof=@(d) M.ws0 + (d>M.hs0).*(M.ws1-M.ws0).*max(d-M.hs0,0)/M.hs1;
wc=wof(dc);                                   % largeur au centre des couches
wi=wof((1:nm-1)*dd_);                         % largeur aux interfaces

% ---- indexation --------------------------------------------------------
TT=@(i) Nsurf+i; TB=@(i) Nsurf+Ns+i; YY=@(i) Nsurf+2*Ns+i;
MM=@(k,j) Nsurf+3*Ns+(k-1)*nm+j;
Ntot=Nsurf+3*Ns+Ns*nm;
A=zeros(Ntot);
    function add(a,b,g), A(a,a)=A(a,a)+g; A(b,b)=A(b,b)+g; A(a,b)=A(a,b)-g; A(b,a)=A(b,a)-g; end

% ---- fer : face -> bec -> corps -> culasse ----------------------------
g_face=mu0*muI*(dth*Rs)*L/(D/2);                       % face de dent -> bec
w_shoe=0.5*((taus*Rs-M.ws0)+(2*pi*(Rs+D)/Ns-M.ws1));   % largeur moyenne de bec
G_shoe=mu0*muI*w_shoe*L/D;                             % bec -> corps de dent
Gt=mu0*muI*M.wst1*L/M.hs2;
Gy=mu0*muI*M.wsy*L/(taus*(M.Rso-M.wsy/2));
for j=1:Nsurf
    if isFe(j), add(j,TT(tooth(j)),g_face); end
end
for i=1:Ns
    add(TT(i),TB(i),G_shoe); add(TB(i),YY(i),Gt); add(YY(i),YY(mod(i,Ns)+1),Gy);
end

% ---- bouche d'encoche (air, maillee) ----------------------------------
for j=1:Nsurf                                          % surface -> 1re couche
    if ~isFe(j)
        add(j,MM(slot(j),1), mu0*(dth*Rs)*L/(dd_/2));
    end
end
for k=1:Ns
    for j=1:nm-1                                       % radial dans la bouche
        add(MM(k,j),MM(k,j+1), mu0*wi(j)*L/dd_);
    end
    for j=1:nm                                         % lateral vers les 2 becs
        Gl=2*mu0*dd_*L/wc(j);
        add(MM(k,j),TT(k),Gl); add(MM(k,j),TT(mod(k,Ns)+1),Gl);
    end
end

A(1:Nsurf,1:Nsurf)=A(1:Nsurf,1:Nsurf)-AG.Y;
ref=YY(1); keep=true(Ntot,1); keep(ref)=false;

% ---- balayage rotor ---------------------------------------------------
phis=linspace(0,span,Np);
nu=AG.nu; brm=AG.brm;
Isrc=zeros(Ntot,Np);
Isrc(1:Nsurf,:)=L*Rs*pi*( AG.Wc.'*(brm.*cos(nu*phis)) + AG.Ws.'*(brm.*sin(nu*phis)) );
U=zeros(Ntot,Np); U(keep,:)=A(keep,keep)\Isrc(keep,:);

R.PhiT=Gt*(U(TB(1:Ns),:)-U(YY(1:Ns),:));
R.PhiY=Gy*(U(YY(1:Ns),:)-U(YY([2:Ns 1]),:));

% ---- couple MST par harmoniques --------------------------------------
Us=U(1:Nsurf,:); Usc=AG.Wc*Us; Uss=AG.Ws*Us;
cphi=cos(nu*phis); sphi=sin(nu*phis);
Brc=AG.bru.*Usc+AG.brmq.*cphi;  Brs=AG.bru.*Uss+AG.brmq.*sphi;
Btc=-AG.btu.*Uss-AG.btmq.*sphi; Bts=AG.btu.*Usc+AG.btmq.*cphi;
T=(L*M.rmid^2/mu0)*pi*sum(Brc.*Btc+Brs.*Bts,1); T=T-mean(T);
Yf=abs(fft(T)); Yf=Yf(2:max(2,floor(Np/2))); [~,km]=max(Yf);
R.T=T; R.phis=phis; R.Tpp=(max(T)-min(T))*1e3; R.order=round(km*(2*pi/span));
R.U=U; R.AG=AG; R.Usurf=Us; R.Nsurf=Nsurf; R.nm=nm; R.Gt=Gt; R.Gy=Gy; R.span=span;

thq=linspace(0,2*pi,3601);
[Br,Bt]=AG.field(U(1:Nsurf,1),0,thq);
R.thq=thq; R.Br=Br; R.Bt=Bt;
R.Bmean=mean(abs(Br)); R.Bpeak=max(Br); R.Btrms=sqrt(mean(Bt.^2));
R.Bg1=hypot(2*mean(Br.*cos(M.p*thq)),2*mean(Br.*sin(M.p*thq)));
end
