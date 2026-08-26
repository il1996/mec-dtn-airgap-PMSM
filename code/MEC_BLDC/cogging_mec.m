function R = cogging_mec(M, Nsurf, numax, Np, muI, kfringe, span, basis)
%
%   basis (optionnel, defaut 'p0') : base de surface, transmise a
%   airgap_magnet. 'p1' = chapeau, projection en 1/nu^2 (tache V1).
%   ATTENTION : numax reste impose a Nsurf/2 quelle que soit la base. La
%   troncature n'est donc PAS un parametre libre de cette chaine -- elle
%   est verrouillee sur la discretisation de surface, contrairement a
%   mec.airgap_dtn_tooth cote machine asynchrone.
%COGGING_MEC  Champ dente + couple de detente : DtN etendu couple au stator.
%
%   R = cogging_mec(M,Nsurf,numax,Np,muI,kfringe,span)
%     Nsurf   : noeuds de surface stator (Nyquist : Nsurf > 2*LCM pour la
%               detente d'ordre 210 -> viser Nsurf >= 840)
%     numax   : IGNORE (impose a Nsurf/2, voir ci-dessous)
%     Np      : positions rotor sur le balayage
%     muI     : permeabilite relative du fer (modele lineaire)
%     kfringe : frange de l'ouverture d'encoche (0 = isolante ; identifie a
%               0.75 par RUN_SLOT_IDENT sur les bandes de denture FEA)
%     span    : etendue du balayage rotor [rad mec] (defaut : 1 pas polaire ;
%               utiliser 2*pi/p pour une periode ELECTRIQUE complete)
%
%   Sorties : R.Tpp/R.T/R.order (detente), R.Br/R.Bt/R.thq (profils),
%   R.U (potentiels), R.PhiT/R.PhiY (flux de dent / de culasse par position),
%   R.Bg1/R.Bmean/R.Bpeak/R.Btrms, R.err_dir (controle croise).

mu0=4*pi*1e-7;
Rs=M.Rsi; L=M.ls; Ns=M.Ns; Nm=M.Nm; taus=2*pi/Ns;
if nargin<7 || isempty(span), span=2*pi/Nm; end
% CONDITIONNEMENT : le potentiel de surface a Nsurf degres de liberte ; si
% numax < Nsurf/2 l'operateur de couronne Y est de RANG DEFICIENT et les
% noeuds d'ouverture d'encoche (non tenus par le reseau quand kfringe=0)
% tombent dans son noyau -> systeme singulier (RCOND ~ 1e-18 constate).
%
% DEVERROUILLAGE (bloc R4). La contrainte porte sur numax INFERIEUR a
% Nsurf/2. Rien n'interdit de monter AU-DESSUS : le rang reste plein, et la
% projection W(n) = (2/n.pi) sin(n.d/2) cos(n.th) s'ANNULE aux multiples de
% Nsurf (sin(k.pi) = 0), donc aucun repliement catastrophique. Monter numax
% a pavage FIXE libere donc la troncature du pavage -- ce que le §3.5 de
% l'Article I annonce comme falsifiable et ne testait pas.
% Le 3e argument, jusqu'ici ignore, est desormais HONORE s'il vaut au moins
% Nsurf/2 ; sinon la valeur sure est imposee, avec avertissement.
nu_lock=floor(Nsurf/2);
if nargin>=3 && ~isempty(numax) && numax>=nu_lock
    numax=floor(numax);
else
    if nargin>=3 && ~isempty(numax) && numax>0 && numax<nu_lock
        warning('cogging_mec:numaxTropBas', ...
            'numax = %d < Nsurf/2 = %d : rang deficient. Valeur imposee.', ...
            numax,nu_lock);
    end
    numax=nu_lock;
end

% ---- grille de surface + classification fer / ouverture ----
dth=2*pi/Nsurf; ths=(0:Nsurf-1)*dth; dths=dth*ones(1,Nsurf);
wtf=taus-M.ws0/Rs;                            % arc de face de dent
tooth=zeros(1,Nsurf); isFe=false(1,Nsurf); nb1=zeros(1,Nsurf); nb2=nb1;
for j=1:Nsurf
    i0=mod(round(ths(j)/taus),Ns);
    dd=angle(exp(1i*(ths(j)-i0*taus)));
    tooth(j)=i0+1;
    if abs(dd)<wtf/2
        isFe(j)=true;
    elseif dd>0
        nb1(j)=i0+1;              nb2(j)=mod(i0+1,Ns)+1;
    else
        nb1(j)=mod(i0-1+Ns,Ns)+1; nb2(j)=i0+1;
    end
end

if nargin<8||isempty(basis), basis='p0'; end
AG=airgap_magnet(M,ths,dths,numax,basis);

% ---- reseau statorique ----
g_face=mu0*muI*(dth*Rs)*L/(M.hs0+M.hs1);
Gt    =mu0*muI*M.wst1*L/M.hs2;
Gy    =mu0*muI*M.wsy*L/(taus*(M.Rso-M.wsy/2));
Gfr   =kfringe*mu0*(dth*Rs)*L/(M.ws0/2);

Ntot=Nsurf+2*Ns; TB=@(i)Nsurf+i; YY=@(i)Nsurf+Ns+i;
A=zeros(Ntot);
    function add(a,b,g)
        A(a,a)=A(a,a)+g; A(b,b)=A(b,b)+g; A(a,b)=A(a,b)-g; A(b,a)=A(b,a)-g;
    end
for j=1:Nsurf
    if isFe(j)
        add(j,TB(tooth(j)),g_face);
    elseif Gfr>0
        add(j,TB(nb1(j)),Gfr); add(j,TB(nb2(j)),Gfr);
    end
end
for i=1:Ns
    add(TB(i),YY(i),Gt); add(YY(i),YY(mod(i,Ns)+1),Gy);
end
A(1:Nsurf,1:Nsurf)=A(1:Nsurf,1:Nsurf)-AG.Y;

ref=YY(1); keep=true(Ntot,1); keep(ref)=false;

% ---- balayage rotor ----
phis=linspace(0,span,Np);
nu=AG.nu; brm=AG.brm;
Isrc=zeros(Ntot,Np);
Isrc(1:Nsurf,:)=L*Rs*pi*( AG.Wc.'*(brm.*cos(nu*phis)) + AG.Ws.'*(brm.*sin(nu*phis)) );
U=zeros(Ntot,Np);
U(keep,:)=A(keep,keep)\Isrc(keep,:);

% ---- flux de dent et de culasse (pertes fer locales, C2) ----
iT=Nsurf+(1:Ns); iY=Nsurf+Ns+(1:Ns);
R.PhiT = Gt*(U(iT,:)-U(iY,:));                       % Ns x Np
R.PhiY = Gy*(U(iY,:)-U(iY([2:Ns 1]),:));

% ---- couple MST par harmoniques (Parseval, exact) ----
Us=U(1:Nsurf,:); Usc=AG.Wc*Us; Uss=AG.Ws*Us;
cphi=cos(nu*phis); sphi=sin(nu*phis);
Brc=AG.bru.*Usc + AG.brmq.*cphi;   Brs=AG.bru.*Uss + AG.brmq.*sphi;
Btc=-AG.btu.*Uss - AG.btmq.*sphi;  Bts=AG.btu.*Usc + AG.btmq.*cphi;
T=(L*M.rmid^2/mu0)*pi*sum(Brc.*Btc + Brs.*Bts, 1);
T=T-mean(T);
Yf=abs(fft(T)); Yf=Yf(2:max(2,floor(Np/2))); [~,km]=max(Yf);

R.T=T; R.phis=phis; R.Tpp=(max(T)-min(T))*1e3; R.cycles=km;
R.order=round(km*(2*pi/span)); R.numax_eff=numax; R.Nsurf=Nsurf;
R.U=U; R.Gt=Gt; R.Gy=Gy; R.span=span; R.AG=AG; R.Usurf=U(1:Nsurf,:);

% ---- metriques et profils de champ a phi=0 ----
thq=linspace(0,2*pi,3601);
[Br,Bt]=AG.field(U(1:Nsurf,1),0,thq);
R.Bmean=mean(abs(Br)); R.Bpeak=max(Br); R.Btrms=sqrt(mean(Bt.^2));
R.Bg1=hypot(2*mean(Br.*cos(M.p*thq)),2*mean(Br.*sin(M.p*thq)));
R.thq=thq; R.Br=Br; R.Bt=Bt;

% ---- controle croise : MST par INTEGRATION DIRECTE (21 positions) ----
nc=min(21,Np); ic=unique(round(linspace(1,Np,nc)));
thd=linspace(0,2*pi,4201); thd(end)=[];
Tdir=zeros(1,numel(ic));
for q=1:numel(ic)
    [Brq,Btq]=AG.field(U(1:Nsurf,ic(q)),phis(ic(q)),thd);
    Tdir(q)=(L*M.rmid^2/mu0)*trapz([thd 2*pi],[Brq.*Btq Brq(1)*Btq(1)]);
end
R.Tdir=Tdir-mean(Tdir); R.Tpar=T(ic)-mean(T(ic)); R.ic=ic;
if max(abs(R.Tpar))>0
    R.err_dir=max(abs(R.Tdir-R.Tpar))/max(abs(R.Tpar))*100;
else, R.err_dir=NaN; end
end
