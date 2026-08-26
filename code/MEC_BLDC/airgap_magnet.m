function AG = airgap_magnet(M, ths, dths, numax, basis)
%
%   BASE DE SURFACE (argument 'basis', defaut 'p0') -- tache V1.
%     'p0' : potentiel CONSTANT par colonne, projection en 1/nu. Le terme
%            d'energie decroit alors en 1/nu : serie harmonique, DIVERGENTE.
%     'p1' : chapeau (lineaire par morceaux), projection en 1/nu^2, terme
%            d'energie en 1/nu^3 : ABSOLUMENT convergent.
%   Meme portage que mec.airgap_fourier cote machine asynchrone. Voir
%   RUN_V1_PMSM_BASE.
%AIRGAP_MAGNET  Operateur DtN a 2 couches (aimant + entrefer) couple au stator.
%
%   AG = airgap_magnet(M, ths, dths, numax) construit le couplage EXACT de la
%   couronne aimant (rmi<=r<=Rro) + entrefer (Rro<=r<=Rsi) vu depuis les
%   NOEUDS DE SURFACE STATOR (arcs ths, largeurs dths, a r=Rsi), la culasse
%   rotor a la masse (phi=0 en rmi). Generalise magnet_subdomain.m : le
%   potentiel de surface stator est porte par le reseau de reluctances
%   statorique -> DENTURE, Bt tangentiel et couple de detente EMERGENT.
%
%   Formulation par harmonique nu = 1..numax (TOUS les ordres entiers : la
%   denture excite des ordres non multiples de p), base sinh/cosh LOCALE a
%   chaque couche (bien conditionnee, pas de r^nu qui deborde). Une inconnue
%   par harmonique : le potentiel d'interface phi_i(nu). Reponse LINEAIRE :
%       Br_bore(nu) = g(nu)*Us_nu + brm(nu)
%   La source d'aimant brm(nu) n'est non nulle qu'aux nu = (impair)*p.
%
%   Sorties : AG.Y (Ms x Ms, flux_entrant_fer = Y*Us + Isrc), AG.Isrc(phi),
%   AG.field(Us,phi,thq) -> [Br,Bt] au mi-entrefer, AG.g, AG.brm, AG.nu.

mu0=4*pi*1e-7;
Rs=M.Rsi; Rro=M.Rro; rmi=M.rmi; mur=M.mu_r; p=M.p; L=M.ls; rmid=M.rmid;
ths=ths(:).'; dths=dths(:).'; Ms=numel(ths);
nu=(1:numax).';

Ua=log(Rs/Rro); Um=log(Rro/rmi); umid=log(rmid/Rro);
Sa=sinh(nu*Ua); Ca=cosh(nu*Ua); coth_a=Ca./Sa;
Sm=sinh(nu*Um); Cm=cosh(nu*Um); coth_m=Cm./Sm;
den=(nu/Rro).*(coth_a + mur*coth_m);

% --- magnetisation : Mn non nul aux nu = (impair)*p ---
Mn=zeros(size(nu));
isM = (mod(nu,p)==0) & (mod(nu/p,2)==1);
nnM = nu(isM)/p;                                   % n impair
Mn(isM)=(4*M.Br./(mu0*nnM*pi)).*sin(nnM*pi*M.embrace/2);
Pn=zeros(size(nu)); Pn(isM)=Mn(isM)./(mur*(1-nu(isM).^2));

% --- reponse du potentiel d'interface phi_i = alphau*Us + betasrc ---
alphau = ((nu/Rro)./Sa)./den;
betasrc= ( mur*(nu/Rro).*Pn.*(Rro*coth_m - rmi./Sm) - mur*Pn + Mn )./den;

% --- admittance de surface et source (Br au bore x=1) ---
g   = -mu0*(nu/Rs)./Sa.*(Ca - alphau);
brm =  mu0*(nu/Rs).*betasrc./Sa;

% --- facteurs de champ au mi-entrefer ---
cu=cosh(nu*umid); su=sinh(nu*umid);
ca=cosh(nu*(Ua-umid)); sa2=sinh(nu*(Ua-umid));
bru = -mu0*(nu/rmid).*(cu - alphau.*ca)./Sa;       % Br / Us_nu (cos)
btu =  mu0*(nu/rmid).*(su + alphau.*sa2)./Sa;      % Bt / Us_nu
brmq=  mu0*(nu/rmid).*(betasrc.*ca)./Sa;           % Br aimant au mi-entrefer
btmq=  mu0*(nu/rmid).*(betasrc.*sa2)./Sa;          % Bt aimant au mi-entrefer

% --- projections de Fourier des segments de surface ---
if nargin<5||isempty(basis), basis='p0'; end
switch lower(basis)
case 'p0', prj=(2./(nu*pi)).*sin(nu*dths/2);            % creneau
case 'p1', prj=(4./(pi*(nu.^2).*dths)).*sin(nu.*dths/2).^2;  % chapeau
otherwise, error('airgap_magnet:basis','base inconnue : %s',basis);
end
Wc=prj.*cos(nu*ths);      % (numax x Ms)
Ws=prj.*sin(nu*ths);

Y = L*Rs*pi*( Wc.'*(g.*Wc) + Ws.'*(g.*Ws) );       % symetrique

AG.Y=Y; AG.g=g; AG.brm=brm; AG.nu=nu; AG.Wc=Wc; AG.Ws=Ws; AG.basis=lower(basis);
AG.Ms=Ms; AG.L=L; AG.Rs=Rs; AG.mu0=mu0; AG.rmid=rmid;
AG.bru=bru; AG.btu=btu; AG.brmq=brmq; AG.btmq=btmq;
% --- grandeurs necessaires au champ DANS la couronne d'aimant -----------
%  Potentiel d'interface : phi_i = alphau*Us_nu + betasrc  (par harmonique).
%  Dans l'aimant (rmi<=r<=Rro) la partie pilotee par Us verifie Laplace avec
%  phi=0 en rmi (culasse rotor) et phi=phi_i en Rro :
%      phi_m(r) = phi_i * sinh(nu*v)/sinh(nu*Um),   v = log(r/rmi)
%  d'ou  B_r = -mu0*mur*(nu/r)*phi_i*cosh(nu*v)/sinh(nu*Um).
%  (La part de MAGNETISATION est STATIQUE dans le repere rotor : elle ne
%   produit aucun dB/dt, donc aucune perte -> seule la part Us compte.)
AG.alphau=alphau; AG.betasrc=betasrc; AG.Um=Um; AG.rmi=rmi; AG.Rro=Rro; AG.mur=mur;
%  Harmoniques de MAGNETISATION exposes : necessaires au champ COMPLET dans
%  la couche d'aimant (B_r = -mu0*mur*dphi/dr + mu0*M_r), donc au facteur de
%  fuite k_l mesure DANS l'aimant comme le fait ANSYS (phi_Bmk).
AG.Mn=Mn; AG.Pn=Pn;
AG.Brm_of=@(r,phi_ic,phi_is,thq) magnet_field(nu,mu0,mur,rmi,Um,r,phi_ic,phi_is,thq);

AG.Isrc = @(phi) L*Rs*pi*( Wc.'*(brm.*cos(nu*phi)) + Ws.'*(brm.*sin(nu*phi)) );

    function [Br,Bt]=field(Us,phi,thq)
        Us=Us(:); thq=thq(:).';
        Usc=Wc*Us; Uss=Ws*Us;
        Cq=cos(nu*thq); Sq=sin(nu*thq);
        Cqp=cos(nu*(thq-phi)); Sqp=sin(nu*(thq-phi));
        % Br : surface (cos/sin) + aimant (cos, tourne avec le rotor)
        Br = (bru.*Usc).'*Cq + (bru.*Uss).'*Sq + brmq.'*Cqp;
        % Bt : surface (cos->+sin, sin->-cos) + aimant (+sin)
        Bt = (btu.*Usc).'*Sq - (btu.*Uss).'*Cq + btmq.'*Sqp;
        Br=Br(:).'; Bt=Bt(:).';
    end
AG.field=@field;
end

% ======================================================================
function Br = magnet_field(nu,mu0,mur,rmi,Um,r,phi_ic,phi_is,thq)
%  Composante radiale DANS la couronne d'aimant, au rayon r, pour un
%  potentiel d'interface de composantes harmoniques (phi_ic, phi_is).
    thq=thq(:).'; v=log(r/rmi);
    kk = -mu0*mur*(nu/r).*cosh(nu*v)./sinh(nu*Um);
    Br = (kk.*phi_ic).'*cos(nu*thq) + (kk.*phi_is).'*sin(nu*thq);
    Br = Br(:).';
end
