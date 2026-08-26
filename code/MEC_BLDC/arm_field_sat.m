function A = arm_field_sat(M, Nsurf, kfringe, Np)
%ARM_FIELD_SAT  Champ d'induit MODULE PAR LA SATURATION TOURNANTE.
%
%   POURQUOI. Le potentiel de surface cree par 1 A de phase etait jusqu'ici
%   calcule UNE FOIS, avec une permeabilite de fer uniforme et fixe. C'est
%   exact pour un reseau lineaire -- mais les dents sont polarisees a ~1.5 T
%   par les aimants, et CETTE POLARISATION TOURNE AVEC LE ROTOR. La
%   permeabilite offerte au champ d'induit varie donc dent par dent au fil
%   de la rotation, et le champ d'induit vu par l'aimant acquiert, dans le
%   repere rotor, des harmoniques asynchrones supplementaires. C'est un
%   mecanisme connu et souvent dominant pour les pertes aimant des machines
%   a bobinage concentre saturees : le produit "champ d'induit x saturation
%   tournante" n'existe dans aucun des deux problemes pris separement.
%
%   CE QUI EST CALCULE. Pour chaque position rotor phi :
%     1. resolution NON LINEAIRE du reseau avec les seuls aimants (courbe
%        B(H) M350-50A reelle) -> permeabilite mu_dent(i,phi) et
%        mu_culasse(i,phi) ;
%     2. a ces permeabilites FIGEES, resolution du probleme d'INDUIT pour
%        1 A dans chaque phase -> potentiels de surface U_k(theta_s ; phi).
%   La polarisation est prise a vide : a 2 A le flux de dent est domine par
%   les aimants (1.5 T contre ~0.1 T d'induit), l'approximation est faite
%   et assumee.
%
%   COUT. Les noeuds de surface portent le DtN, dense mais INVARIANT : on
%   les condense une fois (complement de Schur), apres quoi chaque position
%   ne coute qu'un systeme 30x30 par iteration, plus trois remontees de
%   surface. Le tout tient en une minute.

mu0=4*pi*1e-7; Rs=M.Rsi; L=M.ls; Ns=M.Ns; taus=2*pi/Ns; Ntc=M.Ntc;
PH={[1 -2 -15 3 14],[6 -7 -5 8 4],[11 -12 -10 13 9]};
BH=bh_curve();
if nargin<2||isempty(Nsurf),   Nsurf=1260; end
if nargin<3||isempty(kfringe), kfringe=0.325; end
if nargin<4||isempty(Np),      Np=361;  end

%% ---------- topologie ---------------------------------------------------
numax=floor(Nsurf/2);
dth=2*pi/Nsurf; ths=(0:Nsurf-1)*dth; dths=dth*ones(1,Nsurf);
wtf=taus-M.ws0/Rs;
tooth=zeros(1,Nsurf); isFe=false(1,Nsurf); nb1=zeros(1,Nsurf); nb2=nb1;
for j=1:Nsurf
    i0=mod(round(ths(j)/taus),Ns); dd=angle(exp(1i*(ths(j)-i0*taus)));
    tooth(j)=i0+1;
    if abs(dd)<wtf/2, isFe(j)=true;
    elseif dd>0, nb1(j)=i0+1; nb2(j)=mod(i0+1,Ns)+1;
    else,        nb1(j)=mod(i0-1+Ns,Ns)+1; nb2(j)=i0+1; end
end
AG=airgap_magnet(M,ths,dths,numax); nu=AG.nu;

%% ---------- blocs invariants -------------------------------------------
gf =mu0*M.muI*(dth*Rs)*L/(M.hs0+M.hs1);
Gt0=mu0*M.wst1*L/M.hs2;  Gy0=mu0*M.wsy*L/(taus*(M.Rso-M.wsy/2));
Gfr=kfringe*mu0*(dth*Rs)*L/(M.ws0/2);
gs=zeros(Nsurf,1); Ast=zeros(Nsurf,Ns); dgt=zeros(Ns,1);
for j=1:Nsurf
    if isFe(j)
        gs(j)=gs(j)+gf; Ast(j,tooth(j))=Ast(j,tooth(j))-gf;
        dgt(tooth(j))=dgt(tooth(j))+gf;
    elseif Gfr>0
        gs(j)=gs(j)+2*Gfr;
        Ast(j,nb1(j))=Ast(j,nb1(j))-Gfr; Ast(j,nb2(j))=Ast(j,nb2(j))-Gfr;
        dgt(nb1(j))=dgt(nb1(j))+Gfr;     dgt(nb2(j))=dgt(nb2(j))+Gfr;
    end
end
Ass=diag(gs)-AG.Y; dAss=decomposition(Ass,'lu');
K=Ast.'*(dAss\Ast); Kdg=-K+diag(dgt);

%% ---------- source des aimants, condensee ------------------------------
phis=linspace(0,2*pi,Np);                    % TOUR COMPLET (gcd(Ns,Nm)=1)
Isr=L*Rs*pi*( AG.Wc.'*(AG.brm.*cos(nu*phis)) + AG.Ws.'*(AG.brm.*sin(nu*phis)) );
Cpm=Ast.'*(dAss\Isr);

%% ---------- motifs de FMM de branche par phase --------------------------
Fp=zeros(Ns,3);
for k=1:3, for c=PH{k}, Fp(abs(c),k)=Fp(abs(c),k)+sign(c)*Ntc; end, end

%% ---------- balayage rotor ---------------------------------------------
kp=true(2*Ns,1); kp(Ns+1)=false; iY=Ns+(1:Ns);
nx=[2:Ns 1]; pv=[Ns 1:Ns-1];
idg=sub2ind([Ns Ns],(1:Ns).',(1:Ns).');
iup=sub2ind([Ns Ns],(1:Ns).',nx.'); idn=sub2ind([Ns Ns],nx.',(1:Ns).');
Usurf=zeros(Nsurf,3,Np); MUT=zeros(Ns,Np);
mut=ones(Ns,1)*M.muI; muy=mut;
for q=1:Np
    % --- (1) etat de saturation cree par les AIMANTS a cette position ---
    for it=1:80
        Gt=mut*Gt0; Gy=0.5*(muy+muy(nx))*Gy0;
        Ayy=zeros(Ns); Ayy(idg)=Gt+Gy+Gy(pv);
        Ayy(iup)=Ayy(iup)-Gy; Ayy(idn)=Ayy(idn)-Gy;
        DG=diag(Gt); A2=[Kdg+DG, -DG; -DG, Ayy];
        r2=[-Cpm(:,q); zeros(Ns,1)];
        U2=zeros(2*Ns,1); U2(kp)=A2(kp,kp)\r2(kp);
        Ut=U2(1:Ns); Uy=U2(iY);
        Phi=Gt.*(Ut-Uy); Phy=Gy.*(Uy-Uy(nx));
        Bt=Phi/(M.wst1*L*M.Ki); By=Phy/(M.wsy*L*M.Ki);
        nt=BH.mur(abs(Bt)); nyk=BH.mur(abs(By));
        er=max(max(abs(nt-mut)./mut),max(abs(nyk-muy)./muy));
        w=0.35; mut=exp((1-w)*log(mut)+w*log(nt));
        muy=exp((1-w)*log(muy)+w*log(nyk));
        if er<1e-6, break; end
    end
    MUT(:,q)=mut;
    % --- (2) probleme d'INDUIT a permeabilites FIGEES -------------------
    Gt=mut*Gt0; Gy=0.5*(muy+muy(nx))*Gy0;
    Ayy=zeros(Ns); Ayy(idg)=Gt+Gy+Gy(pv);
    Ayy(iup)=Ayy(iup)-Gy; Ayy(idn)=Ayy(idn)-Gy;
    DG=diag(Gt); A2=[Kdg+DG, -DG; -DG, Ayy];
    for k=1:3
        r2=[-Gt.*Fp(:,k); Gt.*Fp(:,k)];      % source d'induit seule
        U2=zeros(2*Ns,1); U2(kp)=A2(kp,kp)\r2(kp);
        Usurf(:,k,q)=dAss\(-Ast*U2(1:Ns));   % potentiel de surface
    end
end

A.U=Usurf; A.phis=phis; A.mut=MUT; A.Nsurf=Nsurf; A.kfringe=kfringe;
A.mu_min=min(MUT(:)); A.mu_max=max(MUT(:));
end
