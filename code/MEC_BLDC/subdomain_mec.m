function R = subdomain_mec(M, K, numax, muI, phis)
%SUBDOMAIN_MEC  Modele HYBRIDE sous-domaines analytiques <-> reseau de reluctances.
%
%   R = subdomain_mec(M, K, numax, muI, phis)
%     K     : nombre de modes propres par ouverture d'encoche (6-10 suffit)
%     numax : ordre max de la serie de Fourier d'entrefer (>= 3*K*pi/delta)
%     muI   : permeabilite relative du fer (modele lineaire ; [] -> M.muI)
%     phis  : positions rotor [rad mec] (vecteur)
%
%   *** SUPPRESSION DU DERNIER PARAMETRE AJUSTE (kfringe) ***
%   Dans cogging_mec.m l'ouverture d'encoche etait reduite a une conductance
%   FORFAITAIRE Gfr = kfringe*mu0*(dth*Rs)*L/(ws0/2), avec kfringe IDENTIFIE
%   sur les bandes de denture FEA. C'etait le seul degre de liberte ajuste du
%   modele -- et precisement celui qui gouverne le spectre de denture (a8, a22)
%   et donc les pertes aimant. Ici l'ouverture est resolue EXACTEMENT :
%
%   SOUS-DOMAINE D'OUVERTURE (par encoche i, r in [Rs,Rt], x = theta-theta_a
%   in [0,delta], Rt = Rs+hs0, delta = ws0/Rs) :
%     - flancs de dents = surfaces de fer EQUIPOTENTIELLES (phi_F,i, phi_F,i+1) ;
%     - haut de l'ouverture : dphi/dr = 0 (le flux ne penetre pas le corps
%       d'encoche a vide) ;
%     - Laplace resolu par separation :
%         phi = phi_L + (phi_R-phi_L)*x/delta + SUM_k a_k*g_k(r)*sin(k*pi*x/delta)
%         g_k(r) = cosh(lam_k*ln(r/Rt))/cosh(lam_k*h),  h = ln(Rt/Rs), lam_k=k*pi/delta
%       (la partie lineaire verifie Laplace ET les deux flancs ; les modes
%        propres s'annulent sur les flancs -> base complete).
%     - au bore :  Br = (mu0/Rs)*SUM_k a_k*lam_k*tanh(lam_k*h)*sin(lam_k*x)
%     - flux sortant par les FLANCS (= fuite de bec, exacte) :
%         Phi_gauche = mu0*L*[h*(phi_R-phi_L)/delta + SUM_k a_k*tanh(lam_k*h)]
%         Phi_droite = -mu0*L*[h*(phi_R-phi_L)/delta + SUM_k a_k*(-1)^k*tanh(lam_k*h)]
%       (conservation verifiee : Phi_bas = Phi_gauche + Phi_droite)
%
%   RACCORD EN MODES avec la couronne aimant+entrefer (airgap_magnet) :
%   le potentiel de bore est parametre par [phi_F ; a] ; ses coefficients de
%   Fourier sont obtenus par projections ANALYTIQUES (aucune discretisation de
%   surface, donc aucun repliement) ; Br au bore vient du DtN ; on impose
%   l'egalite des projections modales sur chaque ouverture.
%
%   Inconnues : phi_F (Ns) + a (Ns*K) + phi_TB (Ns) + phi_YY (Ns).
%
%   Sorties : R.Br/R.Bt/R.thq (mi-entrefer), R.a8/a22/a23, R.Bg1/Bmean/Bpeak/
%   R.Btrms, R.T (couple MST), R.PhiT/R.PhiY (flux dent/culasse), R.U, R.AG,
%   R.Usc/R.Uss (harmoniques du potentiel de bore), R.balance (controle).

mu0=4*pi*1e-7;
Rs=M.Rsi; L=M.ls; Ns=M.Ns; taus=2*pi/Ns;
if nargin<4 || isempty(muI), muI=M.muI; end
if nargin<5 || isempty(phis), phis=0; end
phis=phis(:).';  Np=numel(phis);

delta = M.ws0/Rs;              % ouverture angulaire
wtf   = taus - delta;          % arc de face de dent
Rt    = Rs + M.hs0;            % haut de l'ouverture
h     = log(Rt/Rs);
lam   = (1:K)'*pi/delta;       % nombres d'onde des modes d'ouverture
th_t  = (0:Ns-1)*taus;                 % centres de dents
th_o  = th_t + wtf/2 + delta/2;        % centres d'ouvertures

if isempty(numax), numax = ceil(3*lam(end)); end
nu=(1:numax)';

AG = airgap_magnet(M, 0, 2*pi, numax);   % ths/dths non utilises ici (on
                                          % fournit nos propres projections)
g=AG.g; brm=AG.brm;

% ================= projections analytiques =================
% convention : Usc(nu) = (1/pi)*Int phi(th) cos(nu th) dth   (idem sin)
segC=@(w,thc) (2./(nu*pi)).*sin(nu*w/2).*cos(nu*thc);   % segment constant=1
segS=@(w,thc) (2./(nu*pi)).*sin(nu*w/2).*sin(nu*thc);

% --- faces de dents ---
Fc=zeros(numax,Ns); Fs=zeros(numax,Ns);
for i=1:Ns, Fc(:,i)=segC(wtf,th_t(i)); Fs(:,i)=segS(wtf,th_t(i)); end

% --- ouvertures : partie constante et rampe ---
Oc=zeros(numax,Ns); Os=zeros(numax,Ns);      % constante = 1 sur l'ouverture
Qc=zeros(numax,Ns); Qs=zeros(numax,Ns);      % partie IMPAIRE s/delta
Qk=(2/delta)*( sin(nu*delta/2)./nu.^2 - (delta/2)*cos(nu*delta/2)./nu );
for i=1:Ns
    Oc(:,i)=segC(delta,th_o(i)); Os(:,i)=segS(delta,th_o(i));
    Qc(:,i)=-(1/pi)*sin(nu*th_o(i)).*Qk;
    Qs(:,i)= (1/pi)*cos(nu*th_o(i)).*Qk;
end
% x/delta = 1/2 + s/delta   ;   1-x/delta = 1/2 - s/delta
RampC = 0.5*Oc + Qc;   RampS = 0.5*Os + Qs;      % coefficient de phi_F,i+1
CompC = 0.5*Oc - Qc;   CompS = 0.5*Os - Qs;      % coefficient de phi_F,i

% --- modes propres sin(lam_k x) sur chaque ouverture ---
Ac=zeros(numax,Ns*K); As=zeros(numax,Ns*K);
for k=1:K
    lk=lam(k);
    d1=lk-nu; d2=lk+nu;
    s1=sin((lk-nu)*delta/2); s2=sin((lk+nu)*delta/2);
    t1=s1./d1; t1(abs(d1)<1e-12)=delta/2;        % limite lam_k = nu
    t2=s2./d2;
    Ick=t1+t2;    % Int cos(lam s)cos(nu s) ds
    Isk=t1-t2;    % Int sin(lam s)sin(nu s) ds
    for i=1:Ns
        col=(i-1)*K+k;
        if mod(k,2)==1
            sg=(-1)^((k-1)/2);
            Ac(:,col)= (sg/pi)*cos(nu*th_o(i)).*Ick;
            As(:,col)= (sg/pi)*sin(nu*th_o(i)).*Ick;
        else
            sg=(-1)^(k/2);
            Ac(:,col)=-(sg/pi)*sin(nu*th_o(i)).*Isk;
            As(:,col)= (sg/pi)*cos(nu*th_o(i)).*Isk;
        end
    end
end

% --- matrices de projection du potentiel de bore : Usc = Pc*Xs ---
nF=Ns; nA=Ns*K; ns=nF+nA;
Pc=zeros(numax,ns); Ps=zeros(numax,ns);
for i=1:Ns
    ip=mod(i,Ns)+1;                      % tooth i+1
    Pc(:,i)  = Pc(:,i)  + Fc(:,i) + CompC(:,i);      % face i + (1-x/d) sur ouv. i
    Ps(:,i)  = Ps(:,i)  + Fs(:,i) + CompS(:,i);
    Pc(:,ip) = Pc(:,ip) + RampC(:,i);                % (x/d) sur ouverture i
    Ps(:,ip) = Ps(:,ip) + RampS(:,i);
end
Pc(:,nF+1:end)=Ac;  Ps(:,nF+1:end)=As;

% ================= assemblage du systeme =================
%  Sections de FER : le foisonnement ki reduit la section utile (et donc
%  augmente l'induction reelle du fer de 1/ki) -- omis auparavant.
ki=M.Ki;
Atip=0.5*(wtf*Rs + M.wst1)*L*ki;    % tete de dent (entonnoir face -> corps)
Abody=M.wst1*L*ki;                  % corps de dent
Ayoke=M.wsy*L*ki;                   % section de culasse
ltip=M.hs0+M.hs1; lbody=M.hs2; lyoke=taus*(M.Rso-M.wsy/2);
NLIN = isstruct(muI) || (ischar(muI)&&strcmpi(muI,'nl'));   % mode non lineaire
if NLIN, BHc=bh_curve('M350'); muv=3000*ones(3*Ns,1);       % [tip;body;yoke]
else,    muv=muI*ones(3*Ns,1); end
g_tip = mu0*muv(1:Ns)      *Atip /ltip;
Gt    = mu0*muv(Ns+(1:Ns)) *Abody/lbody;
Gy    = mu0*muv(2*Ns+(1:Ns))*Ayoke/lyoke;

iF =@(i) i;
iA =@(i,k) nF+(i-1)*K+k;
iTB=@(i) ns+i;
iYY=@(i) ns+Ns+i;
Ntot=ns+2*Ns;

% couplage d'entrefer : flux entrant dans le fer = L*Rs*pi*( W'*(g.*Us + brm) )
GPc = L*Rs*pi*( Fc.'*(g.*Pc) + Fs.'*(g.*Ps) );      % (Ns x ns)  faces
GPa = pi*( Ac.'*(g.*Pc) + As.'*(g.*Ps) );           % (Ns*K x ns) modes

A0=zeros(Ntot);       % partie INDEPENDANTE de mu (entrefer + ouvertures)
Amat=A0;
% --- (1) bilan de flux de face de dent ---
for i=1:Ns
    Amat(iF(i),1:ns) = Amat(iF(i),1:ns) + GPc(i,:);
    % flancs : ouverture i (flanc gauche -> dent i) et ouverture i-1 (flanc droit)
    im=mod(i-2,Ns)+1;  ip=mod(i,Ns)+1;
    cL = mu0*L*h/delta;                    % ouverture i : +h*(phiF,ip-phiF,i)/delta
    Amat(iF(i),iF(ip)) = Amat(iF(i),iF(ip)) + cL;
    Amat(iF(i),iF(i))  = Amat(iF(i),iF(i))  - cL;
    for k=1:K
        Amat(iF(i),iA(i,k)) = Amat(iF(i),iA(i,k)) + mu0*L*tanh(lam(k)*h);
    end
    % ouverture i-1 : flanc droit -> dent i  (=-mu0*L*[h*(phiF,i-phiF,im)/delta + sum a(-1)^k tanh])
    Amat(iF(i),iF(i))  = Amat(iF(i),iF(i))  - cL;
    Amat(iF(i),iF(im)) = Amat(iF(i),iF(im)) + cL;
    for k=1:K
        Amat(iF(i),iA(im,k)) = Amat(iF(i),iA(im,k)) - mu0*L*((-1)^k)*tanh(lam(k)*h);
    end
end
% --- (2) raccord modal des ouvertures ---
kap=zeros(Ns*K,1);
for i=1:Ns, for k=1:K
    col=(i-1)*K+k;
    kap(col)=(mu0/Rs)*lam(k)*tanh(lam(k)*h)*(delta/2);
end, end
for c=1:Ns*K
    Amat(nF+c,1:ns) = Amat(nF+c,1:ns) + GPa(c,:);
    Amat(nF+c,nF+c) = Amat(nF+c,nF+c) - kap(c);
end
A0=Amat;                                   % partie independante de mu

% ================= resolution (non lineaire par position) =================
ref=iYY(1); keep=true(Ntot,1); keep(ref)=false;
X=zeros(Ntot,Np); R.iters=zeros(1,Np); R.muT=zeros(Ns,Np);
    function Am=addiron(gt,Gtv,Gyv)
        Am=A0;
        for i2=1:Ns
            ip2=mod(i2,Ns)+1; im2=mod(i2-2,Ns)+1;
            Am(iF(i2),iF(i2))   = Am(iF(i2),iF(i2))   - gt(i2);
            Am(iF(i2),iTB(i2))  = Am(iF(i2),iTB(i2))  + gt(i2);
            Am(iTB(i2),iF(i2))  = Am(iTB(i2),iF(i2))  + gt(i2);
            Am(iTB(i2),iTB(i2)) = Am(iTB(i2),iTB(i2)) - gt(i2) - Gtv(i2);
            Am(iTB(i2),iYY(i2)) = Am(iTB(i2),iYY(i2)) + Gtv(i2);
            Am(iYY(i2),iTB(i2)) = Am(iYY(i2),iTB(i2)) + Gtv(i2);
            Am(iYY(i2),iYY(i2)) = Am(iYY(i2),iYY(i2)) - Gtv(i2) - Gyv(i2) - Gyv(im2);
            Am(iYY(i2),iYY(ip2))= Am(iYY(i2),iYY(ip2))+ Gyv(i2);
            Am(iYY(i2),iYY(im2))= Am(iYY(i2),iYY(im2))+ Gyv(im2);
        end
    end
    function x=solveq(Am,rhs)
        Ak=Am(keep,keep);
        rsc=1./max(abs(Ak),[],2); rsc(~isfinite(rsc)|rsc==0)=1;
        Ak=diag(rsc)*Ak;
        csc=1./max(abs(Ak),[],1).'; csc(~isfinite(csc)|csc==0)=1;
        Ak=Ak*diag(csc);
        x=zeros(Ntot,1); x(keep)=csc.*(Ak\(rsc.*rhs(keep)));
        R.rcond=rcond(Ak);
    end
for q=1:Np
    ph=phis(q);
    bc=brm.*cos(nu*ph); bs=brm.*sin(nu*ph);
    rhs=zeros(Ntot,1);
    rhs(1:Ns)       = -L*Rs*pi*( Fc.'*bc + Fs.'*bs );    % faces
    rhs(nF+(1:Ns*K))= -pi*( Ac.'*bc + As.'*bs );         % modes
    gt=g_tip; Gtv=Gt; Gyv=Gy;
    if ~NLIN
        X(:,q)=solveq(addiron(gt,Gtv,Gyv),rhs); R.iters(q)=1;
    else
        % --- point fixe amorti sur mu(B) local, PAR DENT et PAR SEGMENT ---
        for it=1:60
            xq=solveq(addiron(gt,Gtv,Gyv),rhs);
            PhiTq=Gtv.*(xq(ns+(1:Ns))-xq(ns+Ns+(1:Ns)));
            PhiYq=Gyv.*(xq(ns+Ns+(1:Ns))-xq(ns+Ns+[2:Ns 1]));
            Phifq=gt.*(xq(1:Ns)-xq(ns+(1:Ns)));
            muT=BHc.mur(abs(PhiTq)/Abody);
            muY=BHc.mur(abs(PhiYq)/Ayoke);
            muP=BHc.mur(abs(Phifq)/Atip);
            gtn=mu0*muP*Atip/ltip; Gtn=mu0*muT*Abody/lbody; Gyn=mu0*muY*Ayoke/lyoke;
            rel=0.5; dmax=max(abs([gtn./gt;Gtn./Gtv;Gyn./Gyv]-1));
            gt=gt+rel*(gtn-gt); Gtv=Gtv+rel*(Gtn-Gtv); Gyv=Gyv+rel*(Gyn-Gyv);
            if dmax<1e-6, break; end
        end
        X(:,q)=xq; R.iters(q)=it; R.muT(:,q)=muT;
    end
    R.gtip=gt; R.Gtv=Gtv; R.Gyv=Gyv;
end
Gt=R.Gtv; Gy=R.Gyv; g_tip=R.gtip;

% ================= post-traitement =================
Xs=X(1:ns,:);
Usc=Pc*Xs; Uss=Ps*Xs;                          % harmoniques du potentiel de bore
cphi=cos(nu*phis); sphi=sin(nu*phis);
Brc=AG.bru.*Usc + AG.brmq.*cphi;  Brs=AG.bru.*Uss + AG.brmq.*sphi;
Btc=-AG.btu.*Uss - AG.btmq.*sphi; Bts=AG.btu.*Usc + AG.btmq.*cphi;

R.T=(L*M.rmid^2/mu0)*pi*sum(Brc.*Btc + Brs.*Bts,1);
R.phis=phis; R.X=X; R.Usc=Usc; R.Uss=Uss; R.AG=AG; R.nu=nu;
R.Brc=Brc; R.Brs=Brs; R.Btc=Btc; R.Bts=Bts;

% flux de dent / de culasse (pertes fer locales)
R.PhiT = Gt.*( X(ns+(1:Ns),:) - X(ns+Ns+(1:Ns),:) );
R.PhiY = Gy.*( X(ns+Ns+(1:Ns),:) - X(ns+Ns+[2:Ns 1],:) );
R.Abody=Abody; R.Ayoke=Ayoke; R.Atip=Atip; R.NLIN=NLIN;
R.Gt=Gt; R.Gy=Gy; R.g_tip=g_tip; R.K=K; R.numax=numax; R.delta=delta;

% profils au mi-entrefer (position 1)
thq=linspace(0,2*pi,3601);
Br=(Brc(:,1).'*cos(nu*thq) + Brs(:,1).'*sin(nu*thq));
Bt=(Btc(:,1).'*cos(nu*thq) + Bts(:,1).'*sin(nu*thq));
R.thq=thq; R.Br=Br; R.Bt=Bt;
R.Bmean=mean(abs(Br)); R.Bpeak=max(Br); R.Btrms=sqrt(mean(Bt.^2));
R.Bg1=hypot(Brc(M.p,1),Brs(M.p,1));
R.a8 =hypot(Brc(8,1), Brs(8,1));
R.a22=hypot(Brc(22,1),Brs(22,1));
R.a23=hypot(Brc(23,1),Brs(23,1));

% controle : conservation du flux dans chaque ouverture
bal=zeros(Ns,1);
for i=1:Ns
    a_i=X(nF+(i-1)*K+(1:K),1);
    ip=mod(i,Ns)+1;
    dphi=X(iF(ip),1)-X(iF(i),1);
    PhiL= mu0*L*( h*dphi/delta + sum(a_i.*tanh(lam*h)) );
    PhiR=-mu0*L*( h*dphi/delta + sum(a_i.*((-1).^(1:K)').*tanh(lam*h)) );
    Phib= mu0*L*sum( a_i.*tanh(lam*h).*(1-(-1).^(1:K)') );
    bal(i)=abs(PhiL+PhiR-Phib)/max(abs(Phib),1e-20);
end
R.balance=max(bal);
end
