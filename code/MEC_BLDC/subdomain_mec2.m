function R = subdomain_mec2(M, K, numax, Nf, muI, phis)
%SUBDOMAIN_MEC2  Sous-domaine d'ouverture + FACE DE DENT SEGMENTEE (bec distribue)
%
%   R = subdomain_mec2(M, K, numax, Nf, muI, phis)
%     Nf : nombre de segments par face de dent (1 = modele equipotentiel de
%          subdomain_mec ; 5-9 pour resoudre le bec)
%     muI: permeabilite relative, ou 'nl' pour mu(B) local NON LINEAIRE
%
%   OBJET. Dans subdomain_mec.m la face de dent est UNE equipotentielle : le
%   fer y est suppose infiniment permeable TANGENTIELLEMENT. Ici la face est
%   scindee en Nf segments relies au corps de dent par le CHEMIN REEL dans le
%   bec :
%       - segment au-dessus du corps (|d| <= wst1/2)  : trajet RADIAL
%         (longueur hs0+hs1, section = arc_segment x L x ki) ;
%       - segment en PORTE-A-FAUX (|d| > wst1/2)      : radial PUIS TANGENTIEL
%         (longueur |d|-wst1/2, section = (hs0+hs1) x L x ki), en serie.
%   Chaque branche porte sa propre mu(B) : les cornes de bec, plus chargees,
%   peuvent saturer et faire DECROCHER le potentiel de la face -> frontiere
%   fer/air GRADUEE au lieu d'une marche nette. Les flancs du sous-domaine
%   d'ouverture sont alors portes par les segments EXTREMES (les cornes), et
%   non plus par un potentiel de dent unique.
%
%   Hypothese testee : cette graduation suffit-elle a expliquer le deficit de
%   -28 % sur la premiere bande de denture ? (cf. README, section diagnostic)
%
%   Sorties : identiques a subdomain_mec + R.phiF (Ns x Nf potentiels de face),
%   R.muShoe (permeabilites des branches de bec).

mu0=4*pi*1e-7;
Rs=M.Rsi; L=M.ls; Ns=M.Ns; taus=2*pi/Ns; ki=M.Ki;
if nargin<4||isempty(Nf), Nf=5; end
if nargin<5||isempty(muI), muI=M.muI; end
if nargin<6||isempty(phis), phis=0; end
phis=phis(:).'; Np=numel(phis);

delta=M.ws0/Rs; wtf=taus-delta; Rt=Rs+M.hs0; h=log(Rt/Rs);
lam=(1:K)'*pi/delta;
th_t=(0:Ns-1)*taus; th_o=th_t+wtf/2+delta/2;
if isempty(numax), numax=ceil(3*lam(end)); end
nu=(1:numax)';
AG=airgap_magnet(M,0,2*pi,numax); g=AG.g; brm=AG.brm;

% ---- decoupage de la face en Nf segments (arcs egaux) ----
wseg=wtf/Nf;                                   % largeur angulaire d'un segment
dseg=(-(Nf-1)/2:(Nf-1)/2)*wseg;                % decalage angulaire / axe de dent
th_f=th_t(:)+dseg;                             % (Ns x Nf) centres de segments
dlin=abs(dseg)*Rs;                             % distance lineaire a l'axe (m)

% ---- projections analytiques ----
segC=@(w,thc)(2./(nu*pi)).*sin(nu*w/2).*cos(nu*thc);
segS=@(w,thc)(2./(nu*pi)).*sin(nu*w/2).*sin(nu*thc);
nF=Ns*Nf; nA=Ns*K; ns=nF+nA;
Pc=zeros(numax,ns); Ps=zeros(numax,ns);
Fc=zeros(numax,nF); Fs=zeros(numax,nF);
iFs=@(i,j)(i-1)*Nf+j;
for i=1:Ns, for j=1:Nf
    c=iFs(i,j);
    Fc(:,c)=segC(wseg,th_f(i,j)); Fs(:,c)=segS(wseg,th_f(i,j));
    Pc(:,c)=Pc(:,c)+Fc(:,c);     Ps(:,c)=Ps(:,c)+Fs(:,c);
end, end
% rampe sur chaque ouverture : entre le DERNIER segment de la dent i et le
% PREMIER segment de la dent i+1
Qk=(2/delta)*( sin(nu*delta/2)./nu.^2 - (delta/2)*cos(nu*delta/2)./nu );
for i=1:Ns
    ip=mod(i,Ns)+1;
    Oc=segC(delta,th_o(i)); Os=segS(delta,th_o(i));
    Qc=-(1/pi)*sin(nu*th_o(i)).*Qk; Qs=(1/pi)*cos(nu*th_o(i)).*Qk;
    Pc(:,iFs(i ,Nf))=Pc(:,iFs(i ,Nf))+0.5*Oc-Qc;   % (1-x/delta)
    Ps(:,iFs(i ,Nf))=Ps(:,iFs(i ,Nf))+0.5*Os-Qs;
    Pc(:,iFs(ip,1 ))=Pc(:,iFs(ip,1 ))+0.5*Oc+Qc;   % (x/delta)
    Ps(:,iFs(ip,1 ))=Ps(:,iFs(ip,1 ))+0.5*Os+Qs;
end
% modes propres d'ouverture
Ac=zeros(numax,nA); As=zeros(numax,nA);
for k=1:K
    lk=lam(k); d1=lk-nu; d2=lk+nu;
    t1=sin(d1*delta/2)./d1; t1(abs(d1)<1e-12)=delta/2;
    t2=sin(d2*delta/2)./d2;
    Ick=t1+t2; Isk=t1-t2;
    for i=1:Ns
        c=(i-1)*K+k;
        if mod(k,2)==1
            sg=(-1)^((k-1)/2);
            Ac(:,c)= (sg/pi)*cos(nu*th_o(i)).*Ick;
            As(:,c)= (sg/pi)*sin(nu*th_o(i)).*Ick;
        else
            sg=(-1)^(k/2);
            Ac(:,c)=-(sg/pi)*sin(nu*th_o(i)).*Isk;
            As(:,c)= (sg/pi)*cos(nu*th_o(i)).*Isk;
        end
    end
end
Pc(:,nF+1:end)=Ac; Ps(:,nF+1:end)=As;

% ---- permeances du BEC : trajet radial (+ tangentiel si porte-a-faux) ----
Arad = wseg*Rs*L*ki;                 % section radiale d'un segment
Atan = (M.hs0+M.hs1)*L*ki;           % section tangentielle du bec
lrad = M.hs0+M.hs1;
ltan = max(dlin - M.wst1/2, 0);      % longueur de porte-a-faux (m)
Abody=M.wst1*L*ki; Ayoke=M.wsy*L*ki; lbody=M.hs2; lyoke=taus*(M.Rso-M.wsy/2);
NLIN=ischar(muI)&&strcmpi(muI,'nl');
if NLIN, BHc=bh_curve('M350'); mu0v=3000; else, mu0v=muI; end

% ---- indices ----
iA=@(i,k)nF+(i-1)*K+k; iTB=@(i)ns+i; iYY=@(i)ns+Ns+i; Ntot=ns+2*Ns;
GPc=L*Rs*pi*( Fc.'*(g.*Pc) + Fs.'*(g.*Ps) );      % (nF x ns)
GPa=pi*( Ac.'*(g.*Pc) + As.'*(g.*Ps) );           % (nA x ns)
kap=zeros(nA,1);
for i=1:Ns, for k=1:K, kap((i-1)*K+k)=(mu0/Rs)*lam(k)*tanh(lam(k)*h)*(delta/2); end, end

A0=zeros(Ntot);
for c=1:nF, A0(c,1:ns)=A0(c,1:ns)+GPc(c,:); end
for c=1:nA, A0(nF+c,1:ns)=A0(nF+c,1:ns)+GPa(c,:); A0(nF+c,nF+c)=A0(nF+c,nF+c)-kap(c); end
% flancs d'ouverture -> segments EXTREMES des deux dents adjacentes
for i=1:Ns
    ip=mod(i,Ns)+1; cL=mu0*L*h/delta;
    a1=iFs(i,Nf); a2=iFs(ip,1);
    A0(a1,a2)=A0(a1,a2)+cL;  A0(a1,a1)=A0(a1,a1)-cL;      % flanc gauche -> dent i
    A0(a2,a2)=A0(a2,a2)-cL;  A0(a2,a1)=A0(a2,a1)+cL;      % flanc droit  -> dent i+1
    for k=1:K
        A0(a1,iA(i,k))=A0(a1,iA(i,k))+mu0*L*tanh(lam(k)*h);
        A0(a2,iA(i,k))=A0(a2,iA(i,k))-mu0*L*((-1)^k)*tanh(lam(k)*h);
    end
end

ref=iYY(1); keep=true(Ntot,1); keep(ref)=false;
X=zeros(Ntot,Np); R.muShoe=zeros(Ns*Nf,Np);
    function Am=addiron(gsh,Gtv,Gyv)
        Am=A0;
        for i2=1:Ns
            ip2=mod(i2,Ns)+1; im2=mod(i2-2,Ns)+1;
            for j2=1:Nf
                c2=iFs(i2,j2); gg=gsh(c2);
                Am(c2,c2)=Am(c2,c2)-gg; Am(c2,iTB(i2))=Am(c2,iTB(i2))+gg;
                Am(iTB(i2),c2)=Am(iTB(i2),c2)+gg;
                Am(iTB(i2),iTB(i2))=Am(iTB(i2),iTB(i2))-gg;
            end
            Am(iTB(i2),iTB(i2))=Am(iTB(i2),iTB(i2))-Gtv(i2);
            Am(iTB(i2),iYY(i2))=Am(iTB(i2),iYY(i2))+Gtv(i2);
            Am(iYY(i2),iTB(i2))=Am(iYY(i2),iTB(i2))+Gtv(i2);
            Am(iYY(i2),iYY(i2))=Am(iYY(i2),iYY(i2))-Gtv(i2)-Gyv(i2)-Gyv(im2);
            Am(iYY(i2),iYY(ip2))=Am(iYY(i2),iYY(ip2))+Gyv(i2);
            Am(iYY(i2),iYY(im2))=Am(iYY(i2),iYY(im2))+Gyv(im2);
        end
    end
    function x=solveq(Am,rhs)
        Ak=Am(keep,keep);
        rsc=1./max(abs(Ak),[],2); rsc(~isfinite(rsc)|rsc==0)=1; Ak=diag(rsc)*Ak;
        csc=1./max(abs(Ak),[],1).'; csc(~isfinite(csc)|csc==0)=1; Ak=Ak*diag(csc);
        x=zeros(Ntot,1); x(keep)=csc.*(Ak\(rsc.*rhs(keep))); R.rcond=rcond(Ak);
    end
    function gs=shoeperm(muv)
        % serie : radial (toujours) + tangentiel (porte-a-faux seulement)
        gs=zeros(Ns*Nf,1);
        for i3=1:Ns, for j3=1:Nf
            c3=iFs(i3,j3); m3=muv(c3);
            Rrad=lrad/(mu0*m3*Arad);
            Rtan=ltan(j3)/(mu0*m3*Atan);
            gs(c3)=1/(Rrad+Rtan);
        end, end
    end

for q=1:Np
    ph=phis(q); bc=brm.*cos(nu*ph); bs=brm.*sin(nu*ph);
    rhs=zeros(Ntot,1);
    rhs(1:nF)      = -L*Rs*pi*( Fc.'*bc + Fs.'*bs );
    rhs(nF+(1:nA)) = -pi*( Ac.'*bc + As.'*bs );
    muv=mu0v*ones(Ns*Nf,1); muT=mu0v*ones(Ns,1); muY=mu0v*ones(Ns,1);
    gsh=shoeperm(muv); Gtv=mu0*muT*Abody/lbody; Gyv=mu0*muY*Ayoke/lyoke;
    nit=1;
    if NLIN
        for it=1:80
            xq=solveq(addiron(gsh,Gtv,Gyv),rhs);
            % flux de chaque branche de bec -> induction dans la section la
            % plus CONTRAINTE (tangentielle si porte-a-faux, sinon radiale)
            Phis=zeros(Ns*Nf,1);
            for i4=1:Ns, for j4=1:Nf
                c4=iFs(i4,j4); Phis(c4)=gsh(c4)*(xq(c4)-xq(iTB(i4)));
            end, end
            Asec=repmat( (ltan>0).*Atan + (ltan<=0)*Arad, Ns,1);
            Asec=reshape(repmat(((ltan(:)>0)*Atan+(ltan(:)<=0)*Arad).',Ns,1).',[],1);
            Bsh=abs(Phis)./Asec;
            PhiT=Gtv.*(xq(ns+(1:Ns))-xq(ns+Ns+(1:Ns)));
            PhiY=Gyv.*(xq(ns+Ns+(1:Ns))-xq(ns+Ns+[2:Ns 1]));
            mun=BHc.mur(Bsh); muTn=BHc.mur(abs(PhiT)/Abody); muYn=BHc.mur(abs(PhiY)/Ayoke);
            rel=0.4;
            muv=muv+rel*(mun-muv); muT=muT+rel*(muTn-muT); muY=muY+rel*(muYn-muY);
            gshn=shoeperm(muv); Gtn=mu0*muT*Abody/lbody; Gyn=mu0*muY*Ayoke/lyoke;
            dd=max(abs([gshn./gsh;Gtn./Gtv;Gyn./Gyv]-1));
            gsh=gshn; Gtv=Gtn; Gyv=Gyn; nit=it;
            if dd<1e-7, break; end
        end
        R.muShoe(:,q)=muv;
    end
    X(:,q)=solveq(addiron(gsh,Gtv,Gyv),rhs); R.iters(q)=nit;
    R.Gtv=Gtv; R.Gyv=Gyv; R.gsh=gsh;
end

% ---- post-traitement ----
Xs=X(1:ns,:); Usc=Pc*Xs; Uss=Ps*Xs;
cphi=cos(nu*phis); sphi=sin(nu*phis);
Brc=AG.bru.*Usc+AG.brmq.*cphi; Brs=AG.bru.*Uss+AG.brmq.*sphi;
Btc=-AG.btu.*Uss-AG.btmq.*sphi; Bts=AG.btu.*Usc+AG.btmq.*cphi;
R.T=(L*M.rmid^2/mu0)*pi*sum(Brc.*Btc+Brs.*Bts,1);
thq=linspace(0,2*pi,3601);
Br=(Brc(:,1).'*cos(nu*thq)+Brs(:,1).'*sin(nu*thq));
Bt=(Btc(:,1).'*cos(nu*thq)+Bts(:,1).'*sin(nu*thq));
R.thq=thq; R.Br=Br; R.Bt=Bt; R.Brc=Brc; R.Brs=Brs;
R.Bmean=mean(abs(Br)); R.Bpeak=max(Br); R.Btrms=sqrt(mean(Bt.^2));
R.Bg1=hypot(Brc(M.p,1),Brs(M.p,1));
R.a8=hypot(Brc(8,1),Brs(8,1)); R.a22=hypot(Brc(22,1),Brs(22,1));
R.a23=hypot(Brc(23,1),Brs(23,1)); R.a37=hypot(Brc(37,1),Brs(37,1));
R.PhiT=R.Gtv.*(X(ns+(1:Ns),:)-X(ns+Ns+(1:Ns),:));
R.phiF=reshape(X(1:nF,1),Nf,Ns).'; R.Nf=Nf; R.X=X; R.ns=ns;
end
