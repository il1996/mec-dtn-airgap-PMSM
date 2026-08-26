function R = inductance_mec2(M, Nsurf, nm, muI)
%INDUCTANCE_MEC2  Inductances avec le SOUS-MODELE 2D de bouche d'encoche.
%
%   TEST DECISIF du sous-modele : avec la bouche MAILLEE, le pont
%   TANGENTIEL dent-a-dent (= la fuite de bec) EMERGE de la geometrie. On ne
%   doit donc PLUS ajouter lambda_tip analytiquement - contrairement au
%   modele a scalaire kfringe, ou il fallait choisir l'un OU l'autre sous
%   peine de double comptage (La : -1.7 % -> +29 %).
%   Si La et M restent justes SANS lambda_tip, le sous-modele est valide :
%   une seule geometrie sert les deux problemes (a vide et en induit).

mu0=4*pi*1e-7;
Rs=M.Rsi; L=M.ls; Ns=M.Ns; taus=2*pi/Ns; Ntc=M.Ntc;
PH={[1 -2 -15 3 14],[6 -7 -5 8 4],[11 -12 -10 13 9]};
numax=floor(Nsurf/2);

dth=2*pi/Nsurf; ths=(0:Nsurf-1)*dth; dths=dth*ones(1,Nsurf);
wtf=taus-M.ws0/Rs;
tooth=zeros(1,Nsurf); isFe=false(1,Nsurf); slot=zeros(1,Nsurf);
for j=1:Nsurf
    i0=mod(round(ths(j)/taus),Ns); dd=angle(exp(1i*(ths(j)-i0*taus)));
    tooth(j)=i0+1;
    if abs(dd)<wtf/2, isFe(j)=true;
    elseif dd>0,      slot(j)=i0+1;
    else,             slot(j)=mod(i0-1+Ns,Ns)+1; end
end
AG=airgap_magnet(M,ths,dths,numax);

D=M.hs0+M.hs1; dd_=D/nm; dc=((1:nm)-0.5)*dd_;
wof=@(d) M.ws0 + (d>M.hs0).*(M.ws1-M.ws0).*max(d-M.hs0,0)/M.hs1;
wc=wof(dc); wi=wof((1:nm-1)*dd_);

TT=@(i) Nsurf+i; TB=@(i) Nsurf+Ns+i; YY=@(i) Nsurf+2*Ns+i;
MM=@(k,j) Nsurf+3*Ns+(k-1)*nm+j;
Ntot=Nsurf+3*Ns+Ns*nm; A=zeros(Ntot);
    function add(a,b,g), A(a,a)=A(a,a)+g; A(b,b)=A(b,b)+g; A(a,b)=A(a,b)-g; A(b,a)=A(b,a)-g; end

g_face=mu0*muI*(dth*Rs)*L/(D/2);
w_shoe=0.5*((taus*Rs-M.ws0)+(2*pi*(Rs+D)/Ns-M.ws1));
G_shoe=mu0*muI*w_shoe*L/D;
Gt=mu0*muI*M.wst1*L/M.hs2;
Gy=mu0*muI*M.wsy*L/(taus*(M.Rso-M.wsy/2));
for j=1:Nsurf, if isFe(j), add(j,TT(tooth(j)),g_face); end, end
for i=1:Ns
    add(TT(i),TB(i),G_shoe); add(TB(i),YY(i),Gt); add(YY(i),YY(mod(i,Ns)+1),Gy);
end
for j=1:Nsurf
    if ~isFe(j), add(j,MM(slot(j),1), mu0*(dth*Rs)*L/(dd_/2)); end
end
for k=1:Ns
    for j=1:nm-1, add(MM(k,j),MM(k,j+1), mu0*wi(j)*L/dd_); end
    for j=1:nm
        Gl=2*mu0*dd_*L/wc(j);
        add(MM(k,j),TT(k),Gl); add(MM(k,j),TT(mod(k,Ns)+1),Gl);
    end
end
A(1:Nsurf,1:Nsurf)=A(1:Nsurf,1:Nsurf)-AG.Y;
ref=YY(1); keep=true(Ntot,1); keep(ref)=false;

% ---- excitation : FMM de bobine en source de branche dans chaque dent ----
I0=1; Fs=zeros(Ns,1);
for c=PH{1}, Fs(abs(c))=sign(c)*Ntc*I0; end
rhs=zeros(Ntot,1);
for i=1:Ns
    rhs(TB(i))=rhs(TB(i))-Gt*Fs(i);
    rhs(YY(i))=rhs(YY(i))+Gt*Fs(i);
end
U=zeros(Ntot,1); U(keep)=A(keep,keep)\rhs(keep);
Phi=Gt*(U(TB(1:Ns))-U(YY(1:Ns))+Fs);
lamk=zeros(1,3);
for k=1:3, P=PH{k}; lamk(k)=Ntc*sum(sign(P(:)).*Phi(abs(P(:)))); end
R.L_gap=lamk(1)/I0; R.M_gap=0.5*(lamk(2)+lamk(3))/I0;

% ---- fuite d'encoche SANS lambda_tip (la bouche le fournit) ------------
ny=4001; y=linspace(0,M.hs2,ny);
b=M.ws2+(M.ws1-M.ws2)*(y/M.hs2);
lam_body=trapz(y,(y/M.hs2).^2./b);
yw=linspace(0,M.hs1,ny); bw=M.ws1+(M.ws0-M.ws1)*(yw/M.hs1);
lam_wedge=trapz(yw,1./bw);
% DOUBLE COMPTAGE A EVITER : la bouche maillee fournit DEJA, pour le trajet
% dent-a-dent, exactement la permeance d'ouverture hs0/ws0 = 0.500 et celle
% du biseau = 0.129 (somme des branches laterales : sum(mu0*dz*L/w(d)) =
% mu0*L*INT dd/w = mu0*L*(lambda_ouv + lambda_biseau), verifie). Il ne reste
% donc a ajouter analytiquement que le CORPS d'encoche.
lam_slot=lam_body;                            % ni ouverture, ni biseau, ni bec
R.lam_slot=lam_slot; R.lam_mouth_net=M.hs0/M.ws0+lam_wedge;
AT=zeros(Ns,3);
for k=1:3
    for c=PH{k}
        i=abs(c); s=sign(c);
        AT(mod(i-2,Ns)+1,k)=AT(mod(i-2,Ns)+1,k)+s*Ntc;
        AT(i,k)=AT(i,k)-s*Ntc;
    end
end
R.L_slot=mu0*L*lam_slot*sum(AT(:,1).^2);
R.M_slot=mu0*L*lam_slot*0.5*(sum(AT(:,1).*AT(:,2))+sum(AT(:,1).*AT(:,3)));
R.La=R.L_gap+R.L_slot; R.M=R.M_gap+R.M_slot; R.Ld=R.La-R.M;
R.nm=nm; R.Nsurf=Nsurf;
end
