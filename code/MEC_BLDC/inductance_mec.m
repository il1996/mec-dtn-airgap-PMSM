function R = inductance_mec(M, Nsurf, muI, kfringe)
%INDUCTANCE_MEC  Inductances propre et mutuelle par MEC couple (amelioration 3).
%
%   ATTENTION - ROLE DE kfringe POUR LE PROBLEME D'INDUIT (a lire avant usage)
%   La conductance de frange de l'ouverture d'encoche joue DEUX roles
%   physiquement DISTINCTS, qu'un seul scalaire ne peut pas porter :
%     (a) A VIDE : entree RADIALE du flux d'entrefer par la bouche d'encoche
%         -> c'est ce que RUN_SLOT_IDENT identifie (kfringe = 0.75) sur les
%            bandes de denture du champ FEA ;
%     (b) EN INDUIT : pont TANGENTIEL dent-a-dent a travers la bouche
%         d'encoche -> c'est la fuite de BEC DE DENT, deja comptee
%            analytiquement par lambda_tip dans la fuite d'encoche.
%   Le pont total du reseau vaut kfringe*mu0*L par encoche (independant de
%   Nsurf, verifie) contre mu0*L*lambda_tip = 0.357*mu0*L pour la formule
%   classique : reutiliser kfringe = 0.75 ici DOUBLE-COMPTE la fuite de bec
%   (La passe de 49.4 a 64.8 mH, soit -1.7 % -> +29 % vs FEA).
%   => Pour l'induit on prend kfringe = 0 (pas de pont dans le reseau) et la
%      fuite de bec reste analytique. Les deux geometries (radiale peu
%      profonde / tangentielle sur ws0) sont differentes : c'est une limite
%      de modele assumee et documentee, pas un recalage.
%
%   Deux contributions, calculees SEPAREMENT et rigoureusement :
%
%   (1) PART D'ENTREFER — le champ d'induit est resolu par le RESEAU couple
%       a l'operateur DtN (aimants eteints) : la FMM de chaque bobine
%       dentaire est injectee comme SOURCE DE BRANCHE dans sa dent
%       (Phi = G*(Ua-Ub+F)), la denture et la frange d'encoche sont donc
%       prises en compte, au lieu d'un entrefer uniforme + facteur de Carter.
%
%   (2) FUITE D'ENCOCHE — FMM DISTRIBUEE, integree sur la VRAIE geometrie
%       trapezoidale (largeur variable ws2 au fond -> ws1 au sommet) :
%           lambda_corps = INT_0^hs2 (y/hs2)^2 / b(y) dy      [y depuis la culasse]
%           lambda_biseau= INT dy/b(y) sur hs1     (FMM totale enlacee)
%           lambda_ouvert= hs0/ws0
%       La formule usuelle a LARGEUR MOYENNE hs2/(3*b_moy) sous-estime le
%       trapeze (le poids (y/h)^2 est maximal la ou la largeur est MINIMALE).
%       La sommation se fait sur la FMM REELLE PAR ENCOCHE : pour un bobinage
%       dentaire, deux cotes de bobine voisins peuvent appartenir a la MEME
%       phase et leurs ampere-tours S'AJOUTENT (encoches 15,1,2,14 ici).
%
%   La mutuelle est obtenue par le meme calcul, ce qui permet de confronter
%   La ET M (deux cibles FEA independantes) et donc L_d = L_a - M.

mu0=4*pi*1e-7;
Rs=M.Rsi; L=M.ls; Ns=M.Ns; taus=2*pi/Ns; Ntc=M.Ntc;
PH={[1 -2 -15 3 14],[6 -7 -5 8 4],[11 -12 -10 13 9]};

%% ---------- (1) part d'entrefer par le reseau couple -------------------
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
AG=airgap_magnet(M,ths,dths,numax);

g_face=mu0*muI*(dth*Rs)*L/(M.hs0+M.hs1);
Gt=mu0*muI*M.wst1*L/M.hs2;
Gy=mu0*muI*M.wsy*L/(taus*(M.Rso-M.wsy/2));
Gfr=kfringe*mu0*(dth*Rs)*L/(M.ws0/2);
Ntot=Nsurf+2*Ns; TB=@(i)Nsurf+i; YY=@(i)Nsurf+Ns+i;
A=zeros(Ntot);
    function add(a,b,g), A(a,a)=A(a,a)+g; A(b,b)=A(b,b)+g; A(a,b)=A(a,b)-g; A(b,a)=A(b,a)-g; end
for j=1:Nsurf
    if isFe(j), add(j,TB(tooth(j)),g_face);
    elseif Gfr>0, add(j,TB(nb1(j)),Gfr); add(j,TB(nb2(j)),Gfr); end
end
for i=1:Ns, add(TB(i),YY(i),Gt); add(YY(i),YY(mod(i,Ns)+1),Gy); end
A(1:Nsurf,1:Nsurf)=A(1:Nsurf,1:Nsurf)-AG.Y;
ref=YY(1); keep=true(Ntot,1); keep(ref)=false;

I0=1;                                        % courant d'essai [A]
%  Les TROIS phases sont excitees separement : on obtient la reponse de
%  denture a 1 A par phase, necessaire (a) aux inductances, (b) au champ
%  d'INDUIT en charge (superposition lineaire dans les pertes fer).
PhiU=zeros(Ns,3);                            % flux de dent par ampere, phase k
Usf=zeros(Nsurf,3);                          % potentiel de SURFACE par ampere
for k=1:3
    Fs=zeros(Ns,1);                          % FMM de branche par dent
    for c=PH{k}, Fs(abs(c))=sign(c)*Ntc*I0; end
    rhs=zeros(Ntot,1);
    for i=1:Ns                               % branche TB(i)->YY(i), source Fs(i)
        rhs(TB(i))=rhs(TB(i))-Gt*Fs(i);
        rhs(YY(i))=rhs(YY(i))+Gt*Fs(i);
    end
    U=zeros(Ntot,1); U(keep)=A(keep,keep)\rhs(keep);
    PhiU(:,k)=Gt*(U(TB(1:Ns))-U(YY(1:Ns))+Fs);
    Usf(:,k)=U(1:Nsurf);
end
Phi=PhiU(:,1);                               % phase A : reference historique
R.PhiU=PhiU/I0;                              % Wb par ampere de phase [Ns x 3]
%  Potentiels de surface d'INDUIT par ampere : source du champ d'induit dans
%  la couronne d'aimant (pertes aimant en charge, via le DtN etendu).
%  Ils ne dependent PAS de la position rotor : la geometrie de l'anneau est
%  invariante, seule la source de magnetisation tourne.
R.Usurf=Usf/I0; R.ths=ths; R.AG=AG;

lamk=zeros(1,3);
for k=1:3, P=PH{k}; lamk(k)=Ntc*sum(sign(P(:)).*Phi(abs(P(:)))); end
R.L_gap=lamk(1)/I0;  R.M_gap=0.5*(lamk(2)+lamk(3))/I0;

%% ---------- (2) fuite d'encoche : FMM distribuee, trapeze exact --------
ny=4001;
y=linspace(0,M.hs2,ny);                      % y = 0 a la culasse, hs2 au sommet
b=M.ws2+(M.ws1-M.ws2)*(y/M.hs2);             % largeur : ws2 (fond) -> ws1 (sommet)
lam_body=trapz(y,(y/M.hs2).^2./b);
yw=linspace(0,M.hs1,ny);
bw=M.ws1+(M.ws0-M.ws1)*(yw/M.hs1);           % biseau : ws1 -> ws0
lam_wedge=trapz(yw,1./bw);
lam_open=M.hs0/M.ws0;
lam_slot=lam_body+lam_wedge+lam_open;
lam_tip=(5*(M.lag/M.ws0))/(5+4*(M.lag/M.ws0));
lam_tot=lam_slot+lam_tip;
R.lam_body=lam_body; R.lam_wedge=lam_wedge; R.lam_slot=lam_slot; R.lam_tot=lam_tot;
R.lam_body_approx=M.hs2/(3*(M.ws1+M.ws2)/2);

% ampere-tours par ENCOCHE : bobine de la dent i -> +s dans l'encoche i-1,
% -s dans l'encoche i (encoche k = entre dent k et dent k+1)
AT=zeros(Ns,3);
for k=1:3
    for c=PH{k}
        i=abs(c); s=sign(c);
        AT(mod(i-2,Ns)+1,k)=AT(mod(i-2,Ns)+1,k)+s*Ntc;
        AT(i,k)            =AT(i,k)            -s*Ntc;
    end
end
R.AT=AT;
R.L_slot=mu0*L*lam_tot*sum(AT(:,1).^2);
R.M_slot=mu0*L*lam_tot*0.5*(sum(AT(:,1).*AT(:,2))+sum(AT(:,1).*AT(:,3)));

%% ---------- totaux ----------------------------------------------------
R.La=R.L_gap+R.L_slot;  R.M=R.M_gap+R.M_slot;  R.Ld=R.La-R.M;
R.Nsurf=Nsurf; R.muI=muI; R.kfringe=kfringe;
end
