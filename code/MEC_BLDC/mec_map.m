function S = mec_map(M, Nsurf, kfringe, Nth, iax, Np_pm)
%MEC_MAP  Cartographie NON LINEAIRE flux/couple du reseau MEC : lambda(th,i), T(th,i).
%
%   POURQUOI UNE CARTOGRAPHIE. Au demarrage l'onduleur impose 500 V sur deux
%   phases en serie et le courant atteint 23 A : la FMM de bobine vaut alors
%   Ntc*i = 3500 At et le flux d'induit represente 4.7 fois le flux des
%   aimants. Dans ce regime la denture est saturee par l'induit et le flux
%   des aimants vu par la phase s'effondre : ni lambda_pm ni L ne sont plus
%   des constantes. Un modele lineaire predit 11 A la ou la FEA en donne 23.
%   On calcule donc, avec le MEME reseau (DtN etendu + denture + culasse)
%   mais resolu avec la vraie courbe B(H) M350-50A, la surface
%       lambda_k(theta, i_alpha, i_beta)   et   T(theta, i_alpha, i_beta)
%   sur laquelle le modele d'entrainement viendra interpoler.
%
%   MISE EN OEUVRE. Les Nsurf noeuds de surface portent l'operateur DtN,
%   dense mais INVARIANT (la geometrie de la couronne ne bouge pas, seule la
%   source de magnetisation tourne) ; la saturation n'affecte que les 2*Ns
%   conductances de dent et de culasse. On condense donc les noeuds de
%   surface UNE SEULE FOIS (complement de Schur K = A_ts*A_ss^-1*A_st) et
%   chaque iteration non lineaire ne coute plus qu'un systeme 30x30. Le bec
%   de dent (hs0+hs1 = 1.5 mm contre hs2 = 19.6 mm) est laisse lineaire pour
%   que A_ss et K restent invariants : approximation assumee et documentee.
%
%   Le couple est obtenu par le TENSEUR DE MAXWELL sur les harmoniques de
%   surface (forme de Parseval, identique a cogging_mec), donc valable en
%   regime sature ou l'expression sum(i*dlambda/dtheta) ne l'est plus.

mu0=4*pi*1e-7; Rs=M.Rsi; L=M.ls; Ns=M.Ns; taus=2*pi/Ns; Ntc=M.Ntc; p=M.p;
PH={[1 -2 -15 3 14],[6 -7 -5 8 4],[11 -12 -10 13 9]};
BH=bh_curve();
if nargin<2||isempty(Nsurf),  Nsurf=720;  end
if nargin<3||isempty(kfringe), kfringe=0.75; end
if nargin<4||isempty(Nth),    Nth=37;     end
if nargin<5||isempty(iax),    iax=linspace(-28,28,15); end
if nargin<6||isempty(Np_pm),  Np_pm=Nth;  end   %#ok<NASGU>

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

%% ---------- blocs invariants (bec de dent + frange : lineaires) ---------
gf =mu0*M.muI*(dth*Rs)*L/(M.hs0+M.hs1);
Gt0=mu0*M.wst1*L/M.hs2;                       % a multiplier par mu_dent
Gy0=mu0*M.wsy*L/(taus*(M.Rso-M.wsy/2));       % a multiplier par mu_culasse
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
K=Ast.'*(dAss\Ast);                           % complement de Schur [Ns x Ns]

%% ---------- source des aimants, condensee ------------------------------
the=linspace(0,2*pi,Nth+1); the(end)=[];      % angle ELECTRIQUE
phi=the/p;                                    % position mecanique
brm=AG.brm;
Isr=L*Rs*pi*( AG.Wc.'*(brm.*cos(nu*phi)) + AG.Ws.'*(brm.*sin(nu*phi)) );
Cpm=Ast.'*(dAss\Isr);                         % [Ns x Nth]

%% ---------- motifs de FMM de branche par phase --------------------------
Fp=zeros(Ns,3);
for k=1:3, for c=PH{k}, Fp(abs(c),k)=Fp(abs(c),k)+sign(c)*Ntc; end, end

%% ---------- fuite d'encoche (air, lineaire, hors reseau) ----------------
ny=4001; y=linspace(0,M.hs2,ny); b=M.ws2+(M.ws1-M.ws2)*(y/M.hs2);
yw=linspace(0,M.hs1,ny); bw=M.ws1+(M.ws0-M.ws1)*(yw/M.hs1);
lam_tot=trapz(y,(y/M.hs2).^2./b)+trapz(yw,1./bw)+M.hs0/M.ws0 ...
        +(5*(M.lag/M.ws0))/(5+4*(M.lag/M.ws0));
AT=zeros(Ns,3);
for k=1:3
    for c=PH{k}
        i=abs(c); s=sign(c);
        AT(mod(i-2,Ns)+1,k)=AT(mod(i-2,Ns)+1,k)+s*Ntc;
        AT(i,k)            =AT(i,k)            -s*Ntc;
    end
end
Lsl=mu0*L*lam_tot*(AT.'*AT);                  % matrice de fuite [3x3]

%% ---------- balayage (theta, i_alpha, i_beta) ---------------------------
na=numel(iax); Bsat=BH.Bsat; muI=M.muI;
LAM=zeros(3,Nth,na,na); TQ=zeros(Nth,na,na);
PHT=zeros(Ns,Nth,na,na); PHY=zeros(Ns,Nth,na,na); LSLQ=zeros(Nth,na,na);
kp=true(2*Ns,1); kp(Ns+1)=false;              % reference : culasse 1
iY=Ns+(1:Ns); nx=[2:Ns 1]; pv=[Ns 1:Ns-1];
idg=sub2ind([Ns Ns],(1:Ns).',(1:Ns).');       % diagonale du bloc culasse
iup=sub2ind([Ns Ns],(1:Ns).',nx.');
idn=sub2ind([Ns Ns],nx.',(1:Ns).');
Kdg=-K+diag(dgt);                             % bloc dent invariant
nit=0;
for qa=1:na
  for qb=1:na
    ia=iax(qa); ib=iax(qb);
    iabc=[ia; -ia/2+sqrt(3)/2*ib; -ia/2-sqrt(3)/2*ib];   % Clarke inverse
    Fs0=Fp*iabc;
    mut=ones(Ns,1)*muI; muy=mut;
    for q=1:Nth
        for it=1:80
            Gt=mut*Gt0;
            Gy=0.5*(muy+muy(nx))*Gy0;         % segment i : culasse i -> i+1
            Ayy=zeros(Ns); Ayy(idg)=Gt+Gy+Gy(pv);
            Ayy(iup)=Ayy(iup)-Gy; Ayy(idn)=Ayy(idn)-Gy;
            DG=diag(Gt);
            A2=[Kdg+DG, -DG; -DG, Ayy];       % Kdg = -K + diag(dgt), invariant
            r2=[-Gt.*Fs0-Cpm(:,q); Gt.*Fs0];
            U2=zeros(2*Ns,1); U2(kp)=A2(kp,kp)\r2(kp);
            Ut=U2(1:Ns); Uy=U2(iY);
            Phi=Gt.*(Ut-Uy+Fs0);
            Phy=Gy.*(Uy-Uy([2:Ns 1]));
            Bt=Phi/(M.wst1*L*M.Ki); By=Phy/(M.wsy*L*M.Ki);
            %  PAS d'ecretage a Bsat : au-dela du coude la tole suit
            %  B = mu0*H + Js, et bh_curve prolonge deja H(B) a la pente
            %  1/mu0. Ecreter reviendrait a bloquer mu_r a 45 alors qu'a
            %  2.25 T il vaut 4.9 -- c'est precisement dans ce regime que
            %  le flux des aimants est chasse de la dent.
            nt=BH.mur(abs(Bt)); nyk=BH.mur(abs(By));
            er=max(max(abs(nt-mut)./mut),max(abs(nyk-muy)./muy));
            %  sous-relaxation en log : mu varie sur 3 decades (3000 -> 5)
            w=0.35; mut=exp((1-w)*log(mut)+w*log(nt));
            muy=exp((1-w)*log(muy)+w*log(nyk));
            nit=nit+1;
            if er<1e-6, break; end
        end
        % --- SATURATION DE LA FUITE D'ENCOCHE -------------------------
        %  La permeance d'encoche classique lambda = h/(3b) suppose un fer
        %  de permeabilite INFINIE de part et d'autre de l'encoche. Elle
        %  pese ici 63.8 mH sur les 106 mH d'inductance de ligne, soit
        %  60 % : la laisser en air pur revient a interdire a l'inductance
        %  de descendre sous 64 mH, quel que soit le courant. Or le flux de
        %  fuite lui-meme traverse les dents : a 10 A il y ajoute ~0.9 T,
        %  la dent depasse 2.3 T et sa permeabilite tombe a ~5. On met donc
        %  la permeance d'air EN SERIE avec la reluctance du fer qui ferme
        %  son trajet, evaluee avec la permeabilite deja calculee par le
        %  reseau non lineaire :
        %      lambda_eff = lambda_air / (1 + lambda_air*mu0*L*R_fer)
        %      R_fer      = (2/3)*hs2/(mu0*mu_dent*wst1*L*Ki)
        %  le facteur 2/3 traduit la meme moyenne de FMM distribuee que le
        %  h/(3b) de la permeance d'encoche (montee dans une dent, descente
        %  dans la voisine).
        mub=0.5*(mut+mut(nx));                   % dents bordant l'encoche
        Rfe=(2/3)*M.hs2./(mu0*mub*M.wst1*L*M.Ki);
        lam_ef=lam_tot./(1+lam_tot*mu0*L*Rfe);   % Ns x 1
        Lslq=mu0*L*(AT.'*(lam_ef.*AT));          % matrice de fuite 3x3
        % --- flux totalise (reseau + fuite d'encoche saturable) --------
        for k=1:3
            LAM(k,q,qa,qb)=Ntc*sum(sign(PH{k}(:)).*Phi(abs(PH{k}(:))));
        end
        LAM(:,q,qa,qb)=LAM(:,q,qa,qb)+Lslq*iabc;
        LSLQ(q,qa,qb)=lam_ef(1)/lam_tot;         % taux de saturation de fuite
        PHT(:,q,qa,qb)=Phi; PHY(:,q,qa,qb)=Phy;
        % --- couple : tenseur de Maxwell sur les harmoniques de surface ---
        Us=dAss\(Isr(:,q)-Ast*Ut);
        Usc=AG.Wc*Us; Uss=AG.Ws*Us;
        cq=cos(nu*phi(q)); sq=sin(nu*phi(q));
        Brc=AG.bru.*Usc+AG.brmq.*cq;  Brs=AG.bru.*Uss+AG.brmq.*sq;
        Btc=-AG.btu.*Uss-AG.btmq.*sq; Bts=AG.btu.*Usc+AG.btmq.*cq;
        TQ(q,qa,qb)=(L*M.rmid^2/mu0)*pi*sum(Brc.*Btc+Brs.*Bts);
    end
  end
end

%% ---------- derivees ----------------------------------------------------
%  d lambda / d theta_m : derivee SPECTRALE (periodicite exacte en theta_e)
%  En regime sature lambda(theta) presente des coudes ; une derivee
%  spectrale brute sonne (phenomene de Gibbs). On tronque donc le spectre
%  aux harmoniques REELLEMENT portes par le bobinage (au-dela de l'ordre 25
%  il ne reste que le bruit de la grille en theta).
PSI=zeros(size(LAM));
if mod(Nth,2)==1, kk=[0:(Nth-1)/2, -(Nth-1)/2:-1];
else,             kk=[0:Nth/2-1, 0, -Nth/2+1:-1]; end
%  Le bobinage ne peut porter que les rangs impairs jusqu'a ~13 (au-dela
%  le facteur de bobinage est negligeable) : tronquer la, c'est retenir
%  toute la physique et rejeter tout le bruit de grille.
kmax=min(13,floor(Nth/2)-1); msk=abs(kk)<=kmax;
for k=1:3, for qa=1:na, for qb=1:na
    F=fft(LAM(k,:,qa,qb));
    PSI(k,:,qa,qb)=p*real(ifft(1i*kk.*F.*msk)); % d/d theta_m
end, end, end
%  d lambda / d i_alpha , d i_beta
LA=zeros(size(LAM)); LB=zeros(size(LAM));
for k=1:3, for q=1:Nth
    Z=squeeze(LAM(k,q,:,:));
    [gb,ga]=gradient(Z,iax,iax);              % dim1 = i_alpha, dim2 = i_beta
    LA(k,q,:,:)=ga; LB(k,q,:,:)=gb;
end, end

%% ---------- correction du DOUBLE COMPTAGE de la fuite de bec -----------
%  La conductance de frange kfringe joue deux roles distincts (voir l'entete
%  d'inductance_mec) : a vide elle porte l'entree RADIALE du flux par la
%  bouche d'encoche -- c'est pour cela qu'elle est identifiee a 0.75 et
%  qu'il faut la garder ici, la cartographie etant d'abord celle du champ
%  des aimants ; en induit elle cree en plus un pont TANGENTIEL dent a dent
%  qui fait DOUBLE EMPLOI avec la fuite de bec lambda_tip deja comptee
%  analytiquement. Sur l'inductance de ligne a courant nul cela represente
%  +26 % (133.7 mH au lieu des 106.1 mH valides a +1.3 % sur la FEA).
%  L'exces est un chemin dans l'AIR : il ne sature pas, la correction est
%  donc SOUSTRACTIVE et de meme forme que la matrice de fuite d'encoche.
ATll=AT(:,1)-AT(:,2);
RIref=inductance_mec(M,Nsurf,M.muI,0);        % reference : kfringe = 0
i0=find(abs(iax)<1e-12,1); if isempty(i0), [~,i0]=min(abs(iax)); end
LAA=squeeze(LA(1,:,i0,i0)); LBA=squeeze(LA(2,:,i0,i0));
LAB=squeeze(LB(1,:,i0,i0)); LBB=squeeze(LB(2,:,i0,i0));
Leq0=mean((LAA-LBA)-(LAB-LBB)/sqrt(3));
alpha=(Leq0-2*RIref.Ld)/(mu0*L*sum(ATll.^2));
dLsl=alpha*mu0*L*(AT.'*AT);
for qa=1:na, for qb=1:na
    ia=iax(qa); ib=iax(qb);
    iabc=[ia; -ia/2+sqrt(3)/2*ib; -ia/2-sqrt(3)/2*ib];
    LAM(:,:,qa,qb)=LAM(:,:,qa,qb)-dLsl*iabc;
end, end
for k=1:3, for q=1:Nth
    Z=squeeze(LAM(k,q,:,:)); [gb,ga]=gradient(Z,iax,iax);
    LA(k,q,:,:)=ga; LB(k,q,:,:)=gb;
end, end
S.alpha=alpha; S.Leq0=Leq0; S.Ldref=RIref.Ld; S.dLsl=dLsl; S.fslot=LSLQ;

S.the=the; S.iax=iax; S.lam=LAM; S.psi=PSI; S.La=LA; S.Lb=LB; S.T=TQ;
S.PhiT=PHT; S.PhiY=PHY; S.Nsurf=Nsurf; S.kfringe=kfringe; S.p=p;
S.Lsl=Lsl; S.nit=nit;
end
