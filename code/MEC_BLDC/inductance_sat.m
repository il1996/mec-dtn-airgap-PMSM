function S = inductance_sat(M, Nsurf, kfringe, ivec)
%INDUCTANCE_SAT  Inductance de LIGNE SATUREE L_ll(i) par MEC non lineaire.
%
%   POURQUOI. Au demarrage l'onduleur applique 500 V sur deux phases en
%   serie ; le courant n'est limite que par 2R et 2L. Avec l'inductance
%   NON SATUREE (L_d = 53 mH) la constante de temps vaut 5.2 ms et la crete
%   de courant plafonne a 11 A, alors que la FEA atteint 23 A en 5.1 ms.
%   L'ecart n'est pas numerique : a 23 A la FMM de bobine vaut Ntc*i = 3500
%   At, la dent est largement saturee et son inductance s'effondre. Le
%   reseau de reluctances sait le calculer -- c'est meme son point fort --
%   a condition d'etre resolu avec la vraie courbe B(H).
%
%   CE QUI EST CALCULE. Pour le motif de conduction reel (phase p en serie
%   avec phase n, courants +i et -i) :
%       lambda_ll(i) = lambda_p(i) - lambda_n(i)   [reseau non lineaire]
%       L_ll(i)      = d lambda_ll / d i           [inductance INCREMENTALE]
%   C'est exactement la grandeur qui intervient dans
%       Vdc = 2 R i + L_ll di/dt + (e_p - e_n).
%   La fuite d'encoche (air) est ajoutee analytiquement et NE sature pas.
%
%   METHODE NUMERIQUE. Les noeuds de surface (Nsurf ~ 1e3) portent
%   l'operateur DtN, dense mais INVARIANT ; seules les 2*Ns conductances de
%   dent et de culasse dependent de la saturation. On condense donc une
%   seule fois les noeuds de surface (complement de Schur K = A_ts*A_ss^-1
%   *A_st), apres quoi chaque iteration non lineaire ne coute qu'un systeme
%   30x30. Le cout total est celui d'UNE factorisation dense.

mu0=4*pi*1e-7; Rs=M.Rsi; L=M.ls; Ns=M.Ns; taus=2*pi/Ns; Ntc=M.Ntc;
PH={[1 -2 -15 3 14],[6 -7 -5 8 4],[11 -12 -10 13 9]};
BH=bh_curve();
if nargin<4||isempty(ivec), ivec=[0.001 0.25 0.5 1 1.5 2 3 4 6 8 10 13 16 20 25 30]; end

%% ---------- topologie de surface (identique a inductance_mec) -----------
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

%% ---------- blocs invariants : surface <-> dents ------------------------
%  Le bec de dent (hs0+hs1) est traite avec la permeabilite de la dent, il
%  sature donc lui aussi ; on l'ecrit g_face = mu_t * gf0.
gf0=mu0*(dth*Rs)*L/(M.hs0+M.hs1);
Gt0=mu0*M.wst1*L/M.hs2;
Gy0=mu0*M.wsy*L/(taus*(M.Rso-M.wsy/2));
Gfr=kfringe*mu0*(dth*Rs)*L/(M.ws0/2);
Ass0=-AG.Y;                                  % partie invariante
Pst=sparse(Nsurf,Ns);                        % motif surface->dent (unitaire)
for j=1:Nsurf
    if isFe(j), Pst(j,tooth(j))=Pst(j,tooth(j))+1; end             %#ok<SPRIX>
end
Pfr=sparse(Nsurf,Ns);
if Gfr>0
    for j=1:Nsurf
        if ~isFe(j)
            Pfr(j,nb1(j))=Pfr(j,nb1(j))+1; Pfr(j,nb2(j))=Pfr(j,nb2(j))+1; %#ok<SPRIX>
        end
    end
end

%% ---------- fuite d'encoche (air : NE sature pas) -----------------------
ny=4001; y=linspace(0,M.hs2,ny); b=M.ws2+(M.ws1-M.ws2)*(y/M.hs2);
lam_body=trapz(y,(y/M.hs2).^2./b);
yw=linspace(0,M.hs1,ny); bw=M.ws1+(M.ws0-M.ws1)*(yw/M.hs1);
lam_tot=lam_body+trapz(yw,1./bw)+M.hs0/M.ws0 ...
        +(5*(M.lag/M.ws0))/(5+4*(M.lag/M.ws0));
AT=zeros(Ns,3);
for k=1:3
    for c=PH{k}
        i=abs(c); s=sign(c);
        AT(mod(i-2,Ns)+1,k)=AT(mod(i-2,Ns)+1,k)+s*Ntc;
        AT(i,k)            =AT(i,k)            -s*Ntc;
    end
end
ATll=AT(:,1)-AT(:,2);                        % motif de LIGNE (a en serie -b)
Lslot_ll=mu0*L*lam_tot*sum(ATll.^2);

%% ---------- motif de FMM de branche pour la conduction a -> b -----------
Fp=zeros(Ns,1);
for c=PH{1}, Fp(abs(c))=Fp(abs(c))+sign(c)*Ntc; end
for c=PH{2}, Fp(abs(c))=Fp(abs(c))-sign(c)*Ntc; end

%% ---------- balayage en courant, resolution non lineaire ----------------
Bsat=BH.Bsat; lam_ll=zeros(size(ivec)); mu_t=zeros(Ns,numel(ivec));
Bt_=zeros(Ns,numel(ivec)); By_=zeros(Ns,numel(ivec));
mut=ones(Ns,1)*M.muI; muy=ones(Ns,1)*M.muI;   % initialisation
for q=1:numel(ivec)
    ii=ivec(q); Fs=Fp*ii;
    for it=1:60
        gface=mut*gf0;                        % une valeur par dent
        Gt=mut*Gt0;                           % dent
        Gy=zeros(Ns,1);                       % culasse : segment i -> i+1
        for i=1:Ns, Gy(i)=0.5*(muy(i)+muy(mod(i,Ns)+1))*Gy0; end
        % --- blocs ---
        gs=zeros(Nsurf,1);
        for j=1:Nsurf
            if isFe(j), gs(j)=gface(tooth(j)); elseif Gfr>0, gs(j)=2*Gfr; end
        end
        Ass=Ass0+diag(gs);
        Ast=-(Pst.*gface(tooth(:)))-Pfr*Gfr;
        tf=tooth(isFe); tf=tf(:);
        Att=diag(accumarray(tf,gface(tf),[Ns 1]))+diag(Gt);
        if Gfr>0, Att=Att+diag(full(sum(Pfr,1)).'*Gfr); end
        Aty=-diag(Gt);
        Ayy=diag(Gt);
        for i=1:Ns
            jn=mod(i,Ns)+1;
            Ayy(i,i)=Ayy(i,i)+Gy(i); Ayy(jn,jn)=Ayy(jn,jn)+Gy(i);
            Ayy(i,jn)=Ayy(i,jn)-Gy(i); Ayy(jn,i)=Ayy(jn,i)-Gy(i);
        end
        K=full(Ast.'*(Ass\full(Ast)));        % complement de Schur
        A2=[Att-K, Aty; Aty.', Ayy];
        r2=[-Gt.*Fs; Gt.*Fs];
        kp=true(2*Ns,1); kp(Ns+1)=false;      % reference : culasse 1
        U2=zeros(2*Ns,1); U2(kp)=A2(kp,kp)\r2(kp);
        Ut=U2(1:Ns); Uy=U2(Ns+1:end);
        Phi=Gt.*(Ut-Uy+Fs);                   % flux de dent
        Phy=zeros(Ns,1);
        for i=1:Ns, Phy(i)=Gy(i)*(Uy(i)-Uy(mod(i,Ns)+1)); end
        Bt=Phi/(M.wst1*L*M.Ki); By=Phy/(M.wsy*L*M.Ki);
        mun_t=BH.mur(min(abs(Bt),Bsat)); mun_y=BH.mur(min(abs(By),Bsat));
        e=max(max(abs(mun_t-mut)./mut),max(abs(mun_y-muy)./muy));
        w=0.35;                               % sous-relaxation
        mut=(1-w)*mut+w*mun_t; muy=(1-w)*muy+w*mun_y;
        if e<1e-6, break; end
    end
    lam_ll(q)=Ntc*(sum(sign(PH{1}(:)).*Phi(abs(PH{1}(:)))) ...
                  -sum(sign(PH{2}(:)).*Phi(abs(PH{2}(:)))));
    mu_t(:,q)=mut; Bt_(:,q)=Bt; By_(:,q)=By;
end

%% ---------- inductance incrementale ------------------------------------
lam_ll=lam_ll+Lslot_ll*ivec;                 % fuite d'encoche, lineaire
Linc=gradient(lam_ll,ivec);                  % d lambda / d i
S.i=ivec; S.lam=lam_ll; S.L=Linc; S.Lapp=lam_ll./ivec;
S.L0=Linc(1); S.Lslot=Lslot_ll; S.mu_t=mu_t; S.Bt=Bt_; S.By=By_;
S.Bt_max=max(abs(Bt_),[],1);
end
