function R = onload_mec(M, Nsurf, kfringe, phis, Iabc, opt)
%ONLOAD_MEC  Champ EN CHARGE : aimants + FMM d'induit + fer NON LINEAIRE.
%
%   Leve la limite dominante du modele : jusqu'ici le couple etait obtenu en
%   multipliant la back-EMF A VIDE (aimants seuls, fer lineaire) par le
%   courant. On resout ici le champ REEL en charge :
%     * source d'aimant   : injection nodale Isrc(phi) sur la surface (DtN) ;
%     * FMM des bobines   : source de branche F = +/- Ntc*i dans chaque dent
%                           (Phi = A*B((dU+F)/l)) -> REACTION D'INDUIT ;
%     * fer NON LINEAIRE  : reluctance differentielle exacte (Newton amorti),
%                           courbe B(H) de bh_curve -> SATURATION.
%   Le couple sort du tenseur de Maxwell analytique de l'operateur DtN : le
%   champ d'entrefer est entierement determine par les potentiels de surface
%   (qui portent la reaction d'induit) et par la magnetisation.
%
%   Entrees : phis (positions rotor, rad mec), Iabc (3 x Np courants [A]).
%   Sorties : R.T (couple par position), R.Tavg, R.PhiT/PhiY (flux de dent /
%             culasse), R.Us, R.Bt_max, R.iter.

if nargin<6, opt=struct(); end
tol   = getdef(opt,'tol',1e-10);
itmax = getdef(opt,'itmax',60);
mu0=4*pi*1e-7;
Rs=M.Rsi; L=M.ls; Ns=M.Ns; taus=2*pi/Ns; Ntc=M.Ntc;
PH={[1 -2 -15 3 14],[6 -7 -5 8 4],[11 -12 -10 13 9]};
BH=bh_curve();
numax=floor(Nsurf/2);

% ---- surface : fer (face de dent) / air (ouverture) -------------------
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
TB=@(i) Nsurf+i; YY=@(i) Nsurf+Ns+i; Ntot=Nsurf+2*Ns;

% ---- liste de branches : a, b, type, longueur, section, permeance ------
ba=[]; bb=[]; iron=[]; bl=[]; bA=[]; bP=[]; bcoil=[];
for j=1:Nsurf                                   % face de dent (fer)
    if isFe(j)
        ba(end+1)=j; bb(end+1)=TB(tooth(j)); iron(end+1)=1;
        bl(end+1)=M.hs0+M.hs1; bA(end+1)=dth*Rs*L; bP(end+1)=0; bcoil(end+1)=0;
    elseif kfringe>0                            % frange d'ouverture (air)
        for nb=[nb1(j) nb2(j)]
            ba(end+1)=j; bb(end+1)=TB(nb); iron(end+1)=0;
            bl(end+1)=0; bA(end+1)=0;
            bP(end+1)=kfringe*mu0*(dth*Rs)*L/(M.ws0/2); bcoil(end+1)=0;
        end
    end
end
for i=1:Ns                                      % corps de dent (fer, + FMM)
    ba(end+1)=TB(i); bb(end+1)=YY(i); iron(end+1)=1;
    bl(end+1)=M.hs2; bA(end+1)=M.wst1*L; bP(end+1)=0; bcoil(end+1)=i;
end
for i=1:Ns                                      % culasse (fer)
    ba(end+1)=YY(i); bb(end+1)=YY(mod(i,Ns)+1); iron(end+1)=1;
    bl(end+1)=taus*(M.Rso-M.wsy/2); bA(end+1)=M.wsy*L; bP(end+1)=0; bcoil(end+1)=0;
end
Nb=numel(ba); iron=logical(iron(:));
ba=ba(:); bb=bb(:); bl=bl(:); bA=bA(:); bP=bP(:); bcoil=bcoil(:);
rows=(1:Nb).';
Ainc=sparse([rows;rows],[ba;bb],[ones(Nb,1);-ones(Nb,1)],Nb,Ntot);

Yf=zeros(Ntot); Yf(1:Nsurf,1:Nsurf)=AG.Y;
ref=YY(1); keep=true(Ntot,1); keep(ref)=false; kk=find(keep);
Ared=Ainc(:,kk); Yr=Yf(kk,:); Ykk=Yf(kk,kk);

% ---- signe de bobine par dent (pour repartir la FMM) ------------------
sg=zeros(Ns,3);
for k=1:3, for c=PH{k}, sg(abs(c),k)=sign(c); end, end

Np=numel(phis);
R.T=zeros(1,Np); R.PhiT=zeros(Ns,Np); R.PhiY=zeros(Ns,Np);
R.Us=zeros(Nsurf,Np); R.iter=zeros(1,Np); R.Btmax=zeros(1,Np);
U=zeros(Ntot,1);
for q=1:Np
    F=zeros(Nb,1);
    Fi=Ntc*(sg*Iabc(:,q));                       % FMM par dent [A]
    F(bcoil>0)=Fi(bcoil(bcoil>0));
    Is=zeros(Ntot,1); Is(1:Nsurf)=AG.Isrc(phis(q));
    Isr=Is(kk);
    % ---------------- Newton amorti ----------------
    it=0;
    while it<itmax
        it=it+1;
        % Bilan nodal : flux sortant par le reseau = flux entrant de la
        % couronne (= Y*Us + Isrc), d'ou (K - Y)*U = Isrc — MEME convention
        % de SIGNE que cogging_mec (A(1:Nsurf,1:Nsurf) -= AG.Y).
        dU=Ainc*U+F;
        [Phi,Pd]=branch_flux(dU,iron,bl,bA,bP,BH);
        r=Ared.'*Phi - Yr*U - Isr;
        if max(abs(r))<max(tol, 1e-9*max(abs(Isr))), break; end
        J=Ared.'*spdiags(Pd,0,Nb,Nb)*Ared - Ykk;
        dx=-(J\r); lam=1; rn0=norm(r);
        for ls=1:30
            Ut=U; Ut(kk)=U(kk)+lam*dx;
            [Pt,~]=branch_flux(Ainc*Ut+F,iron,bl,bA,bP,BH);
            rt=Ared.'*Pt - Yr*Ut - Isr;
            if norm(rt)<rn0 || lam<1e-6, U=Ut; break; end
            lam=lam/2;
        end
    end
    R.iter(q)=it;
    dU=Ainc*U+F; [Phi,~]=branch_flux(dU,iron,bl,bA,bP,BH);
    it0=find(bcoil>0); R.PhiT(:,q)=Phi(it0);
    iy=Nb-Ns+1:Nb;     R.PhiY(:,q)=Phi(iy);
    R.Us(:,q)=U(1:Nsurf);
    R.Btmax(q)=max(abs(Phi(it0)))/(M.wst1*L);
    % ---- couple : tenseur de Maxwell analytique (DtN) ----
    nu=AG.nu; Usc=AG.Wc*U(1:Nsurf); Uss=AG.Ws*U(1:Nsurf);
    cp=cos(nu*phis(q)); sp=sin(nu*phis(q));
    Brc=AG.bru.*Usc+AG.brmq.*cp;  Brs=AG.bru.*Uss+AG.brmq.*sp;
    Btc=-AG.btu.*Uss-AG.btmq.*sp; Bts=AG.btu.*Usc+AG.btmq.*cp;
    R.T(q)=(L*M.rmid^2/mu0)*pi*sum(Brc.*Btc+Brs.*Bts);
end
R.Tavg=mean(R.T); R.phis=phis; R.AG=AG; R.Nsurf=Nsurf;
end

% ======================================================================
function [Phi,Pd]=branch_flux(dU,iron,bl,bA,bP,BH)
Phi=zeros(size(dU)); Pd=Phi;
H=dU(iron)./bl(iron); B=BH.Bof(H);
Phi(iron)=B.*bA(iron);
Pd(iron)=bA(iron)./(bl(iron).*BH.dHdB(B));
Phi(~iron)=bP(~iron).*dU(~iron);
Pd(~iron)=bP(~iron);
end

function v=getdef(s,f,d)
if isfield(s,f)&&~isempty(s.(f)), v=s.(f); else, v=d; end
end
