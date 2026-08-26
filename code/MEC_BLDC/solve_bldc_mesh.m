function S = solve_bldc_mesh(ME, tol, itmax, U0)
%SOLVE_BLDC_MESH  Resolution du reseau 2D raffine + DtN etendu.
%
%   Systeme en POTENTIEL SCALAIRE aux noeuds : pour chaque branche
%   Phi = G*(Ua - Ub + F), avec G = mu0*mu_r*A/l pour le fer et G = P pour
%   l'air. La loi des noeuds donne A*U = rhs. Les Ms noeuds de SURFACE
%   portent en plus l'operateur DtN de la couronne entrefer + aimant, dense
%   mais exact, et la source de magnetisation.
%
%   Non lineaire : la permeabilite de chaque branche de FER est reprise sur
%   la courbe B(H) reelle a partir de son induction, avec sous-relaxation
%   logarithmique (mu varie sur trois decades entre 0 et 2.3 T).

%   METHODE. Le point fixe sur mu_r stagne : dans les cellules de bec la
%   permeabilite saute de 6700 a 65 d'une iteration a l'autre et l'iteration
%   oscille. On resout donc par NEWTON EXACT sur la reluctance
%   DIFFERENTIELLE. La branche de fer est une resistance non lineaire
%       DeltaU = H(B) * l ,   B = Phi/A
%   qui s'INVERSE directement -- B = Bof(DeltaU/l) -- puisque bh_curve
%   fournit B(H). Le jacobien utilise la permeance incrementale
%       G_diff = A / (l * dH/dB)
%   et non A*mu_r/l. La convergence devient quadratique.

if nargin<2||isempty(tol),   tol=1e-6; end
if nargin<3||isempty(itmax), itmax=60; end
mu0=4*pi*1e-7; BH=bh_curve();
nb=numel(ME.a); N=ME.N;
nlin=ischar(ME.muI)||isstring(ME.muI);
mur=ones(nb,1); if ~nlin, mur(:)=ME.muI; else, mur(:)=M_default(); end
mur(~ME.iron)=1;
Ysurf=ME.AG.Y; ids=ME.Ng0+(1:ME.Ms);
ref=1;                                          % noeud de reference

kp=true(N,1); kp(ref)=false;
Fe=ME.iron; Ai=ME.A; li=ME.l; Pk=ME.P; Ek=ME.E;
asm=@(g)assemble(g,ME.a,ME.b,N,Ysurf,ids);

%% ---------- solution LINEAIRE (aussi le point de depart de Newton) ------
G=zeros(nb,1);
G(Fe)=mu0*mur(Fe).*Ai(Fe)./li(Fe); G(~Fe)=Pk(~Fe);
A=asm(G);
rhs=ME.Isrc+accumarray(ME.a,-G.*Ek,[N 1])+accumarray(ME.b,G.*Ek,[N 1]);
U=zeros(N,1); U(kp)=A(kp,kp)\rhs(kp);
%  DEMARRAGE A CHAUD. Sur un balayage rotor, la solution de la position
%  precedente est un bien meilleur point de depart que la solution lineaire :
%  Newton passe de 9 iterations a 2 ou 3.
if nargin>=4 && ~isempty(U0) && numel(U0)==N, U=U0; end

if nlin
    %% ------------------ NEWTON EXACT ---------------------------------
    nrm0=[]; er=inf;
    for it=1:itmax
        [f,Phi,B,Gd]=resid(U,ME,BH,mu0,Ysurf,ids,N);
        nf=norm(f(kp));
        if isempty(nrm0), nrm0=max(nf,eps); end
        er=nf/nrm0;
        if er<tol, break; end
        J=asm(Gd);
        d=zeros(N,1); d(kp)=-J(kp,kp)\f(kp);
        %  recherche lineaire : le pas complet peut sortir du domaine de
        %  la courbe B(H) au premier abord
        al=1; nf2=inf;
        for ls=1:12
            f2=resid(U+al*d,ME,BH,mu0,Ysurf,ids,N);
            nf2=norm(f2(kp));
            if nf2<nf, break; end
            al=al/2;
        end
        U=U+al*d;
    end
    [~,Phi,B]=resid(U,ME,BH,mu0,Ysurf,ids,N);
    mur=ones(nb,1); m=Fe & abs(B)>1e-9;
    mur(m)=abs(B(m))./(mu0*max(BH.Hof(abs(B(m))),1e-9));
else
    Phi=G.*(U(ME.a)-U(ME.b)+Ek);
    B=zeros(nb,1); m=Ai>0; B(m)=Phi(m)./Ai(m); er=0; it=1;
end
S.U=U; S.Phi=Phi; S.B=B; S.mur=mur; S.iter=it; S.err=er; S.Usurf=U(ids);
end

%% ------------------------------------------------------------------------
function A=assemble(G,a,b,N,Y,ids)
    A=sparse([a;b;a;b],[a;b;b;a],[G;G;-G;-G],N,N);
    A(ids,ids)=A(ids,ids)-Y;                    %#ok<SPRIX>
end

function [f,Phi,B,Gd]=resid(U,ME,BH,mu0,Y,ids,N)
%  Residu nodal et permeances INCREMENTALES.
    dU=U(ME.a)-U(ME.b)+ME.E;
    Fe=ME.iron; nb=numel(dU);
    Phi=zeros(nb,1); B=zeros(nb,1); Gd=zeros(nb,1);
    %  air : lineaire
    Phi(~Fe)=ME.P(~Fe).*dU(~Fe); Gd(~Fe)=ME.P(~Fe);
    %  fer : B = Bof(DeltaU/l), inversion directe de H(B)*l = DeltaU
    H=dU(Fe)./ME.l(Fe);
    Bf=BH.Bof(H);
    Phi(Fe)=ME.A(Fe).*Bf; B(Fe)=Bf;
    Gd(Fe)=ME.A(Fe)./(ME.l(Fe).*max(BH.dHdB(Bf),1e-9));
    m=~Fe & ME.A>0; B(m)=Phi(m)./ME.A(m);
    f=accumarray(ME.a,Phi,[N 1])-accumarray(ME.b,Phi,[N 1]);
    f(ids)=f(ids)-Y*U(ids)-ME.Isrc(ids);
end

function v=M_default(), v=3000; end
