function [P,D] = pm_loss(M, Nsurf, Np, kfringe, n_rpm, ext)
%
%   SOURCE DES POTENTIELS DE SURFACE (argument 'ext', tache T18).
%   Par defaut la chaine est pilotee par le reseau a UNE DENT (cogging_mec).
%   En passant ext = struct('AG',..,'Usurf',..,'phis',..) on la pilote par un
%   AUTRE modele -- typiquement le maillage polaire -- sans rien changer a la
%   formulation des pertes. C'est le seul point qui doit varier entre les
%   colonnes 'Lumped' et 'Mesh' du tableau : reecrire la chaine reintroduirait
%   les pieges qu'elle documente (division par nu, balayage sur un PAS
%   D'ENCOCHE, masque d'embrasse, derivee selon la bonne dimension).
%PM_LOSS  Pertes par courants de Foucault dans les aimants (C3).
%
%   Au synchronisme, les aimants tournent AVEC le fondamental : seule la
%   MODULATION DE DENTURE (f = Ns*n/60) cree un dB/dt dans l'aimant.
%
%   Formulation 2D rigoureuse (celle d'ANSYS 2D transitoire) : courants
%   induits axiaux a courant NET NUL par aimant,
%       J_z = -sigma*(dA/dt - <dA/dt>),   P = sigma*L*INT (dA/dt-<dA/dt>)^2 dS
%   avec A obtenu du champ radial par B_r = (1/r)*dA/dtheta. On n'utilise PAS
%   la formule forfaitaire (w^2/12rho)(dB/dt)^2 V : elle suppose B uniforme
%   sur la largeur d'aimant, faux ici (longueur d'onde de denture ~ arc
%   d'aimant). Seule la part pilotee par les potentiels de surface stator
%   intervient (la magnetisation est statique dans le repere rotor).
%
%   Balayage sur UN PAS D'ENCOCHE : le motif du repere rotor y est periodique.

if nargin<5||isempty(n_rpm), n_rpm=M.FEA.n_nl; end
om=n_rpm*2*pi/60;
if nargin>=6 && ~isempty(ext)
    AG=ext.AG; phis=ext.phis(:).'; Us=ext.Usurf; R=ext;
    Np=numel(phis);
else
    R=cogging_mec(M,Nsurf,0,Np,1500,kfringe,2*pi/M.Ns);
    AG=R.AG; phis=R.phis; Us=R.Usurf;
end
nu=AG.nu;
Usc=AG.Wc*Us; Uss=AG.Ws*Us;

nr=7; rq=linspace(AG.rmi,AG.Rro,nr);
nth=1440; thr=linspace(0,2*pi,nth+1); thr(end)=[]; dth=thr(2)-thr(1);

A=zeros(nr,nth,Np);
for q=1:Np
    pic=AG.alphau.*Usc(:,q); pis=AG.alphau.*Uss(:,q);
    for j=1:nr
        r=rq(j); v=log(r/AG.rmi);
        kk=-AG.mu0*AG.mur*(nu/r).*cosh(nu*v)./sinh(nu*AG.Um);
        th=thr+phis(q);
        A(j,:,q)=r*( ((kk.*pic)./nu).'*sin(nu*th) - ((kk.*pis)./nu).'*cos(nu*th) );
    end
end
% ATTENTION : gradient(X,h) sur une MATRICE derive selon la DIM 2. Ici
% squeeze(A(j,:,:)) est (nth x Np) : la dim 2 est bien phi -> NE PAS
% transposer (une transposition ferait deriver selon theta, avec en plus le
% mauvais pas : erreur d'un facteur ~1.7 sur la perte).
dAdt=zeros(size(A));
for j=1:nr, dAdt(j,:,:)=om*gradient(squeeze(A(j,:,:)),phis(2)-phis(1)); end

%  COURANT NET NUL SUR LA SECTION ENTIERE. La condition physique est
%  INT_S J dS = 0 sur toute la section de l'aimant, soit UNE constante par
%  aimant et par instant. L'imposer rayon par rayon reviendrait a supposer
%  chaque couronne elementaire isolee des autres et retirerait aussi la
%  variation RADIALE de dA/dt -- forte ici, le champ de denture decroissant
%  en cosh(n ln(r/rmi)) a travers l'aimant.
wr=zeros(nr,1); wr(1)=(rq(2)-rq(1))/2; wr(end)=(rq(end)-rq(end-1))/2;
if nr>2, wr(2:end-1)=(rq(3:end)-rq(1:end-2))/2; end
wS=wr(:).*rq(:)*dth;
arc=M.embrace*2*pi/M.Nm; P=0;
for k=0:M.Nm-1
    m=abs(angle(exp(1i*(thr-k*2*pi/M.Nm))))<=arc/2;
    Sk=sum(wS)*sum(m);
    for q=1:Np
        num=0;
        for j=1:nr, num=num+wS(j)*sum(dAdt(j,m,q)); end
        K=num/Sk;
        for j=1:nr, P=P+wS(j)*sum((dAdt(j,m,q)-K).^2)/Np; end
    end
end
P=M.sigma_pm*M.ls*P;
D.fslot=M.Ns*n_rpm/60;
D.delta=sqrt(2/(2*pi*D.fslot*4*pi*1e-7*M.sigma_pm));
D.att8=(AG.rmi/AG.Rro)^8; D.att22=(AG.rmi/AG.Rro)^22;
D.R=R;
end
