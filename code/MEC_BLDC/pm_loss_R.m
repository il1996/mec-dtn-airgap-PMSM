function P = pm_loss_R(M, R, n_rpm)
%PM_LOSS_R  Pertes Foucault aimant a partir d'une solution DEJA calculee.
%   R doit provenir de cogging_mec ou cogging_mec2 avec span = 2*pi/Ns
%   (periodicite du motif dans le repere rotor).
%   Formulation 2D a courant net nul : P = sigma*L*INT (dA/dt-<dA/dt>)^2 dS.
%   Seule la part pilotee par les potentiels de surface intervient (la
%   magnetisation est statique dans le repere rotor).

if nargin<3||isempty(n_rpm), n_rpm=M.FEA.n_nl; end
om=n_rpm*2*pi/60;
AG=R.AG; nu=AG.nu; phis=R.phis; Np=numel(phis);
Usc=AG.Wc*R.Usurf; Uss=AG.Ws*R.Usurf;

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
% gradient(X,h) sur une MATRICE derive selon la DIM 2 : squeeze(A(j,:,:)) est
% (nth x Np), dim 2 = phi -> NE PAS transposer.
dAdt=zeros(size(A));
for j=1:nr, dAdt(j,:,:)=om*gradient(squeeze(A(j,:,:)),phis(2)-phis(1)); end

arc=M.embrace*2*pi/M.Nm; P=0;
for k=0:M.Nm-1
    m=abs(angle(exp(1i*(thr-k*2*pi/M.Nm))))<=arc/2;
    for q=1:Np
        w=zeros(1,nr);
        for j=1:nr
            d=dAdt(j,m,q); w(j)=sum((d-mean(d)).^2)*rq(j)*dth;
        end
        P=P+trapz(rq,w)/Np;
    end
end
P=M.sigma_pm*M.ls*P;
end
