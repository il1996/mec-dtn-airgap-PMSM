%% RUN_PMLOSS  Pertes par courants de Foucault dans les AIMANTS (C3)
%
%  PHYSIQUE. Au synchronisme et a vide, les aimants tournent AVEC le
%  fondamental : ils n'en voient aucune variation. Seule la MODULATION DE
%  DENTURE (passage des dents, f = Ns*n/60 = 422 Hz ici) cree un dB/dt dans
%  l'aimant. C'est exactement le champ que le modele couple sait produire.
%
%  FORMULATION 2D RIGOUREUSE (celle d'ANSYS en 2D transitoire) : les courants
%  induits sont axiaux, avec courant NET NUL dans chaque aimant :
%      J_z = -sigma*(dA/dt - <dA/dt>)        <.> = moyenne sur la section
%      P   = sigma*L*INT (dA/dt - <dA/dt>)^2 dS
%  Le potentiel vecteur s'obtient du champ radial par B_r = (1/r) dA/dtheta.
%  On N'UTILISE PAS la formule forfaitaire p = (w^2/12rho)(dB/dt)^2 V : elle
%  suppose B UNIFORME sur la largeur d'aimant, ce qui est faux ici (la
%  longueur d'onde de denture ~ un pas d'encoche 24 deg est du meme ordre que
%  l'arc d'aimant 24.2 deg).
%
%  Seule la part du champ pilotee par les POTENTIELS DE SURFACE STATOR
%  intervient : la part de magnetisation est STATIQUE dans le repere rotor
%  (elle tourne avec l'aimant) et ne produit aucune perte.
%
%  Validite : epaisseur de peau dans le NdFeB a 422 Hz = 29 mm >> dimensions
%  de l'aimant (3.2 x 12.6 mm) -> regime RESISTIF, pas de reaction d'induit
%  des courants de Foucault (approximation basse frequence justifiee).
%
%  Reference FEA : SolidLoss = 0.333 W a vide (1688 tr/min).
clear; clc;
M=machine_bldc();
if isfile('kfringe_ident.mat'), S=load('kfringe_ident.mat'); kfr=S.kbest; else, kfr=0.75; end
Ns=M.Ns; Nm=M.Nm; n=M.FEA.n_nl;
om=n*2*pi/60;  fslot=Ns*n/60;
fprintf('=== PERTES AIMANT PAR COURANTS DE FOUCAULT (C3) ===\n');
fprintf('  %g tr/min -> frequence de passage de dents = %.0f Hz\n',n,fslot);
delta=sqrt(2/(2*pi*fslot*4*pi*1e-7*M.sigma_pm));
fprintf('  Epaisseur de peau NdFeB (rho=%.1f uOhm.m) = %.1f mm >> aimant %.1f mm -> regime resistif\n',...
        M.rho_pm*1e6,delta*1e3,M.hm*1e3);

%% ---- balayage sur UN PAS D'ENCOCHE (periodicite du repere rotor) -------
Np=181;
R=cogging_mec(M,840,0,Np,1500,kfr,2*pi/Ns);
AG=R.AG; nu=AG.nu; phis=R.phis;

%% ---- champ AC dans l'aimant, repere ROTOR ------------------------------
%  phi_i(Us) = alphau .* Us_nu   (la part source est statique -> exclue)
Usc=AG.Wc*R.Usurf; Uss=AG.Ws*R.Usurf;                 % (numax x Np)
nr=9; rq=linspace(AG.rmi,AG.Rro,nr);                  % quadrature radiale
nth=2160; thr=linspace(0,2*pi,nth+1); thr(end)=[];    % angles ROTOR
dth=thr(2)-thr(1);

% A(theta_rotor, phi) par rayon : B_r = (1/r) dA/dtheta -> A = r*INT B_r dtheta
A=zeros(nr,nth,Np);
for q=1:Np
    pic=AG.alphau.*Usc(:,q); pis=AG.alphau.*Uss(:,q);
    for j=1:nr
        r=rq(j); v=log(r/AG.rmi);
        kk=-AG.mu0*AG.mur*(nu/r).*cosh(nu*v)./sinh(nu*AG.Um);
        % champ evalue aux angles STATOR theta = theta_rotor + phi
        th=thr+phis(q);
        bc=kk.*pic; bs=kk.*pis;
        % A = r * INT B_r dtheta  (par harmonique : cos->sin/nu, sin->-cos/nu)
        A(j,:,q)=r*( (bc./nu).'*sin(nu*th) - (bs./nu).'*cos(nu*th) );
    end
end

%% ---- dA/dt dans le repere rotor ---------------------------------------
% ATTENTION : gradient(X,h) sur une MATRICE derive selon la DIM 2. Ici
% squeeze(A(j,:,:)) est (nth x Np), dim 2 = phi -> NE PAS transposer.
dAdphi=zeros(size(A));
for j=1:nr, dAdphi(j,:,:)=gradient(squeeze(A(j,:,:)),phis(2)-phis(1)); end
dAdt=om*dAdphi;

%% ---- integration par aimant -------------------------------------------
arc=M.embrace*2*pi/Nm;                                % arc d'un aimant
P=0;
for k=0:Nm-1
    c=k*2*pi/Nm;
    m=abs(angle(exp(1i*(thr-c))))<=arc/2;             % secteur de l'aimant k
    for q=1:Np
        w=zeros(1,nr);
        for j=1:nr
            d=dAdt(j,m,q);
            dm=sum(d.*rq(j)*dth)/sum(rq(j)*dth*ones(size(d)));   % <dA/dt> pondere
            w(j)=sum((d-dm).^2 .* rq(j)*dth);          % INT (.)^2 r dtheta
        end
        P=P+trapz(rq,w)/Np;                            % moyenne temporelle
    end
end
P=M.sigma_pm*M.ls*P;

fprintf('\n[1] Geometrie : %d aimants, arc %.2f deg, epaisseur %.2f mm\n',...
        Nm,arc*180/pi,M.hm*1e3);
fprintf('    Volume total d''aimant : %.2f cm3\n',Nm*arc*(AG.Rro^2-AG.rmi^2)/2*M.ls*1e6);

fprintf('\n[2] RESULTAT\n');
fprintf('    Pertes aimant MEC : %.4f W\n',P);
fprintf('    Pertes aimant FEA : %.4f W  (SolidLoss a vide)\n',M.FEA.Psolid_nl);
fprintf('    ECART : %+.1f %%\n',(P-M.FEA.Psolid_nl)/M.FEA.Psolid_nl*100);

fprintf('\n[3] Sensibilite a la resistivite (donnee materiau incertaine)\n');
for rho=[1.2 1.4 1.6]*1e-6
    fprintf('    rho = %.1f uOhm.m -> %.4f W (%+.1f %%)\n',rho*1e6,P*M.rho_pm/rho,...
            (P*M.rho_pm/rho-M.FEA.Psolid_nl)/M.FEA.Psolid_nl*100);
end

%% ---- [4] ATTRIBUTION DE L'ECART -------------------------------------
%  La perte varie en B^2 et n'est pilotee QUE par les harmoniques de denture
%  (ceux qui ne tournent pas avec le rotor). L'harmonique dominant vu par
%  l'aimant est nu = |Ns-p| = 8 (le moins attenue radialement : (rmi/Rro)^nu).
fea=M.FEA.dir;
d4=readmatrix(fullfile(fea,'magnetostique(Magnetic_loading)','Calculator Expressions Plot 4.tab'),'FileType','text','NumHeaderLines',1);
angF=d4(:,2)*pi/180; BrF=d4(:,3);
thu=linspace(0,2*pi,3601); thu(end)=[];
BrFu=interp1([angF;2*pi],[BrF;BrF(1)],thu,'linear','extrap');
amp=@(B,th,k) 2*abs(mean(B.*exp(-1i*k*th)));
thm=R.thq(1:end-1); Brm=R.Br(1:end-1);
a8m=amp(Brm,thm,8);   a8f=amp(BrFu,thu,8);
a22m=amp(Brm,thm,22); a22f=amp(BrFu,thu,22);
fprintf('\n[4] ATTRIBUTION DE L''ECART (la perte varie en B^2)\n');
fprintf('    Attenuation radiale a travers l''aimant : (rmi/Rro)^nu\n');
fprintf('      nu = 8  -> %.3f  (harmonique DOMINANT dans l''aimant)\n',(AG.rmi/AG.Rro)^8);
fprintf('      nu = 22 -> %.3f  (fortement attenue)\n',(AG.rmi/AG.Rro)^22);
fprintf('    Harmonique 8 du champ d''entrefer : MEC %.4f T | FEA %.4f T (%+.1f %%)\n',...
        a8m,a8f,(a8m-a8f)/a8f*100);
fprintf('    Facteur attendu sur la perte : (%.4f/%.4f)^2 = %.2f\n',a8f,a8m,(a8f/a8m)^2);
fprintf('    -> perte corrigee de ce seul deficit : %.4f W  (%+.1f %% vs FEA)\n',...
        P*(a8f/a8m)^2,(P*(a8f/a8m)^2-M.FEA.Psolid_nl)/M.FEA.Psolid_nl*100);
fprintf('\n    BILAN HONNETE : le deficit connu de l''harmonique 8 explique environ la\n');
fprintf('    MOITIE de l''ecart (%.0f %% -> %.0f %%). Le reste tient a deux incertitudes\n',...
        (P-M.FEA.Psolid_nl)/M.FEA.Psolid_nl*100,(P*(a8f/a8m)^2-M.FEA.Psolid_nl)/M.FEA.Psolid_nl*100);
fprintf('    de la REFERENCE, non levables depuis les .tab :\n');
fprintf('      (a) la resistivite des aimants reglee dans ANSYS est inconnue ; la\n');
fprintf('          plage usuelle du NdFeB (0.6-1.6 uOhm.m) couvre un facteur >2 ;\n');
fprintf('      (b) "SolidLoss" agrege TOUS les conducteurs massifs : si la culasse\n');
fprintf('          rotorique est modelisee en acier massif conducteur, ses propres\n');
fprintf('          courants induits y sont comptes, et pas seulement les aimants.\n');
fprintf('    Ce qui EST etabli : la chaine de calcul (formulation 2D a courant net\n');
fprintf('    nul, champ de denture dans l''aimant) donne le bon ORDRE DE GRANDEUR et\n');
fprintf('    la bonne PHYSIQUE (regime resistif, pilotage par la denture a 422 Hz).\n');

figure('Name','Pertes aimant','Color','w','Position',[70 70 1000 400]);
subplot(1,2,1);
jm=ceil(nr/2);
plot(thr*180/pi,squeeze(dAdt(jm,:,1)),'b','LineWidth',1.2); grid on; xlim([0 60]);
xlabel('\theta rotor (deg)'); ylabel('dA/dt (V/m)');
title('dA/dt au rayon moyen d''aimant (modulation de denture)');
subplot(1,2,2);
bar([P M.FEA.Psolid_nl]); grid on; set(gca,'XTickLabel',{'MEC','FEA'});
ylabel('pertes aimant (W)');
title(sprintf('%.3f vs %.3f W (%+.0f %%)',P,M.FEA.Psolid_nl,...
      (P-M.FEA.Psolid_nl)/M.FEA.Psolid_nl*100));
saveas(gcf,fullfile(fileparts(mfilename('fullpath')),'FIG_pmloss.png'));
