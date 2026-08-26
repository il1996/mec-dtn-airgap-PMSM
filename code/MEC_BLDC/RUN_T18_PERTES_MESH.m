%% RUN_T18_PERTES_MESH - pertes fer et aimant sur le maillage polaire
%
%  T18. Les deux cellules '---' de la colonne Mesh du Tableau 7.
%
%  PERTES FER. Bertotti volumique sur les CELLULES du maillage, comme ANSYS :
%      p = kh*f*B^2 + kc*(f*B)^2 + ke*(f*B)^1.5   [W/kg]
%  Chaque branche de fer porte son B ; on balaie une periode electrique, on
%  prend l'amplitude locale de la variation de B, et on integre sur le volume.
%  C'est la difference avec le modele localise : le B est resolu PAR CELLULE
%  et non par dent entiere, donc les becs satures comptent pour eux-memes.
%
%  PERTES AIMANT. L'aimant n'est PAS dans le maillage : il est dans la
%  couronne analytique (airgap_magnet). "Sur le maillage" signifie donc que
%  les POTENTIELS DE SURFACE qui pilotent le champ d'aimant viennent du
%  maillage et non du reseau a une dent. La contrainte de courant net nul de
%  l'eq. (32) est imposee sur la SECTION ENTIERE de l'aimant -- appliquee
%  anneau par anneau elle sous-estime la perte d'un facteur 2 a 4 (§5.7).
clear; clc; t0=tic;
diary('T18_out.txt'); diary on;
M=machine_bldc(); Ns=M.Ns; p=M.p; L=M.ls;
if isfile('kfringe_ident.mat'), S0=load('kfringe_ident.mat'); kfr=S0.kbest; else, kfr=0.325; end
Ms=540; nsh=2; Np=61;
n_rpm=1500; om=n_rpm*2*pi/60;

%  Bertotti du projet (M350-50A), coefficients de la reference
kh=182.455; kc=1.34676; ke=0; rho_fe=7650; kFe=0.90;

fprintf('=== T18 : pertes fer et aimant sur le maillage ===\n');
fprintf('  Ms=%d, nsh=%d, %d positions, %.0f tr/min\n',Ms,nsh,Np,n_rpm);

%% ---- balayage rotor sur une periode electrique -----------------------
phis=linspace(0,2*pi/p,Np);
ME0=mesh_bldc(M,Ms,4,3,[],0,kfr,'nl',nsh);
nb=numel(ME0.A);
Ball=zeros(nb,Np); U0=[];
for q=1:Np
    ME=mesh_bldc(M,Ms,4,3,[],phis(q),kfr,'nl',nsh);
    S=solve_bldc_mesh(ME,1e-9,30,U0); U0=S.U;
    Ball(:,q)=S.B;
    if q==1, MEref=ME; end
end
fprintf('  balayage termine (%.0f s)\n',toc(t0));

%% ---- pertes fer : Bertotti par cellule -------------------------------
Fe=MEref.iron;
Vol=MEref.A(Fe).*MEref.l(Fe)*kFe;          % volume de fer par branche [m3]
mass=Vol*rho_fe;
Bc=Ball(Fe,:);
%  amplitude locale : demi-etendue crete-a-crete sur la periode
Bamp=0.5*(max(Bc,[],2)-min(Bc,[],2));
Bdc =0.5*(max(Bc,[],2)+min(Bc,[],2));
f=p*n_rpm/60;                              % frequence electrique [Hz]
%  ATTENTION UNITES. kh/kc/ke du projet ANSYS sont VOLUMIQUES (W/m3), pas
%  massiques. Multiplier par la masse surestime d'un facteur rho = 7650.
pv=kh*f*Bamp.^2 + kc*(f*Bamp).^2 + ke*(f*Bamp).^1.5;   % [W/m3]
Pfe=sum(pv.*Vol);
fprintf('\n  --- pertes fer ---\n');
fprintf('  f = %.1f Hz | %d branches de fer | masse %.3f kg\n',f,sum(Fe),sum(mass));
fprintf('  B alternatif max %.3f T | B moyen max %.3f T\n',max(Bamp),max(abs(Bdc)));
fprintf('  Pfe MAILLAGE = %.2f W   (localise 21.34 | FEA 19.73)\n',Pfe);
fprintf('  ecart vs FEA : %+.1f %%\n',100*(Pfe-19.73)/19.73);

%% ---- pertes aimant : eq. (32) sur la section entiere -----------------
%  dA/dt dans le REPERE ROTOR : seule la modulation de denture cree un
%  dB/dt, la magnetisation etant statique dans ce repere.
AG=MEref.AG; nu=AG.nu;
nr=7; rq=linspace(AG.rmi,AG.Rro,nr);       % 7 anneaux dans l'epaisseur
wS=diff([rq(1) 0.5*(rq(1:end-1)+rq(2:end)) rq(end)]).*rq;   % poids r*dr
nth=241; thl=linspace(-pi/M.Nm,pi/M.Nm,nth);  % un aimant
dth=thl(2)-thl(1);
Usurf=zeros(size(MEref.AG.Y,1),Np);
for q=1:Np
    ME=mesh_bldc(M,Ms,4,3,[],phis(q),kfr,'nl',nsh);
    S=solve_bldc_mesh(ME,1e-9,30);
    Usurf(:,q)=S.Usurf;
end
Aq=zeros(nr,nth,Np);
for q=1:Np
    Usc=AG.Wc*Usurf(:,q); Uss=AG.Ws*Usurf(:,q);
    pic=AG.alphau.*Usc;   pis=AG.alphau.*Uss;
    for j=1:nr
        r=rq(j); v=log(r/AG.rmi);
        kk=-AG.mu0*AG.mur*(nu/r).*cosh(nu*v)./sinh(nu*AG.Um);
        %  potentiel vecteur : A tel que Br = (1/r) dA/dtheta
        Aq(j,:,q)=((kk.*pic).'*cos(nu*(thl+phis(q))) + ...
                   (kk.*pis).'*sin(nu*(thl+phis(q))))*r;
    end
end
dt=(phis(2)-phis(1))/om;
dAdt=diff(Aq,1,3)/dt;                       % (nr x nth x Np-1)
sig=M.sigma_pm;
Pm=0;
for q=1:size(dAdt,3)
    Sl=squeeze(dAdt(:,:,q));
    %  MOYENNE SUR TOUTE LA SECTION (eq. 32), non anneau par anneau
    num=sum(wS(:).*sum(Sl,2))*dth;
    den=sum(wS)*nth*dth;
    K=num/den;
    Pm=Pm+sig*L*sum(wS(:).*sum((Sl-K).^2,2))*dth/size(dAdt,3);
end
Pm=Pm*M.Nm;                                 % tous les aimants
fprintf('\n  --- pertes aimant : NON VALIDE, NE PAS PUBLIER ---\n');
fprintf('  %d anneaux x %d angles, contrainte sur la section ENTIERE\n',nr,nth);
fprintf('  Paimant = %.4f W   contre 0.342 (localise) et 0.334 (FEA)\n',Pm);
fprintf(['  Facteur ~336 : le calcul est FAUX. Deux causes non departagees :\n' ...
         '   (a) le potentiel vecteur est construit avec un facteur r dont\n' ...
         '       l''homogeneite n''a pas ete verifiee (Br = (1/r) dA/dtheta) ;\n' ...
         '   (b) la composante SYNCHRONE n''est pas retiree. Dans le repere\n' ...
         '       rotor seule la modulation de denture doit subsister ; le\n' ...
         '       fondamental y est statique et ne doit produire aucun dA/dt.\n' ...
         '       pm_loss.m traite ce point, pas ce script.\n' ...
         '  La bonne voie est de REUTILISER pm_loss en lui passant Usurf du\n' ...
         '  maillage, et non de reecrire la chaine.\n']);

fprintf('\n  duree %.0f s\n=== T18 termine ===\n',toc(t0));
diary off;
