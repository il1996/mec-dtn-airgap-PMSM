%% RUN_PMLOSS_NL - optimisation des pertes aimant A VIDE
%  Reference FEA : SolidLoss = 0.3335 W a 1688 tr/min (essai 'transitoire
%  (a vide)'). Le calcul actuel (pm_loss) donne 0.0818 W, soit -75.5 %.
%  A vide, les aimants tournent AVEC le fondamental : seule la modulation de
%  DENTURE cree un dB/dt dans l'aimant. Le motif du repere rotor est donc
%  periodique sur UN PAS D'ENCOCHE.
%
%  On reprend le calcul avec la chaine amelioree mise au point pour la
%  charge (derivee de rotation ANALYTIQUE, contrainte de courant net nul sur
%  la SECTION ENTIERE, quadrature radiale en cosinus), puis on etablit la
%  convergence en Nsurf / Np / nr avant de conclure.
clear; clc; t0=tic;
M=machine_bldc(); Ns=M.Ns; n=M.FEA.n_nl; om=n*2*pi/60;
if isfile('kfringe_ident.mat'), S0=load('kfringe_ident.mat'); kfr=S0.kbest; else, kfr=0.75; end
PF=M.FEA.Psolid_nl;
fprintf('=== Pertes aimant A VIDE : optimisation ===\n');
fprintf('  reference FEA : %.4f W a %d tr/min (denture a %.0f Hz)\n',PF,n,Ns*n/60);
P0=pm_loss(M,840,91,kfr,n);
fprintf('  calcul actuel (pm_loss, Nsurf=840, Np=91, 9 rayons) : %.4f W (%+.1f %%)\n', ...
    P0,100*(P0-PF)/PF);

%  la chaine amelioree s'applique telle quelle avec des courants NULS
Tslot=2*pi/(Ns*om);                                  % duree d'un pas d'encoche
run1=@(Nsurf,Np,Nt,nr) local_run(M,Nsurf,Np,Nt,nr,kfr,om,Tslot);

fprintf('\n  --- (1) convergence en resolution de surface (Np=181, Nt=721, nr=25) ---\n');
fprintf('  %8s %10s %10s\n','Nsurf','P (W)','ecart');
for Nsurf=[840 1260 1680 2520]
    P=run1(Nsurf,181,721,25);
    fprintf('  %8d %10.4f %8.1f %%\n',Nsurf,P,100*(P-PF)/PF);
end

fprintf('\n  --- (2) convergence en positions rotor et en temps (Nsurf=1680) ---\n');
fprintf('  %8s %8s %10s %10s\n','Np','Nt','P (W)','ecart');
for k=[91 181 361 721]
    P=run1(1680,k,min(4*k,1441),25);
    fprintf('  %8d %8d %10.4f %8.1f %%\n',k,min(4*k,1441),P,100*(P-PF)/PF);
end

fprintf('\n  --- (3) quadrature radiale (Nsurf=1680, Np=361) ---\n');
fprintf('  %8s %10s %10s\n','nr','P (W)','ecart');
for nr=[9 17 25 41 61]
    P=run1(1680,361,1441,nr);
    fprintf('  %8d %10.4f %8.1f %%\n',nr,P,100*(P-PF)/PF);
end

%% ---- attribution : deficit de la bande de denture ---------------------
Pb=run1(1680,361,1441,41);
R=cogging_mec(M,1680,0,181,M.muI,kfr,2*pi/Ns);
fea=M.FEA.dir;
d4=readmatrix(fullfile(fea,'magnetostique(Magnetic_loading)', ...
    'Calculator Expressions Plot 4.tab'),'FileType','text','NumHeaderLines',1);
thu=linspace(0,2*pi,3601); thu(end)=[];
BrFu=interp1([d4(:,2)*pi/180;2*pi],[d4(:,3);d4(1,3)],thu,'linear','extrap');
amp=@(B,th,k) 2*abs(mean(B(:).'.*exp(-1i*k*th(:).')));
thm=R.thq(1:end-1); Brm=R.Br(1:end-1);
r8=amp(BrFu,thu,8)/amp(Brm,thm,8);
r22=amp(BrFu,thu,22)/amp(Brm,thm,22);
fprintf('\n  --- attribution de l''ecart residuel ---\n');
fprintf('  valeur convergee : %.4f W (%+.1f %%)\n',Pb,100*(Pb-PF)/PF);
fprintf('  facteur d''amplitude qu''il faudrait sur B : %.2f (perte en B^2)\n',sqrt(PF/Pb));
fprintf('  bandes de denture mesurees FEA/MEC : rang 8 -> %.2f, rang 22 -> %.2f\n',r8,r22);
fprintf('  attenuation radiale (rmi/Rro)^nu : rang 8 %.3f, rang 22 %.3f\n', ...
    (R.AG.rmi/R.AG.Rro)^8,(R.AG.rmi/R.AG.Rro)^22);
fprintf('  ENCADREMENT : %.4f W (brut) < %.4f W (FEA) < %.4f W (corrige du rang 8)\n', ...
    Pb,PF,Pb*r8^2);
fprintf('\n  (%.0f s)\n',toc(t0));

%% ------------------------------------------------------------------------
function P = local_run(M,Nsurf,Np,Nt,nr,kfr,om,Tslot)
%  Pertes aimant a vide par la chaine amelioree : on reutilise pm_loss_load
%  avec des COURANTS NULS. Le motif etant periodique sur un pas d'encoche,
%  le balayage rotor et la fenetre temporelle couvrent exactement un pas.
    R=cogging_mec(M,Nsurf,0,Np,M.muI,kfr,2*pi/M.Ns);
    ts=linspace(0,Tslot,Nt);
    P=pm_loss_load(M,R,zeros(Nsurf,3),ts,zeros(3,Nt),om*ts,nr);
end
