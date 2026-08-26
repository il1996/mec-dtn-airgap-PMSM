%% RUN_DIAG_FEM - la FEM de la FEA en charge est-elle celle de l'essai a vide ?
%  Le modele MEC donne un courant 30 % plus faible que la FEA au meme point.
%  Comme R, Vdc et L concordent, la seule variable restante est la FEM.
%  On l'extrait des sondes FEA : e_ab = (Ie_a - Ie_b) - L_ll di/dt.
clear; clc;
M=machine_bldc(); p=M.p; fea=M.FEA.dir; ol=fullfile(fea,'transitoire (en charge)');
rd=@(f)readmatrix(f,'FileType','text','Delimiter','\t','NumHeaderLines',1);
dV=rd(fullfile(ol,'NodeVoltage Plot 2.tab')); dI=rd(fullfile(ol,'BranchCurrent Plot 1.tab'));
dS=rd(fullfile(ol,'Speed Plot 1.tab')); dL=rd(fullfile(ol,'Plot 1_loss.tab'));
dFLo=rd(fullfile(ol,'Output Variables Plot 1.tab'));
t=dV(:,1)*1e-9; Ie=dV(:,2:4)/1e3; i=interp1(dI(:,1)*1e-3,dI(:,2:4),t,'linear','extrap');
n=interp1(dS(:,1)*1e-3,dS(:,2),t,'linear','extrap'); om=n*2*pi/60;
di=[gradient(i(:,1),t) gradient(i(:,2),t) gradient(i(:,3),t)];

RI=inductance_mec(M,1260,M.muI,0); Lll=2*RI.Ld;
m=t>0.75*t(end) & abs(i(:,3))<1e-6 & abs(i(:,1))>0.3;   % conduction a/b pure
eab=(Ie(:,1)-Ie(:,2))-Lll*(di(:,1));                    % di_a = -di_b
fprintf('=== FEM de ligne en charge, extraite des sondes FEA ===\n');
fprintf('  %d points de conduction a/b pure sur %d\n',sum(m),numel(t));
fprintf('  tension aux bornes (Ie_a - Ie_b)   : %7.1f V (moyenne)\n',mean(Ie(m,1)-Ie(m,2)));
fprintf('  chute inductive L_ll di/dt         : %7.1f V\n',mean(Lll*di(m,1)));
fprintf('  -> FEM e_ab en charge              : %7.1f V a %.0f tr/min\n',mean(eab(m)),mean(n(m)));
fprintf('     soit psi_ab = %.3f Wb/rad\n',mean(eab(m))/mean(om(m)));

%% ---- reference : FEM A VIDE ramenee a la meme vitesse ------------------
dE=rd(fullfile(fea,'transitoire (Back_emf)','Output Variables Plot.tab'));
env=max(dE(:,3)); nnl=M.FEA.n_nl;
fprintf('\n  enveloppe a vide FEA : %.1f V a %d tr/min\n',env,nnl);
fprintf('  ramenee a %.0f tr/min : %.1f V   -> psi_ab = %.3f Wb/rad\n', ...
    mean(n(m)),env*mean(n(m))/nnl,env/(nnl*2*pi/60));
fprintf('  => chute de FEM en charge : %+.1f %%\n', ...
    100*(mean(eab(m))/(env*mean(n(m))/nnl)-1));

%% ---- meme extraction sur le FLUX TOTALISE ------------------------------
FL=[dFLo(:,2) dFLo(:,4) dFLo(:,3)];                     % a, b, c
tf=dFLo(:,1)*1e-3; FLi=interp1(tf,FL,t,'linear','extrap');
lam_pm=FLi-[i(:,1) i(:,2) i(:,3)]*RI.Ld;                % retrait de L*i
fprintf('\n  flux totalise a vide reconstitue (FL - L_d*i) :\n');
fprintf('     crete phase a en charge = %.4f Wb\n',max(abs(lam_pm(t>0.5*t(end),1))));
d4=rd(fullfile(fea,'transitoire (Back_emf)','Output Variables Plot 4.tab'));
fprintf('     crete phase a a vide    = %.4f Wb\n',max(abs(d4(:,3))));

%% ---- verification de la coherence du couple ---------------------------
Pw=sum(Ie.*i,2); Pconv=Pw-M.Rph*0;                      % Ie est deja apres Rph
Pem=interp1(dL(:,1)*1e-9,dL(:,7),t,'linear','extrap');
mm=t>0.75*t(end);
fprintf('\n  puissance dans les enroulements  : %7.1f W\n',mean(Pw(mm)));
fprintf('  Pem declare par ANSYS            : %7.1f W\n',mean(Pem(mm)));
fprintf('  couple deduit des enroulements   : %7.3f N.m\n',mean(Pw(mm))/mean(om(mm)));
fprintf('  couple declare par ANSYS         : %7.3f N.m\n',mean(Pem(mm))/mean(om(mm)));
