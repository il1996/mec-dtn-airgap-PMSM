%% RUN_MAP - construction et controle de la cartographie non lineaire MEC
clear; clc;
M=machine_bldc();
%  ROLE (b) DE LA FRANGE : pont TANGENTIEL dent-a-dent. La cartographie
%  sert le modele d'entrainement ; elle garde donc la valeur historique
%  0.75, et NON la valeur 0.325 identifiee sur les bandes de denture, qui
%  releve du role (a) -- l'entree RADIALE du flux par la bouche d'encoche.
%  Verifie : appliquer 0.325 ici degrade le couple par ampere sature de
%  31 % et l'ensemble du balayage en vitesse.
kfr=0.75;
t0=tic;
%  Grille de courant CONCENTREE PRES DE ZERO : le regime etabli travaille a
%  2 A (il faut y resoudre finement dL/di) alors que le demarrage monte a
%  23 A (il faut couvrir la plage). Un maillage quadratique fait les deux.
u=linspace(-1,1,19); iax=28*sign(u).*abs(u).^2;
S=mec_map(M,720,kfr,61,iax);
fprintf('=== Cartographie MEC non lineaire ===\n');
fprintf('  grille : %d positions x %d x %d courants (%d iterations, %.1f s)\n', ...
    numel(S.the),numel(S.iax),numel(S.iax),S.nit,toc(t0));
save('mec_map.mat','S','-v7.3');

%% --- controle 1 : a courant nul on doit retrouver le flux a vide --------
i0=find(abs(S.iax)<1e-9,1);
lam0=squeeze(S.lam(1,:,i0,i0));
fprintf('\n  a i = 0 : lambda_a crete = %.4f Wb\n',max(abs(lam0)));
fprintf('            FEM crete a 1500 tr/min = %.1f V\n', ...
    max(abs(squeeze(S.psi(1,:,i0,i0))))*M.speed*2*pi/60);

%% --- controle 2 : effondrement du flux avec le courant ------------------
fprintf('\n  --- effondrement de la FEM sous saturation ---\n');
fprintf('  %10s %14s %14s %12s\n','i_ligne(A)','lambda_ab c-c','psi_ab max','L_inc (mH)');
for ii=[0 2 5 10 15 20 25]
    ia=ii; ib=-ii/sqrt(3);                    % motif de conduction a+ / b-
    la=interpn(1:3,S.the,S.iax,S.iax,S.lam,1,S.the,ia,ib,'linear');
    lb=interpn(1:3,S.the,S.iax,S.iax,S.lam,2,S.the,ia,ib,'linear');
    pa=interpn(1:3,S.the,S.iax,S.iax,S.psi,1,S.the,ia,ib,'linear');
    pb=interpn(1:3,S.the,S.iax,S.iax,S.psi,2,S.the,ia,ib,'linear');
    Laa=interpn(1:3,S.the,S.iax,S.iax,S.La,1,S.the,ia,ib,'linear');
    Lba=interpn(1:3,S.the,S.iax,S.iax,S.La,2,S.the,ia,ib,'linear');
    Lab=interpn(1:3,S.the,S.iax,S.iax,S.Lb,1,S.the,ia,ib,'linear');
    Lbb=interpn(1:3,S.the,S.iax,S.iax,S.Lb,2,S.the,ia,ib,'linear');
    Leq=mean((Laa-Lba)-(Lab-Lbb)/sqrt(3));
    fprintf('  %10.1f %14.4f %14.4f %12.2f\n', ...
        ii,max(la-lb)-min(la-lb),max(abs(pa-pb)),Leq*1e3);
end

%% --- controle 3 : couple ------------------------------------------------
fprintf('\n  --- couple moyen (MST) selon le courant de ligne ---\n');
fprintf('  %10s %14s %14s\n','i_ligne(A)','T moyen (N.m)','T/i (N.m/A)');
for ii=[1 2 5 10 15 20 25]
    ia=ii; ib=-ii/sqrt(3);
    Tq=interpn(S.the,S.iax,S.iax,S.T,S.the,ia,ib,'linear');
    fprintf('  %10.1f %14.3f %14.4f\n',ii,mean(abs(Tq)),mean(abs(Tq))/ii);
end
fprintf('\n-> mec_map.mat (%.1f s)\n',toc(t0));
