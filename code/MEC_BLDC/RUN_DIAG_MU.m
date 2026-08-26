%% RUN_DIAG_MU  -  D'ou vient le deficit de denture ? Sensibilite a mu du fer
%  Le sous-domaine exact d'ouverture (RUN_SUBDOMAIN) donne les MEMES
%  harmoniques que le modele forfaitaire : la cause n'est pas l'ouverture.
%  Ici on balaie la permeabilite du fer dans le modele EXACT, et on situe
%  l'induction reelle de dent sur la courbe M350 (foisonnement inclus).
clear; clc;
M=machine_bldc(); BH=bh_curve('M350');

D=readmatrix(fullfile(M.FEA.dir,'magnetostique(Magnetic_loading)', ...
    'Calculator Expressions Plot 4.tab'),'FileType','text','Delimiter','\t','NumHeaderLines',1);
thF=D(:,2)*pi/180; BrF=D(:,3); [thF,ix]=sort(thF); BrF=BrF(ix);
fu=@(y,t,k)hypot(2*mean(y(:).*cos(k*t(:))),2*mean(y(:).*sin(k*t(:))));
Fa8=fu(BrF,thF,8); Fa22=fu(BrF,thF,22); Fa23=fu(BrF,thF,23); FBg1=fu(BrF,thF,M.p);

K=8; nm=ceil(4*K*pi/(M.ws0/M.Rsi));
fprintf('=== Sensibilite du spectre de denture a mu_r du fer (modele EXACT) ===\n');
fprintf('%8s | %8s %8s %8s %8s | %7s %9s\n','mu_r','Bg1','a8','a22','a23','a22/a8','rcond');
for mu=[200 500 1000 3000 10000 1e5 1e7]
    R=subdomain_mec(M,K,nm,mu,0);
    fprintf('%8g | %8.4f %8.5f %8.5f %8.5f | %7.3f %9.1e\n',mu,R.Bg1,R.a8,R.a22,R.a23,R.a22/R.a8,R.rcond);
end
fprintf('%8s | %8.4f %8.5f %8.5f %8.5f | %7.3f\n','FEA',FBg1,Fa8,Fa22,Fa23,Fa22/Fa8);

% --- ou se situe la dent sur la courbe M350 ? ---
R=subdomain_mec(M,K,nm,3000,0);
taus=2*pi/M.Ns; pitch=taus*M.Rsi;
PhiT_max=max(abs(R.PhiT));
Bt_app = PhiT_max/(M.wst1*M.ls);          % induction APPARENTE de dent
Bt_fe  = Bt_app/M.Ki;                      % induction dans le FER (foisonnement)
fprintf('\n--- Induction de dent (mu=3000) ---\n');
fprintf('flux de dent max      = %.4e Wb\n',PhiT_max);
fprintf('B dent apparente      = %.3f T   (Phi/(wst1*L))\n',Bt_app);
fprintf('B dans le FER (/ki)   = %.3f T   (ki=%.2f)\n',Bt_fe,M.Ki);
fprintf('  -> H(M350)          = %.0f A/m\n',BH.Hof(Bt_fe));
fprintf('  -> mu_r local REEL  = %.0f      <-- a comparer a 3000 !\n',BH.mur(Bt_fe));
fprintf('B sat de la M350      = %.3f T\n',BH.Bsat);

% --- culasse ---
PhiY_max=max(abs(R.PhiY));
By_app=PhiY_max/(M.wsy*M.ls); By_fe=By_app/M.Ki;
fprintf('B culasse apparente   = %.3f T  -> fer %.3f T -> mu_r %.0f\n', ...
    By_app,By_fe,BH.mur(By_fe));
