%% RUN_DIAG_CARTER  -  Le modele a-t-il la BONNE amplitude de denture ?
%  Test decisif AVANT d'ajouter de la physique : on compare l'effet MOYEN de
%  denture du sous-domaine exact au coefficient de CARTER analytique, qui est
%  une reference independante et classique.
%     - si le modele retrouve Carter -> son traitement de l'encoche est juste,
%       et l'ecart de -28 % vs FEA vient d'ailleurs (fer equipotentiel => la
%       face segmentee ne changerait rien) ;
%     - sinon -> il y a une erreur de modele a corriger, et c'est la priorite.
clear; clc;
M=machine_bldc(); mu0=4*pi*1e-7;
K=8; mu=1e7;                       % fer infiniment permeable : cas pur Laplace
geff=M.lag + M.hm/M.mu_r;          % entrefer magnetique equivalent
taus=2*pi*M.Rsi/M.Ns;              % pas d'encoche (m)

fprintf('=== Effet MOYEN de denture : modele vs coefficient de Carter ===\n');
fprintf('entrefer magnetique equivalent g'''' = %.3f mm ; pas d''encoche = %.3f mm\n\n', ...
    geff*1e3,taus*1e3);

% reference "lisse" du MEME code : ouverture quasi nulle
Ms=M; Ms.ws0=0.02e-3; nm=ceil(4*K*pi/(Ms.ws0/Ms.Rsi));
Rs0=subdomain_mec(Ms,K,min(nm,4000),mu,0);
fprintf('%8s | %10s %10s | %9s %9s | %8s\n', ...
    'ws0[mm]','Bg1 modele','Bg1/Bg1_0','kC modele','kC Carter','ecart');
for w=[0.5 1 1.5 2 2.5 3 4]
    M2=M; M2.ws0=w*1e-3; nm=ceil(4*K*pi/(M2.ws0/M2.Rsi));
    R=subdomain_mec(M2,K,nm,mu,0);
    kC_mod=Rs0.Bg1/R.Bg1;                       % Carter "vu" par le modele
    u=(w*1e-3)/geff;
    sg=(2/pi)*( atan(u/2) - (1/u)*log(1+(u/2)^2) );
    kC_car=taus/(taus - sg*(w*1e-3));
    fprintf('%8.2f | %10.4f %10.5f | %9.5f %9.5f | %+7.1f%%\n', ...
        w,R.Bg1,R.Bg1/Rs0.Bg1,kC_mod,kC_car,100*(kC_mod-1)/(kC_car-1)-100);
end
fprintf('\n(ecart = ecart sur la PART de denture (kC-1), la grandeur physique)\n');
fprintf('Bg1 lisse (ws0->0) = %.4f T\n',Rs0.Bg1);
