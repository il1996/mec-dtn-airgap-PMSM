%% RUN_OPT_H8 - l'harmonique 8 commande 4 entrees du scorecard
%  Les plus gros ecarts restants ont TOUS la meme origine : le rang
%  nu = |Ns-p| = 8 du champ d'entrefer, sous-estime de 40 %. C'est lui qui
%  porte (a) la dispersion de FMM d'un aimant a l'autre (la denture vue par
%  le rotor), (b) les pertes aimant a vide, (c) les pertes aimant en charge
%  -- ces deux dernieres en B^2. Une correction du rang 8 vaudrait donc
%  quatre lignes du scorecard d'un coup.
%
%  Le seul parametre qui gouverne la MODULATION d'encoche est kfringe : la
%  conductance de frange de la bouche d'encoche. kfringe eleve = ouverture
%  "transparente" = modulation faible ; kfringe nul = ouverture isolante =
%  modulation maximale. Il avait ete identifie a 0.75 sur les bandes de
%  denture du champ FEA, mais par un critere GLOBAL (ecart RMS) qui peut
%  avoir sacrifie la bande laterale au profit du fondamental.
clear; clc; t0=tic;
M=machine_bldc(); Ns=M.Ns; p=M.p; fea=M.FEA.dir;
rd=@(f)readmatrix(f,'FileType','text','Delimiter','\t','NumHeaderLines',1);
d4=rd(fullfile(fea,'magnetostique(Magnetic_loading)','Calculator Expressions Plot 4.tab'));
thu=linspace(0,2*pi,3601); thu(end)=[];
BrFu=interp1([d4(:,2)*pi/180;2*pi],[d4(:,3);d4(1,3)],thu,'linear','extrap');
amp=@(B,th,k) 2*abs(mean(B(:).'.*exp(-1i*k*th(:).')));
aF=@(k)amp(BrFu,thu,k);

fprintf('=== Sensibilite du champ d''entrefer a kfringe ===\n');
fprintf('  cible FEA : Bg1 = %.5f T | rang 8 = %.5f T | rang 22 = %.5f T\n', ...
    aF(p),aF(8),aF(22));
fprintf('  %8s %10s %9s %10s %9s %10s %9s %10s\n', ...
    'kfringe','Bg1 (T)','ecart','rang 8','ecart','rang 22','ecart','RMS (T)');
kv=[0 0.15 0.30 0.45 0.60 0.75 0.90];
B1=zeros(size(kv)); B8=B1; B22=B1; RM=B1;
for q=1:numel(kv)
    R=cogging_mec(M,1260,0,3,M.muI,kv(q),1e-6);
    th=R.thq(1:end-1); Br=R.Br(1:end-1);
    B1(q)=amp(Br,th,p); B8(q)=amp(Br,th,8); B22(q)=amp(Br,th,22);
    Bi=interp1([th 2*pi],[Br Br(1)],thu,'linear','extrap');
    %  alignement optimal avant l'ecart RMS (le repere d'origine differe)
    e=inf;
    for s=0:5:355
        e=min(e,sqrt(mean((interp1([thu 2*pi],[Bi Bi(1)], ...
            mod(thu+s*pi/180,2*pi))-BrFu).^2)));
    end
    RM(q)=e;
    fprintf('  %8.2f %10.5f %8.1f %% %10.5f %8.1f %% %10.5f %8.1f %% %10.5f\n', ...
        kv(q),B1(q),100*(B1(q)-aF(p))/aF(p),B8(q),100*(B8(q)-aF(8))/aF(8), ...
        B22(q),100*(B22(q)-aF(22))/aF(22),RM(q));
end

%% ---- que gagnerait-on sur le scorecard ? ------------------------------
[~,qb]=min(abs(B8-aF(8)));
fprintf('\n  --- meilleur compromis sur le rang 8 : kfringe = %.2f ---\n',kv(qb));
fprintf('  Bg1 %+.1f %% | rang 8 %+.1f %% | rang 22 %+.1f %% | RMS %.5f T\n', ...
    100*(B1(qb)-aF(p))/aF(p),100*(B8(qb)-aF(8))/aF(8), ...
    100*(B22(qb)-aF(22))/aF(22),RM(qb));
fprintf('  gain attendu sur les pertes aimant (en B^2) : x%.2f\n',(B8(qb)/B8(6))^2);
fprintf('  -> a vide  %.4f W -> %.4f W (cible %.4f)\n', ...
    0.1679,0.1679*(B8(qb)/B8(6))^2,0.3335);
fprintf('  -> en charge %.3f W -> %.3f W (cible %.3f)\n', ...
    1.529,1.529*(B8(qb)/B8(6))^2,3.183);
fprintf('\n  (%.0f s)\n',toc(t0));
