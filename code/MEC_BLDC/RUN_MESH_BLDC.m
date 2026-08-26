%% RUN_MESH_BLDC - premiere mise en route du maillage raffine BLDC
%  Controle de validite AVANT toute carte : le reseau raffine doit
%  retrouver le champ d'entrefer du reseau a une dent (deja valide a +0.3 %
%  sur Bg1), sinon la carte serait jolie et fausse.
clear; clc;
M=machine_bldc(); Ns=M.Ns; p=M.p;
if isfile('kfringe_ident.mat'), S0=load('kfringe_ident.mat'); kfr=S0.kbest; else, kfr=0.325; end
rd=@(f)readmatrix(f,'FileType','text','Delimiter','\t','NumHeaderLines',1);
d4=rd(fullfile(M.FEA.dir,'magnetostique(Magnetic_loading)','Calculator Expressions Plot 4.tab'));
thu=linspace(0,2*pi,3601); thu(end)=[];
BrF=interp1([d4(:,2)*pi/180;2*pi],[d4(:,3);d4(1,3)],thu,'linear','extrap');
amp=@(B,th,k) 2*abs(mean(B(:).'.*exp(-1i*k*th(:).')));

fprintf('=== Maillage raffine BLDC : mise en route ===\n');
for Ms=[180 360 540]
    t0=tic;
    ME=mesh_bldc(M,Ms,4,3,[],0,kfr,M.muI);
    S=solve_bldc_mesh(ME);
    [Br,~]=ME.AG.field(S.Usurf,0,thu);
    nfe=sum(ME.iron); nair=sum(~ME.iron);
    fprintf('  Ms=%3d : %d noeuds, %d branches (%d fer / %d air), %.1f s\n', ...
        Ms,ME.N,numel(ME.a),nfe,nair,toc(t0));
    fprintf('           Bg1 = %.4f T (FEA %.4f, %+.1f %%) | rang 8 %+.1f %% | B moyen %+.1f %%\n', ...
        amp(Br,thu,p),amp(BrF,thu,p),100*(amp(Br,thu,p)-amp(BrF,thu,p))/amp(BrF,thu,p), ...
        100*(amp(Br,thu,8)-amp(BrF,thu,8))/amp(BrF,thu,8), ...
        100*(mean(abs(Br))-mean(abs(BrF)))/mean(abs(BrF)));
end

%% ---- non lineaire ----
fprintf('\n  --- resolution NON LINEAIRE (courbe B(H) M350-50A) ---\n');
ME=mesh_bldc(M,360,4,3,[],0,kfr,'nl');
S=solve_bldc_mesh(ME);
[Br,~]=ME.AG.field(S.Usurf,0,thu);
Bfe=abs(S.B(ME.iron));
fprintf('  %d iterations, residu %.2e\n',S.iter,S.err);
fprintf('  Bg1 = %.4f T (%+.1f %%) | B fer max %.3f T | mu_r de %.0f a %.0f\n', ...
    amp(Br,thu,p),100*(amp(Br,thu,p)-amp(BrF,thu,p))/amp(BrF,thu,p), ...
    max(Bfe),min(S.mur(ME.iron)),max(S.mur(ME.iron)));
save('mesh_bldc_nl.mat','ME','S','-v7.3');
fprintf('\n-> mesh_bldc_nl.mat\n');
