%% RUN_DIAG_H8_MESH - pourquoi le rang 8 s'inverse entre lineaire et Newton
%
%  CONSTAT. Au meme maillage (Ms=540), la bande de denture rang 8 = |Ns-p|
%  passe de -13.3 % en LINEAIRE a +22.3 % en NEWTON.
%
%  HYPOTHESE. Le sens du basculement est PHYSIQUE : la saturation des cornes
%  de bec fait chuter leur permeabilite, ce qui CREUSE la permeance
%  d'entrefer au droit de l'ouverture et RENFORCE donc la modulation
%  d'encoche. Un modele lineaire ne peut que la sous-estimer -- d'ou -13 %.
%  Si c'est bien cela, l'ampleur (+22 %) doit dependre de la finesse avec
%  laquelle le BEC est maille : la saturation s'y concentre sur quelques
%  dixiemes de millimetre, et une maille trop grossiere fait saturer tout
%  le bloc au lieu de la seule corne -> modulation surestimee.
%
%  TEST. Convergence du rang 8 en (Ms, nsh) sur le cas NON LINEAIRE, avec le
%  cas lineaire en temoin. Si le rang 8 non lineaire descend vers la FEA
%  quand on raffine le bec, l'hypothese est confirmee et l'ecart est
%  numerique. S'il stagne a +22 %, elle est fausse.
clear; clc;
M=machine_bldc(); p=M.p;
if isfile('kfringe_ident.mat'), S0=load('kfringe_ident.mat'); kfr=S0.kbest; else, kfr=0.325; end
rd=@(f)readmatrix(f,'FileType','text','Delimiter','\t','NumHeaderLines',1);
d4=rd(fullfile(M.FEA.dir,'magnetostique(Magnetic_loading)','Calculator Expressions Plot 4.tab'));
thu=linspace(0,2*pi,3601); thu(end)=[];
BrF=interp1([d4(:,2)*pi/180;2*pi],[d4(:,3);d4(1,3)],thu,'linear','extrap');
amp=@(B,th,k) 2*abs(mean(B(:).'.*exp(-1i*k*th(:).')));
t8=amp(BrF,thu,8); t1=amp(BrF,thu,p);
fprintf('=== Rang 8 : lineaire vs Newton, convergence en maillage de bec ===\n');
fprintf('  cible FEA : rang 8 = %.5f T | Bg1 = %.4f T\n\n',t8,t1);
fprintf('  %5s %5s %5s | %9s %8s | %9s %8s | %9s\n', ...
    'Ms','nsh','Ls','rang8 lin','ecart','rang8 NL','ecart','Bfe max');
for Ms=[360 540 900]
    for nsh=[1 2 4]
        ME=mesh_bldc(M,Ms,4,3,[],0,kfr,M.muI,nsh);
        S=solve_bldc_mesh(ME);
        [BrL,~]=ME.AG.field(S.Usurf,0,thu);
        MEn=mesh_bldc(M,Ms,4,3,[],0,kfr,'nl',nsh);
        Sn=solve_bldc_mesh(MEn,1e-9,30);
        [BrN,~]=MEn.AG.field(Sn.Usurf,0,thu);
        fprintf('  %5d %5d %5d | %9.5f %7.1f %% | %9.5f %7.1f %% | %9.3f\n', ...
            Ms,nsh,ME.Ls,amp(BrL,thu,8),100*(amp(BrL,thu,8)-t8)/t8, ...
            amp(BrN,thu,8),100*(amp(BrN,thu,8)-t8)/t8,max(abs(Sn.B(MEn.iron))));
    end
end

%% ---- ou se sature le fer ? -------------------------------------------
ME=mesh_bldc(M,900,4,3,[],0,kfr,'nl',4); S=solve_bldc_mesh(ME,1e-9,30);
B=abs(S.B); Fe=ME.iron;
fprintf('\n  --- repartition de la saturation (Ms=900, nsh=4) ---\n');
fprintf('  %-28s %8s %8s %8s\n','region','B max','B moy','mu_r min');
lab={'isthme+biseau (bec)','corps de dent','culasse'};
sel={ME.lay<=2*4, ME.lay>2*4 & ME.lay<=2*4+4, ME.lay>2*4+4};
for k=1:3
    m=Fe & sel{k}(:);
    fprintf('  %-28s %8.3f %8.3f %8.0f\n',lab{k},max(B(m)),mean(B(m)),min(S.mur(m)));
end
fprintf('\n  Bg1 = %.4f T (%+.1f %%) | rang 8 = %.5f T (%+.1f %%)\n', ...
    amp(BrN,thu,p),100*(amp(BrN,thu,p)-t1)/t1,amp(BrN,thu,8),100*(amp(BrN,thu,8)-t8)/t8);
