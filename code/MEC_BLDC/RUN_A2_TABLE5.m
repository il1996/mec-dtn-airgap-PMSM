%% RUN_A2_TABLE5 - Tables 5(a) et 5(b) reconstruites
%
%  A2. Le tableau publie melange trois provenances : la colonne "unknowns"
%  vient de n_sh = 2 (4320/6480/10800), tandis que B_g1 et nu=8 viennent de
%  n_sh = 1, et L_d de n_sh = 2 contre une reference PERIMEE (52,24 mH au
%  lieu de 52,342). On reconstruit donc :
%
%    Table 5(a) : ENTIEREMENT a n_sh = 1, inconnues 1800/3600/5400/9000,
%                 contre L_d,FEA relu du projet.
%    Table 5(b) : la serie n_sh = 2, solveurs LINEAIRE et NEWTON cote a cote.
%                 Elle montre que le raffinement radial du bec DEPLACE la
%                 bande de denture au lieu de la faire converger.
%
%  Colonne "peak iron B" ajoutee aux deux panneaux : c'est elle qui signe la
%  singularite -- un champ qui croit sans se stabiliser sous raffinement.
clear; clc; t0=tic;
diary('A2_table5_out.txt'); diary on;
M=machine_bldc(); p=M.p;
if isfile('kfringe_ident.mat'), S0=load('kfringe_ident.mat'); kfr=S0.kbest; else, kfr=0.325; end
rd=@(f)readmatrix(f,'FileType','text','Delimiter','\t','NumHeaderLines',1);
fea=M.FEA.dir;
dL=rd(fullfile(fea,'magnetostique(Armature-Field)','Output Variables Table 1.tab'));
LaF=dL(1,2); MF=(dL(1,3)+dL(1,4))/2; LdF=LaF-MF;
d4=rd(fullfile(fea,'magnetostique(Magnetic_loading)','Calculator Expressions Plot 4.tab'));
thu=linspace(0,2*pi,2001); thu(end)=[];
BrF=interp1([d4(:,2)*pi/180;2*pi],[d4(:,3);d4(1,3)],thu,'linear','extrap');
amp=@(y,th,k)abs((2/numel(y))*sum(y(:).'.*exp(-1i*k*th(:).')));
Bg1F=amp(BrF,thu,p); r8F=amp(BrF,thu,8);

fprintf('=== A2 : Tables 5(a) et 5(b) reconstruites ===\n');
fprintf('  references EF relues du projet :\n');
fprintf('    Bg1 = %.5f T | nu=8 = %.6f T | Ld = %.4f mH\n',Bg1F,r8F,LdF*1e3);
fprintf('  (le tableau publie utilisait Ld_FEA = 52.24 mH -- PERIME)\n');

%  NOTE. Les fonctions locales d'un SCRIPT ne partagent pas l'espace de
%  travail (contrairement aux fonctions imbriquees d'une fonction) : tout
%  ce dont pt a besoin lui est passe explicitement.
pt=@(Ms,nsh,sol) point5(M,Ms,nsh,sol,kfr,thu,p);

%% ---- Table 5(a) : n_sh = 1, solveur lineaire ------------------------
fprintf('\n  ---- Table 5(a) : n_sh = 1, solveur LINEAIRE ----\n');
fprintf('  %6s %10s %10s %10s %11s %12s\n', ...
    'Ms','inconnues','Bg1','Ld','nu=8','peak iron B');
fprintf('  %6s %10s %10s %10s %11s %12s\n','','','ecart','ecart','ecart','(T)');
for Ms=[180 360 540 900]
    [b1,r8,Ld,pk,N]=pt(Ms,1,M.muI);
    fprintf('  %6d %10d %9.1f %% %9.1f %% %10.1f %% %12.4f\n', ...
        Ms,N,100*(b1-Bg1F)/Bg1F,100*(Ld-LdF)/LdF,100*(r8-r8F)/r8F,pk);
end

%% ---- Table 5(b) : n_sh = 2, lineaire ET Newton ----------------------
fprintf('\n  ---- Table 5(b) : n_sh = 2, LINEAIRE et NEWTON ----\n');
fprintf('  %6s %10s %12s %12s %12s %12s\n', ...
    'Ms','inconnues','nu=8 lin.','nu=8 Newton','peak B lin.','peak B NL');
for Ms=[360 540 900]
    [~,r8l,~,pkl,N]=pt(Ms,2,M.muI);
    [~,r8n,~,pkn,~]=pt(Ms,2,'nl');
    fprintf('  %6d %10d %11.1f %% %11.1f %% %12.4f %12.4f\n', ...
        Ms,N,100*(r8l-r8F)/r8F,100*(r8n-r8F)/r8F,pkl,pkn);
end

fprintf(['\n  LECTURE. En 5(a) le raffinement ANGULAIRE fait remonter nu=8\n' ...
         '  (-42 -> -10 %%). En 5(b) le raffinement RADIAL du bec le fait\n' ...
         '  redescendre. Les deux directions bougent la bande en sens\n' ...
         '  OPPOSES, et le champ de fer crete croit sans se stabiliser :\n' ...
         '  c''est la signature d''une singularite, non d''une convergence.\n']);
fprintf('\n  duree %.0f s\n=== A2 termine ===\n',toc(t0));
diary off;

% ======================================================================
    function [b1,r8,Ld,pkB,N] = point5(M,Ms,nsh,sol,kfr,thu,p)
        amp=@(y,th,k)abs((2/numel(y))*sum(y(:).'.*exp(-1i*k*th(:).')));
        MEf=mesh_bldc(M,Ms,4,3,[],0,kfr,sol,nsh);
        Sf=solve_bldc_mesh(MEf,1e-9,30);
        Br=MEf.AG.field(Sf.Usurf,0,thu); Br=Br(:).';
        b1=amp(Br,thu,p); r8=amp(Br,thu,8);
        pkB=max(abs(Sf.B(MEf.iron)));
        MEa=mesh_bldc(M,Ms,4,3,[1;0;0],0,kfr,sol,nsh); MEa.Isrc(:)=0;
        Sa=solve_bldc_mesh(MEa,1e-9,30);
        Fb=mesh_bldc(M,Ms,4,3,[0;1;0],0,kfr,sol,nsh).E;
        Fc=mesh_bldc(M,Ms,4,3,[0;0;1],0,kfr,sol,nsh).E;
        La=sum(MEa.E.*Sa.Phi); Mu=0.5*(sum(Fb.*Sa.Phi)+sum(Fc.*Sa.Phi));
        Ld=(La-Mu); N=MEf.N;
    end
