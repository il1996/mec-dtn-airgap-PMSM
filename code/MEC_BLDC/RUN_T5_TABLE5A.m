%% RUN_T5_TABLE5A - Table 5(a) regeneree HOMOGENE
%
%  T4/T5. Le tableau publie melange deux maillages radiaux : le nombre
%  d'inconnues donne le nombre de couches, et 1800/180 = 10 alors que
%  4320/360 = 6480/540 = 10800/900 = 12. Or nsh=1 produit 10 couches et
%  nsh=2 en produit 12. La ligne Ms=180 -- celle qui porte le -42,3 % cite
%  en evidence sur nu=8 -- est donc seule a nsh=1.
%
%  On regenere ici les QUATRE lignes dans les DEUX configurations, contre
%  les MEMES references, pour que le choix entre regenerer et retirer la
%  ligne Ms=180 se fasse sur des chiffres et non sur une hypothese.
clear; clc; t0=tic;
diary('table5a_out.txt'); diary on;
M=machine_bldc(); p=M.p;
if isfile('kfringe_ident.mat'), S0=load('kfringe_ident.mat'); kfr=S0.kbest; else, kfr=0.325; end
rd=@(f)readmatrix(f,'FileType','text','Delimiter','\t','NumHeaderLines',1);
dL=rd(fullfile(M.FEA.dir,'magnetostique(Armature-Field)','Output Variables Table 1.tab'));
LaF=dL(1,2); MF=(dL(1,3)+dL(1,4))/2; LdF=LaF-MF;

%  References de champ : celles de RUN_SCORE_MESH, qui les lit du projet EF.
Bg1F=1.0746; r8F=0.01966;

thu=linspace(0,2*pi,2001); thu(end)=[];
amp=@(y,th,k)abs((2/numel(y))*sum(y(:).'.*exp(-1i*k*th(:).')));

fprintf('=== Table 5(a) regeneree, solveur LINEAIRE ===\n');
fprintf('  references : Bg1 %.4f T | nu=8 %.5f T | Ld %.3f mH\n',Bg1F,r8F,LdF*1e3);
fprintf('  (Ld_FEA courant 52.342 mH ; le tableau publie utilise 52.24 -- perime)\n\n');
fprintf('  %5s %4s %8s %10s %9s %10s %9s %10s\n', ...
    'Ms','nsh','inconn.','Bg1 (T)','ec.Bg1','Ld (mH)','ec.Ld','ec. nu=8');
for nsh=[1 2]
  for Ms=[180 360 540 900]
    %  --- champ a vide : aimants seuls ---
    ME=mesh_bldc(M,Ms,4,3,[],0,kfr,M.muI,nsh);
    S=solve_bldc_mesh(ME);
    Br=ME.AG.field(S.Usurf,0,thu); Br=Br(:).';
    b1=amp(Br,thu,p); r8=amp(Br,thu,8);
    %  --- inductance : FMM de bobinage, aimants eteints (recette RUN_IND_MESH) ---
    MEa=mesh_bldc(M,Ms,4,3,[1;0;0],0,kfr,M.muI,nsh); MEa.Isrc(:)=0;
    Sa=solve_bldc_mesh(MEa);
    Fb=mesh_bldc(M,Ms,4,3,[0;1;0],0,kfr,M.muI,nsh).E;
    Fc=mesh_bldc(M,Ms,4,3,[0;0;1],0,kfr,M.muI,nsh).E;
    La=sum(MEa.E.*Sa.Phi); Mm=0.5*(sum(Fb.*Sa.Phi)+sum(Fc.*Sa.Phi));
    Ld=La-Mm;
    fprintf('  %5d %4d %8d %10.4f %8.1f %% %10.3f %8.1f %% %8.1f %%\n', ...
        Ms,nsh,ME.N,b1,100*(b1-Bg1F)/Bg1F,Ld*1e3,100*(Ld-LdF)/LdF, ...
        100*(r8-r8F)/r8F);
  end
  fprintf('\n');
end
fprintf('  duree %.0f s\n=== termine ===\n',toc(t0));
diary off;
