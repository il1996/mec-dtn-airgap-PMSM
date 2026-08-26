%% RUN_C4C_TABLE8 - Table 8 (spectre spatial de Br) a PLEINE PRECISION
%
%  C4 (reste). La Table 8 liste TOUS les ordres spatiaux que le modele
%  calcule, et non le sous-ensemble qui soutient l'argument. Ses colonnes
%  d'ecart n'avaient pas de source a pleine precision.
%
%  Ordres : p, |Ns-p|=8, 3p=21, Ns+p=22, |2Ns-p|=23, 5p=35, 2Ns+p=37.
%  Les ordres 23 et 37 sont les SECONDES bandes laterales, gouvernees par
%  la singularite de coin (§8.1) : ils sont rapportes comme indicateurs de
%  resolution de coin, non comme sorties validees.
clear; clc; t0=tic;
diary('C4c_table8_out.txt'); diary on;
M=machine_bldc(); Ns=M.Ns; p=M.p;
if isfile('kfringe_ident.mat'), S0=load('kfringe_ident.mat'); kfr=S0.kbest; else, kfr=0.325; end
rd=@(f)readmatrix(f,'FileType','text','Delimiter','\t','NumHeaderLines',1);
fea=M.FEA.dir;
d4=rd(fullfile(fea,'magnetostique(Magnetic_loading)','Calculator Expressions Plot 4.tab'));
thu=linspace(0,2*pi,3601); thu(end)=[];
BrF=interp1([d4(:,2)*pi/180;2*pi],[d4(:,3);d4(1,3)],thu,'linear','extrap');
amp=@(B,th,k) 2*abs(mean(B(:).'.*exp(-1i*k*th(:).')));

fprintf('=== C4c : Table 8, spectre spatial de Br, pleine precision ===\n');
fprintf('  PMSM 15/14 | Ns = %d, p = %d | kfringe = %.4f\n',Ns,p,kfr);
fprintf('  reference EF relue du projet, aucune transcription\n\n');

%  ---- reseau localise (colonne "Network" de la Table 8) ----
R=cogging_mec(M,1260,0,3,M.muI,kfr,1e-6);
Br=R.AG.field(R.U(1:1260,1),0,thu); Br=Br(:).';

ord=[p 8 3*p Ns+p abs(2*Ns-p) 5*p 2*Ns+p];
lab={'p = 7 (fondamental)','|Ns - p| = 8','3p = 21','Ns + p = 22', ...
     '|2Ns - p| = 23','5p = 35','2Ns + p = 37'};
note={'','1re bande','harmonique de travail','1re bande', ...
      '2e bande (coin)','harmonique de travail','2e bande (coin)'};

fprintf('  %-24s %14s %14s %12s  %s\n','ordre','reseau','EF','ecart','nature');
V=nan(numel(ord),3);
for k=1:numel(ord)
    a=amp(Br,thu,ord(k)); f=amp(BrF,thu,ord(k));
    V(k,:)=[a f 100*(a-f)/f];
    fprintf('  %-24s %14.7f %14.7f %11.4f %%  %s\n',lab{k},a,f,V(k,3),note{k});
end

fprintf('\n  ---- lecture ----\n');
fprintf('  harmoniques de TRAVAIL (3p, 5p) : |ecart| max %.3f %%\n', ...
    max(abs(V([3 6],3))));
fprintf('  PREMIERES bandes (8, 22)        : %.3f %% et %.3f %%\n',V(2,3),V(4,3));
fprintf('  SECONDES bandes (23, 37)        : %.3f %% et %.3f %%\n',V(5,3),V(7,3));
fprintf(['  Les secondes bandes echantillonnent la region de coin plus\n' ...
         '  finement que les premieres -- leur longueur d''onde angulaire\n' ...
         '  vaut environ la moitie de l''ouverture d''encoche -- donc la meme\n' ...
         '  singularite non resolue y produit une erreur relative plus\n' ...
         '  grande. A rapporter comme INDICATEUR de resolution de coin.\n']);
save('C4c_table8.mat','ord','lab','V');
fprintf('\n  duree %.0f s\n=== C4c termine ===\n',toc(t0));
diary off;
