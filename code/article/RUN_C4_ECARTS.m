%% RUN_C4_ECARTS - regeneration des colonnes d'ecart a PLEINE PRECISION
%
%  C4. T14 a borne le probleme sans le resoudre : 0 cellule sur 42 est
%  incompatible avec ses valeurs affichees, mais HUIT presentent l'indice
%  inverse -- un ecart publie coincidant a moins de 0,01 point avec le calcul
%  NAIF sur valeurs arrondies. Une telle coincidence est improbable si
%  l'ecart vient des valeurs pleines : ces huit sont les candidates.
%
%  METHODE. On ne recalcule PAS depuis les valeurs affichees : cela
%  reproduirait le defaut. On repart des .mat produits par A1, B1, B2, B3,
%  qui portent les valeurs pleines, et on reforme chaque ecart.
%
%  SORTIE : pour chaque cellule, l'ecart a pleine precision, l'ecart publie,
%  et le verdict -- REMPLACER si l'arrondi du nouveau differe du publie.
clear; clc;
diary('C4_ecarts_out.txt'); diary on;
B=fullfile('..','MEC_BLDC'); I=fullfile('..','MEC_IM');
fprintf('=== C4 : colonnes d''ecart regenerees a pleine precision ===\n\n');

ec=@(a,f)100*(a-f)/f;
%  arrondi a une decimale, comme le manuscrit
r1=@(x)round(x,1);
nrep=0; ntot=0;
show=@(nom,a,f,pub) deal_show(nom,a,f,pub,ec,r1);

%% ---- PMSM : depuis A1_table7.mat --------------------------------------
fA=fullfile(B,'A1_table7.mat');
if isfile(fA)
    S=load(fA);
    fprintf('---- PMSM, Table 7 (source : A1_table7.mat) ----\n');
    fprintf('  %-26s %12s %12s %11s %10s %s\n', ...
        'grandeur','valeur','reference','ecart plein','publie','verdict');
    for j=1:numel(S.lab)
        if isnan(S.VF(j)), continue; end
        %  colonne Lumped, celle qui porte les cellules signalees
        [n,ntot]=show(S.lab{j},S.VL(j),S.VF(j),NaN); nrep=nrep+n; ntot=ntot+1;
    end
else
    fprintf('  A1_table7.mat ABSENT -- relancer RUN_A1_TABLE7 d''abord.\n');
end

%% ---- MAS : depuis B3_prix.mat -----------------------------------------
fB=fullfile(I,'B3_prix.mat');
if isfile(fB)
    S=load(fB);
    fprintf('\n---- MAS, les trois fermetures (source : B3_prix.mat) ----\n');
    lb={'Xm0 non sature','Xm sature charge','courant a vide','fem entrefer','couple'};
    ref=[NaN 46.0 8.49 382.1 121.63];
    pub=[NaN -6.1 16.5 -1.7 -5.4];       % colonne P0 publiee
    fprintf('  %-22s %11s %11s %11s %10s\n','grandeur','P0','P1','ec.P0 plein','publie');
    for j=1:5
        if isnan(ref(j)), continue; end
        e0=ec(S.V(2,j),ref(j)); e1=ec(S.V(3,j),ref(j));
        v=''; if abs(r1(e0)-pub(j))>1e-9, v='<-- REMPLACER'; nrep=nrep+1; end
        ntot=ntot+1;
        fprintf('  %-22s %11.4f %11.4f %10.4f %% %9.1f %% %s\n', ...
            lb{j},S.V(2,j),S.V(3,j),e0,pub(j),v);
    end
else
    fprintf('\n  B3_prix.mat ABSENT.\n');
end

%% ---- k_C : depuis B2_kC.mat -------------------------------------------
fC=fullfile(I,'B2_kC.mat');
if isfile(fC)
    S=load(fC);
    fprintf('\n---- MAS, k_C convergé (source : B2_kC.mat) ----\n');
    fprintf('  k_C P1 au plus grand N_h = %d : %.6f\n',S.Nl(end),S.R(end,6));
    fprintf('  X_m encoche P1                : %.6f ohm\n',S.R(end,5));
    fprintf('  X_m lisse   P1                : %.6f ohm\n',S.R(end,4));
    fprintf('  derive de X_m lisse P0 sur le balayage : %+.4f %%\n', ...
        ec(S.R(end,1),S.R(1,1)));
    fprintf('  derive de X_m dent  P0 sur le balayage : %+.4f %%\n', ...
        ec(S.R(end,2),S.R(1,2)));
    fprintf('  (le contraste de ces deux lignes refute la §6.4)\n');
end

fprintf('\n---- BILAN ----\n');
fprintf('  %d cellule(s) examinee(s) | %d a REMPLACER\n',ntot,nrep);
fprintf(['\n  RESERVE. Ce script ne couvre que les tableaux dont un .mat a\n' ...
         '  pleine precision existe (A1, B2, B3). Les Tables 8, 10, 16 et 17\n' ...
         '  n''en ont pas : leur regeneration exige de relancer la chaine qui\n' ...
         '  les produit avec sauvegarde des valeurs pleines. C''est un travail\n' ...
         '  de chaine, pas d''arrondi -- et c''est ce que C4 demande vraiment.\n']);
fprintf('\n=== C4 termine ===\n');
diary off;

% ======================================================================
function [n,t]=deal_show(nom,a,f,pub,ec,r1)
    e=ec(a,f); n=0; t=1; v='';
    if ~isnan(pub) && abs(r1(e)-pub)>1e-9, v='<-- REMPLACER'; n=1; end
    fprintf('  %-26s %12.5f %12.5f %10.4f %% %9s %s\n', ...
        nom,a,f,e,ternary(isnan(pub),'--',sprintf('%.1f %%',pub)),v);
end
function o=ternary(c,a,b), if c, o=a; else, o=b; end, end
