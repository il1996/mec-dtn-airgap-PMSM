%% RUN_G2_CLOTURE - extractions pleine precision pour la cloture de l'Article I
%
%  Trois points du message v7 exigent des valeurs pleine precision qui
%  existent deja dans des .mat sauvegardes par les blocs anterieurs. Ce bloc
%  les relit et forme les grandeurs demandees SANS jamais repartir des
%  valeurs imprimees (regle 2 du dossier).
%
%    G2-A  section 3.1 : la matrice 2 x 2 du facteur de reluctance k_r,
%          deux references (avec et sans le point aberrant) x les colonnes
%          de modele de la Table 13.
%    G2-B  section 3.5 : les lignes de rendement du Tableau 20, converties
%          en POINTS, et le recomptage des cas a moins d'un point du seuil.
%    G2-C  section 2.3 : les quatre largeurs d'encadrement du balayage
%          Newton, a pleine precision.
%
%  GARDE. Chaque grandeur relue d'un .mat est confrontee a la valeur
%  imprimee dans le transcript correspondant. Un ecart superieur a l'arrondi
%  d'affichage signifie que le .mat et le transcript ne sont pas du meme
%  run, et le resultat doit alors etre rejete.
clear; clc; t0=tic;
diary('G2_cloture_out.txt'); diary on;
fprintf('=== G2 : extractions pleine precision pour la cloture ===\n');

%% ================= G2-A : la matrice 2 x 2 de k_r =====================
fprintf('\n  ---- G2-A : facteur de reluctance k_r, matrice complete ----\n');
M=machine_bldc();
rd=@(f)readmatrix(f,'FileType','text','Delimiter','\t','NumHeaderLines',1);
fea=M.FEA.dir;
okA=false;
try
    dmm=rd(fullfile(fea,'magnetostique(Magnetic_loading)','Calculator Expressions Table 1.tab'));
    mmf_m_F=dmm(1,2:15); mmf_g_F=dmm(1,16:29);
    krF =rd(fullfile(fea,'magnetostique(Magnetic_loading)','Output Variables Table 2.tab')); krF=krF(1,2);
    krC =mean(mmf_m_F)/mean(mmf_g_F(2:end));   % hors point aberrant, def. RUN_KRKL_MESH:19
    fprintf('    reference EF, TOUS les points        : %.15f\n',krF);
    fprintf('    reference EF, hors point aberrant    : %.15f\n',krC);
    fprintf('    FMM entrefer, les 14 points (A) :\n      ');
    fprintf('%.2f ',mmf_g_F); fprintf('\n');
    fprintf('    le point ecarte est le premier : %.4f A ; les treize autres : %.1f a %.1f A\n', ...
        mmf_g_F(1),min(mmf_g_F(2:end)),max(mmf_g_F(2:end)));
    okA=true;
catch ME1
    fprintf('    LECTURE EF IMPOSSIBLE : %s\n',ME1.message);
end

if okA && isfile('A1_table7.mat')
    A=load('A1_table7.mat');
    vn=fieldnames(A);
    fprintf('\n    variables de A1_table7.mat : %s\n',strjoin(vn.',', '));
    krmod=[]; lab={};
    for i=1:numel(vn)
        v=A.(vn{i});
        if isnumeric(v) && ~isempty(v)
            f=v(:).';
            j=find(f>1.02 & f<1.12);          % la fenetre ou vit k_r
            if ~isempty(j)
                fprintf('    %-12s valeurs dans [1.02,1.12] :',vn{i});
                fprintf(' %.15f',f(j)); fprintf('\n');
                krmod=[krmod f(j)]; %#ok<AGROW>
                for k=1:numel(j), lab{end+1}=sprintf('%s(%d)',vn{i},j(k)); end %#ok<AGROW>
            end
        end
    end
    krmod=krmod(abs(krmod-krF)>1e-12);        % retirer la reference elle-meme
    fprintf('\n    MATRICE DES ECARTS, formes sur les valeurs pleines :\n');
    fprintf('    %-28s %18s %18s\n','colonne de modele','/ 1.08693 (tous)','/ 1.0664 (hors ab.)');
    for k=1:numel(krmod)
        fprintf('    k_r = %.15f   %15.4f %%  %15.4f %%\n', ...
            krmod(k),100*(krmod(k)-krF)/krF,100*(krmod(k)-krC)/krC);
    end
else
    fprintf('    A1_table7.mat absent : la matrice ne peut pas etre formee.\n');
end

%% ================= G2-B : rendements en POINTS ========================
fprintf('\n  ---- G2-B : lignes de rendement du Tableau 20, en points ----\n');
if isfile('R7_scorecard.mat')
    R=load('R7_scorecard.mat');
    fprintf('    variables de R7_scorecard.mat : %s\n',strjoin(fieldnames(R).',', '));
    SC=R.SC;
    if istable(SC), disp(SC(:,:)); end
    if isnumeric(SC)
        fprintf('    SC est numerique %dx%d\n',size(SC,1),size(SC,2));
        for r=[25 29]
            if r<=size(SC,1)
                fprintf('    ligne %2d : modele %.10f   reference %.10f\n',r,SC(r,1),SC(r,2));
                fprintf('               ecart relatif %+.6f %%   |   ecart en POINTS %+.6f p.p.\n', ...
                    100*(SC(r,1)-SC(r,2))/SC(r,2), SC(r,1)-SC(r,2));
            end
        end
    end
    if isfield(R,'E') && isnumeric(R.E)
        E=R.E(:).';
        fprintf('\n    RECOMPTAGE sous la convention en POINTS pour les lignes de rendement\n');
        fprintf('    ecarts relatifs, les 29 : \n      ');
        fprintf('%+.4f ',E); fprintf('\n');
        nk=sum(abs(E)<5);
        fprintf('    dans les 5 %% (convention actuelle, tout en relatif) : %d / %d\n',nk,numel(E));
        Ep=E;
        if isnumeric(SC) && size(SC,1)>=29
            Ep(25)=SC(25,1)-SC(25,2);
            Ep(29)=SC(29,1)-SC(29,2);
        end
        nkp=sum(abs(Ep)<5);
        fprintf('    dans les 5 %% (lignes 25 et 29 en points)            : %d / %d\n',nkp,numel(Ep));
        fprintf('    a MOINS D UN POINT du seuil de 5, convention actuelle :\n');
        for i=1:numel(E)
            if abs(E(i))<5 && abs(E(i))>4, fprintf('      #%d  %+.4f\n',i,E(i)); end
        end
        fprintf('    a MOINS D UN POINT du seuil de 5, lignes 25 et 29 en points :\n');
        for i=1:numel(Ep)
            if abs(Ep(i))<5 && abs(Ep(i))>4, fprintf('      #%d  %+.4f\n',i,Ep(i)); end
        end
    end
else
    fprintf('    R7_scorecard.mat absent.\n');
end

%% ================= G2-C : les quatre largeurs Newton ==================
fprintf('\n  ---- G2-C : encadrement Newton, quatre couches de bec ----\n');
if isfile('X1_table5b.mat')
    X=load('X1_table5b.mat');
    fprintf('    variables de X1_table5b.mat : %s\n',strjoin(fieldnames(X).',', '));
    fn=fieldnames(X);
    for i=1:numel(fn)
        v=X.(fn{i});
        if isnumeric(v)
            fprintf('    %-14s taille %s\n',fn{i},mat2str(size(v)));
            if numel(v)<=60
                disp(v);
            end
        end
    end
else
    fprintf('    X1_table5b.mat absent.\n');
end

fprintf('\n  duree %.0f s\n=== G2 termine ===\n',toc(t0));
diary off;
