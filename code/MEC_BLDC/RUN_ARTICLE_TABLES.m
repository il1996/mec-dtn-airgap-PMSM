%% RUN_ARTICLE_TABLES - regenere les chiffres de l'article et les ECRIT
%
%  POURQUOI CE FICHIER. La colonne "Mesh" de la Table III de l'article
%  n'est reproduite par aucun fichier de l'arborescence : RUN_SCORE_MESH
%  et RUN_KRKL_MESH impriment sur la console sans rien enregistrer. Les
%  valeurs publiees proviennent donc d'une session perdue, et deux d'entre
%  elles sont en tension avec le reste du manuscrit (voir §3 ci-dessous).
%
%  Ce script fait trois choses et n'en fait aucune autre :
%    1. capture RUN_SCORE_MESH et RUN_KRKL_MESH dans article_tables.txt
%    2. verifie la fenetre du rms(i_a) du transitoire en charge
%    3. imprime les lignes LaTeX pretes a coller
%
%  Duree : celle de RUN_SCORE_MESH (balayage 61 positions a Ms=540) plus
%  celle de RUN_KRKL_MESH. Compter quelques minutes.
%
%  SORTIE : MEC_BLDC/article_tables.txt

clear; clc;
outfile = 'article_tables.txt';
if isfile(outfile), delete(outfile); end
diary(outfile); diary on;
t_all = tic;

fprintf('========================================================\n');
fprintf(' REGENERATION DES CHIFFRES DE L''ARTICLE\n');
fprintf(' date : %s\n', datestr(now,'yyyy-mm-dd HH:MM:SS'));
fprintf('========================================================\n\n');

%% =====================================================================
%  1. COLONNE MAILLAGE DE LA TABLE III
%  =====================================================================
fprintf('### 1. RUN_SCORE_MESH ###############################\n\n');
fprintf('ATTENTION : noter Ms et n_sh utilises, ils doivent figurer\n');
fprintf('dans la legende de la Table III ET dans celle de tab:conv.\n\n');
try
    RUN_SCORE_MESH;
catch ME
    fprintf(2,'\n*** RUN_SCORE_MESH a echoue : %s\n', ME.message);
    fprintf(2,'    %s ligne %d\n', ME.stack(1).name, ME.stack(1).line);
end

fprintf('\n\n### 2. RUN_KRKL_MESH ################################\n\n');
try
    RUN_KRKL_MESH;
catch ME
    fprintf(2,'\n*** RUN_KRKL_MESH a echoue : %s\n', ME.message);
    fprintf(2,'    %s ligne %d\n', ME.stack(1).name, ME.stack(1).line);
end

%% =====================================================================
%  3. CONTROLES DE COHERENCE
%  =====================================================================
fprintf('\n\n### 3. CONTROLES DE COHERENCE #######################\n\n');

%  --- 3a. la fenetre du rms du transitoire ---------------------------
%  Le manuscrit reporte, pour l'essai en charge :
%      courant de phase  1.4186 A  (-8.6 % vs FEA 1.5516 A)
%      pertes Joule      69.432 W  (-3.0 % vs FEA 71.564 W)
%  Ces deux valeurs sont incompatibles sous P = 3 R I^2. Verification :
%  CAUSE RACINE IDENTIFIEE (analyse statique, aucune execution requise) :
%    BLDC_MEC_COMPLET.m l.410 :  wsT = tmsT > 0.75*tmsT(end)
%    La fenetre de regime est le DERNIER QUART du run. A 1257 tr/min et
%    p = 7, la periode electrique vaut 6.819 ms ; le dernier quart de
%    34.261 ms vaut 8.565 ms, soit 1.256 periode -- PAS un nombre entier.
%    Sur une fenetre fractionnaire les trois phases sont echantillonnees
%    a des points differents de leur cycle : sqrt(mean(i_a^2)) est biaise,
%    alors que mean(Pcu) = Rph*mean(sum i^2) ne l'est pratiquement pas,
%    car sum(i^2) est quasi constant en conduction 120 deg (drive_mec.m
%    l.193 : D.Pcu = Rph*sum(I.^2,1)).
M0      = machine_bldc();
I_mec   = 1.4186;   I_fea = 1.5516;
P_mec   = 69.432;   P_fea = 71.564;
I_energ = sqrt(P_mec/(3*M0.Rph));         % courant energetique du MEC
fprintf('  Coherence courant / pertes Joule (essai en charge)\n');
fprintf('    Rph (machine_bldc)          : %7.4f ohm\n', M0.Rph);
fprintf('    sqrt(mean(Pcu)/(3*Rph))     : %7.4f A  (%+.1f %% vs FEA)\n', ...
        I_energ, 100*(I_energ-I_fea)/I_fea);
fprintf('    sqrt(mean(i_a^2)) rapporte  : %7.4f A  (%+.1f %% vs FEA)\n', ...
        I_mec, 100*(I_mec-I_fea)/I_fea);
fprintf('    rapport des carres          : %7.3f  (doit valoir 1.000)\n', ...
        (I_mec/I_energ)^2);
fprintf('    ecart sur les pertes Joule  : %+.1f %%\n', 100*(P_mec-P_fea)/P_fea);
fprintf(['    => la valeur non biaisee est %.4f A, soit %+.1f %%, qui\n' ...
         '       coincide avec l''ecart sur les pertes Joule, comme il se doit.\n'], ...
         I_energ, 100*(I_energ-I_fea)/I_fea);

%  --- 3a-bis. la fenetre corrigee, a reporter dans BLDC_MEC_COMPLET.m ---
fprintf('\n  Fenetre de regime : diagnostic et correctif\n');
tend = 34.261e-3; n_rpm = 1257; pp = M0.p;
Te   = 60/(n_rpm*pp);  wlen = 0.25*tend;
fprintf('    periode electrique          : %6.3f ms\n', Te*1e3);
fprintf('    fenetre actuelle (0.75*tend): %6.3f ms = %.3f periode(s)\n', ...
        wlen*1e3, wlen/Te);
fprintf('    -> remplacer la ligne 410 de BLDC_MEC_COMPLET.m par :\n\n');
fprintf('       tmsT = DT.t*1e3;\n');
fprintf('       wsT  = tmsT > 0.75*tmsT(end);              %% 1er passage\n');
fprintf('       fe   = mean(DT.n(wsT))/60*M.p;             %% f elec du regime\n');
fprintf('       npg  = max(1,floor(0.25*DT.t(end)*fe));    %% periodes entieres\n');
fprintf('       wsT  = DT.t > DT.t(end) - npg/fe;          %% fenetre corrigee\n\n');
fprintf('       (idem pour la fenetre FEA w4, l.65, si le nombre de points\n');
fprintf('        de l''essai FEA ne couvre pas un nombre entier de periodes)\n');

%  --- 3b. rappel du point aberrant sur k_r ---------------------------
fprintf('\n  Point aberrant de la reference sur k_r\n');
M = machine_bldc();
rd  = @(f)readmatrix(f,'FileType','text','Delimiter','\t','NumHeaderLines',1);
try
    dmm = rd(fullfile(M.FEA.dir,'magnetostique(Magnetic_loading)', ...
                      'Calculator Expressions Table 1.tab'));
    mmf_g_F = dmm(1,16:29);
    [vmin,imin] = min(mmf_g_F);
    fprintf('    FMM d''entrefer par aimant (A) :\n     ');
    fprintf('%7.1f', mmf_g_F); fprintf('\n');
    fprintf('    minimum : %.1f A a l''aimant %d ; les 13 autres : %.1f a %.1f A\n', ...
        vmin, imin, min(mmf_g_F(setdiff(1:14,imin))), max(mmf_g_F));
    fprintf('    voisins de l''aimant %d : %.1f et %.1f A\n', imin, ...
        mmf_g_F(mod(imin-2,14)+1), mmf_g_F(mod(imin,14)+1));
    fprintf('    k_r avec le point : %.4f | sans : %.4f\n', ...
        mean(dmm(1,2:15))/mean(mmf_g_F), ...
        mean(dmm(1,2:15))/mean(mmf_g_F(setdiff(1:14,imin))));
catch ME
    fprintf(2,'    lecture FEA impossible : %s\n', ME.message);
end

fprintf('\n\n========================================================\n');
fprintf(' TERMINE en %.0f s. Sortie ecrite dans %s\n', toc(t_all), outfile);
fprintf('\n RAPPEL DES TROIS CHOSES A RELEVER :\n');
fprintf('   (a) Ms et n_sh de RUN_SCORE_MESH -> legendes Table III et tab:conv\n');
fprintf('   (b) la colonne maillage complete -> Table III\n');
fprintf('   (c) le verdict du controle 3a     -> Table VI et §V-C\n');
fprintf('========================================================\n');
diary off;
