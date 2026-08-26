%% RUN_C4_ECARTS - bloc C4 (SPEC_CLAUDE_CODE_v3 §8, priorite 4)
%
%  L'EXIGENCE. « Regenerer TOUTES les colonnes d'ecart depuis les sorties a
%  pleine precision, en commencant par les huit cellules signalees par
%  l'indice inverse. » T14 avait borne le probleme sans le resoudre : son
%  second bloc ne trouve aucune incompatibilite, mais il le dit lui-meme --
%  « une correction manuelle vers une valeur plausible reste dans
%  l'intervalle. Le test BORNE le probleme ; il ne remplace pas la
%  regeneration. »
%
%  CE QUE FAIT CE BLOC. Pour chaque ecart destine aux deux articles, il
%  calcule DEUX nombres :
%    - PLEIN   : l'ecart forme sur les valeurs a pleine precision ;
%    - ARRONDI : le meme ecart forme sur les valeurs telles qu'elles sont
%                AFFICHEES dans le tableau (donc arrondies).
%  Un ecart dont les deux versions different de plus de 0,05 point est
%  SENSIBLE a l'arrondi : il doit imperativement etre copie depuis la
%  colonne PLEIN.
%
%  C'est le test de l'indice inverse de T14 mene a l'endroit : au lieu de
%  chercher si un chiffre publie TRAHIT un calcul sur arrondis, on FOURNIT
%  la valeur pleine et on mesure ce que l'arrondi couterait.
%
%  SOURCES. Les fichiers .mat des blocs, qui portent les valeurs pleines.
%  Aucun chiffre n'est saisi ici : tout est charge.
clear; clc; t0=tic;
diary('C4_ecarts_out.txt'); diary on;
B=fileparts(mfilename('fullpath')); if isempty(B), B=pwd; end
I=fullfile(fileparts(B),'MEC_IM');

fprintf('=== C4 : regeneration des colonnes d''ecart a pleine precision ===\n');
fprintf('  sources : fichiers .mat des blocs, MEC_BLDC et MEC_IM\n');
fprintf('  critere : |plein - arrondi| > 0.05 point  ->  cellule SENSIBLE\n\n');

%  L : {section, libelle, valeur, reference, dec. valeur, dec. reference}
L={};

%% ---- 1. les huit cellules signalees par T14 --------------------------
fprintf('  ---- 1. les huit cellules de l''indice inverse (T14) ----\n');
fprintf(['  T14 les signalait parce que l''ecart publie coincidait avec un\n' ...
         '  calcul sur valeurs arrondies a moins de 0,01 point. Statut\n' ...
         '  aujourd''hui, apres les blocs de regeneration :\n\n']);
c8={ 'T8  ordre 21',      'tab:spectrum, Article I',      'a remplir depuis la chaine du spectre'
     'T10 vitesse',       'transitoire en charge, Art. I','regenere par BLDC_MEC_COMPLET (6 aout)'
     'T10 pertes fer',    'transitoire en charge, Art. I','SUPERSEDE par X2 et A4'
     'T15 couple',        'ancien tableau MAS',           'SUPERSEDE par B1 et B3'
     'T15 Bg1',           'ancien tableau MAS',           'SUPERSEDE par B1 sec.4'
     'T16 couple charge', 'tab:im_tests, Article II',     'SUPERSEDE par B1 sec.3'
     'T16 I barre charge','tab:im_tests, Article II',     'SUPERSEDE par B1 sec.3 et B6'
     'T17 culasse st vide','Table 17, Article II',        'SUPERSEDE par B7' };
fprintf('  %-22s %-30s %s\n','cellule','destination','statut');
for k=1:size(c8,1), fprintf('  %-22s %-30s %s\n',c8{k,1},c8{k,2},c8{k,3}); end
fprintf(['\n  SIX des huit sont deja produites par une chaine qui forme ses\n' ...
         '  ecarts sur les valeurs pleines. Les deux autres appartiennent a\n' ...
         '  des tableaux encore vides. L''indice inverse de T14 ne survit donc\n' ...
         '  dans aucune cellule destinee a la publication.\n']);

%% ---- 2. Article I ------------------------------------------------------
S=load(fullfile(B,'X1_table5b.mat'));
ib=find(S.NSH==S.nb0);
for k=1:numel(S.Mss)
    L(end+1,:)={'Art. I / tab:convb',sprintf('nu=8 lin., Ms=%d',S.Mss(k)), ...
        S.R8L(ib,k),S.r8F,5,5}; %#ok<*SAGROW>
    L(end+1,:)={'Art. I / tab:convb',sprintf('nu=8 Newton, Ms=%d',S.Mss(k)), ...
        S.R8N(ib,k),S.r8F,5,5};
end

A=load(fullfile(B,'A5bis_mutual.mat'));
L(end+1,:)={'Art. I / tab:pmsm','L_a maillage n_sh=1',A.LA(1),A.LaF,5,5};
L(end+1,:)={'Art. I / tab:pmsm','M   maillage n_sh=1',A.MU(1),A.MF,5,5};
L(end+1,:)={'Art. I / tab:pmsm','L_d maillage n_sh=1',A.LD(1),A.LdF,5,5};
L(end+1,:)={'Art. I / A5bis','L_d si M exact',A.LA(1)-A.MF,A.LdF,5,5};
L(end+1,:)={'Art. I / tab:pmsm','L_a localise',A.LA(3),A.LaF,5,5};
L(end+1,:)={'Art. I / tab:pmsm','M   localise',A.MU(3),A.MF,5,5};
L(end+1,:)={'Art. I / tab:pmsm','L_d localise',A.LD(3),A.LdF,5,5};

F=load(fullfile(B,'A4_ironloss.mat'));
kk=find(~isnan(F.RES(:,5)));
for q=kk(:).'
    L(end+1,:)={'Art. I / sec.8.5',sprintf('P_fe maillage+rotor, n_sh=%d',F.RES(q,1)), ...
        F.RES(q,5)+F.Lrot,F.PfeF,3,2};
end

%% ---- 3. Article II -----------------------------------------------------
P=load(fullfile(I,'B7_probes.mat'));
for j=1:4
    L(end+1,:)={'Art. II / Table 17',sprintf('%s, a vide',P.FEA.lbl{j}), ...
        P.P0m(1,j),P.FEA.vide(j),4,3};
end
for j=1:4
    L(end+1,:)={'Art. II / Table 17',sprintf('%s, en charge',P.FEA.lbl{j}), ...
        P.P1m(1,j),P.FEA.chg(j),4,3};
end

R=load(fullfile(I,'B6_anneau.mat'));
L(end+1,:)={'Art. II / Table 16','rapport reference en charge',R.Ir_F/R.Ib_F,R.kid,4,6};
L(end+1,:)={'Art. II / Table 16','rapport reference calage',R.Ir_Fb/R.Ib_Fb,R.kid,4,6};

K=load(fullfile(I,'B4_barskin.mat'));
L(end+1,:)={'Art. II / tab:im_ec','R''_r calage / ref. 104,31',K.RES(end,7),K.Rr_from_lit,5,5};
L(end+1,:)={'Art. II / tab:im_ec','R''_r calage / ref. 98,98',K.RES(end,7),K.Rr_from_ref,5,5};

%% ---- 4. impression et comptage, en UNE passe --------------------------
fprintf('\n  ---- 2. ecarts regeneres a pleine precision ----\n');
fprintf('  %-20s %-30s %11s %11s %9s %9s %8s\n', ...
    'destination','grandeur','valeur','reference','PLEIN %','ARRONDI %','ecart pt');
nsens=0; sect='';
for k=1:size(L,1)
    a=L{k,3}; b=L{k,4}; da=L{k,5}; db=L{k,6};
    ep=100*(a-b)/b;
    ea=100*(round(a,da)-round(b,db))/round(b,db);
    d=abs(ep-ea);
    if d>0.05, mk='  <-- SENSIBLE'; nsens=nsens+1; else, mk=''; end
    if ~strcmp(sect,L{k,1}), sect=L{k,1}; s1=sect; else, s1=''; end
    fprintf('  %-20s %-30s %11.5f %11.5f %8.3f %9.3f %8.3f%s\n', ...
        s1,L{k,2},a,b,ep,ea,d,mk);
end

%% ---- 5. .mat dont les variables restent a cartographier ---------------
fprintf('\n  ---- 3. fichiers .mat non encore exploites ----\n');
lst={fullfile(I,'B1_im_p1.mat'),fullfile(I,'B2_kC.mat'), ...
     fullfile(I,'B3_prix.mat'),fullfile(B,'A1_table7.mat'), ...
     fullfile(B,'A5_lumped.mat'),fullfile(B,'V1_pmsm.mat')};
for k=1:numel(lst)
    if isfile(lst{k})
        w=load(lst{k}); [~,nm]=fileparts(lst{k});
        fprintf('  %-16s : %s\n',nm,strjoin(fieldnames(w)',', '));
    end
end
fprintf(['\n  Ces fichiers portent les grandeurs brutes, non les colonnes\n' ...
         '  d''ecart : leurs ecarts sont deja formes sur les valeurs pleines\n' ...
         '  dans leur propre .txt. Ils sont listes pour que la verification\n' ...
         '  finale des tableaux puisse s''y adosser.\n']);

%% ---- 6. verdict -------------------------------------------------------
fprintf('\n  ---- 4. verdict ----\n');
fprintf('  %d ecarts regeneres a pleine precision ; %d sensibles a l''arrondi.\n', ...
    size(L,1),nsens);
if nsens==0
    fprintf(['  Aucune cellule ne change de plus de 0,05 point selon qu''on\n' ...
             '  forme l''ecart sur les valeurs pleines ou sur les valeurs\n' ...
             '  affichees. Cela ne dispense PAS de copier depuis la colonne\n' ...
             '  PLEIN -- c''est la regle -- mais cela borne le risque residuel\n' ...
             '  sur les tableaux deja remplis.\n']);
else
    fprintf(['  Les cellules marquees SENSIBLE doivent imperativement etre\n' ...
             '  copiees depuis la colonne PLEIN : les former sur les valeurs\n' ...
             '  affichees introduirait une erreur superieure a 0,05 point.\n']);
end
fprintf(['\n  REGLE POUR LA SUITE. Tout tableau des deux articles se remplit\n' ...
         '  depuis un .txt de bloc, jamais depuis un autre tableau. Les\n' ...
         '  colonnes d''ecart ne sont jamais recalculees a la main a partir\n' ...
         '  de colonnes de valeurs deja arrondies.\n']);

save('C4_ecarts.mat','L','nsens');
fprintf('\n  duree %.0f s\n=== C4 termine ===\n',toc(t0));
diary off;
