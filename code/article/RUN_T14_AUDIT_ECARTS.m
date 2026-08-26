%% RUN_T14_AUDIT_ECARTS - audit des colonnes d'ecart du manuscrit v2
%
%  T14. Trois cellules ont ete corrigees A LA MAIN (Table 8 nu=22, Table 15
%  Bt condense, Table 17 dent rotor). Il faut regenerer TOUTES les colonnes
%  d'ecart depuis les valeurs NON ARRONDIES.
%
%  METHODE -- et elle est le coeur de la tache. Recalculer un ecart depuis
%  les valeurs AFFICHEES est FAUX : celles-ci sont arrondies, alors que
%  l'ecart publie est (correctement) forme sur les valeurs pleines. Une
%  divergence est donc ATTENDUE et ne prouve rien.
%
%  Le test valide est un test d'INTERVALLE. Une valeur affichee v a d
%  decimales represente le vrai nombre dans [v-h, v+h] avec h = 0.5*10^-d.
%  L'ecart (a-f)/f admet alors une plage :
%      min = (a-h_a - (f+h_f)) / (f+h_f)
%      max = (a+h_a - (f-h_f)) / (f-h_f)
%  Un ecart publie DANS cette plage est compatible avec les valeurs pleines.
%  Un ecart HORS plage ne peut pas provenir des memes nombres : c'est soit
%  une correction manuelle, soit une source differente.
%
%  LIMITE A DECLARER. Ce test detecte l'incompatibilite, pas la correction
%  manuelle : une cellule corrigee a la main vers une valeur plausible reste
%  dans la plage. Il BORNE le probleme, il ne le resout pas. La regeneration
%  depuis les sources numeriques reste necessaire.
clear; clc;
diary('T14_audit_out.txt'); diary on;

%  {nom, a_affiche, f_affiche, ecart_publie_%, decimales_a, decimales_f}
T = {
'T8  ordre p=7',      1.0776, 1.0746,  +0.3, 4,4
'T8  ordre 8',        0.0181, 0.0197,  -8.1, 4,4
'T8  ordre 21',       0.2593, 0.2657,  -2.4, 4,4
'T8  ordre 22',       0.0352, 0.0328,  +7.6, 4,4
'T8  ordre 23',       0.0326, 0.0202, +61.1, 4,4
'T8  ordre 35',       0.1047, 0.1067,  -1.9, 4,4
'T8  ordre 37',       0.0427, 0.0271, +57.6, 4,4
'T10 vitesse',        1257,   1340,    -6.2, 0,0
'T10 puiss. arbre',   651.4,  683.0,   -4.6, 1,1
'T10 couple',         4.949,  4.867,   +1.7, 3,3
'T10 pertes cuivre',  69.43,  71.56,   -3.0, 2,2
'T10 courant phase',  1.419,  1.552,   -8.6, 3,3
'T10 pertes fer',     19.08,  18.47,   +3.3, 2,2
'T10 pertes aimant',  1.593,  3.183,  -50.0, 3,3
'T15 Xm en charge',   43.23,  46.00,   -6.0, 2,2
'T15 courant a vide', 9.89,   8.49,   +16.5, 2,2
'T15 fem entrefer',   375.8,  382.1,   -1.7, 1,1
'T15 couple',         115.30, 121.63,  -5.2, 2,2
'T15 Bg1',            0.909,  0.920,   -1.2, 3,3
'T15 Bt maillage',    0.103,  0.131,  -21.3, 3,3
'T15 Bt condense',    0.114,  0.131,  -12.6, 3,3
'T16 couple charge',  115.30, 121.63,  -5.2, 2,2
'T16 I1 charge',      18.95,  19.73,   -3.9, 2,2
'T16 I barre charge', 298.4,  324.7,   -8.1, 1,1
'T16 I anneau charge',1046,   1091,    -4.1, 0,0
'T16 P entree',       18809,  20369,   -7.7, 0,0
'T16 pertes fer chg', 228.3,  232.6,   -1.9, 1,1
'T16 I mag a vide',   9.89,   8.49,   +16.5, 2,2
'T16 fem a vide',     375.8,  382.1,   -1.7, 1,1
'T16 pertes fer vide',236.7,  249.3,   -5.0, 1,1
'T16 couple calage',  99.50,  104.31,  -4.6, 2,2
'T16 I1 calage',      107.87, 108.21,  -0.3, 2,2
'T16 I barre calage', 1944,   1929,    +0.7, 0,0
'T16 I anneau calage',6810,   6530,    +4.3, 0,0
'T17 culasse st vide',1.84,   1.84,    +0.0, 2,2
'T17 dent st vide',   1.94,   1.63,   +19.0, 2,2
'T17 dent rot vide',  2.12,   1.93,   +10.1, 2,2
'T17 culasse rot vide',1.67,  1.67,    -0.1, 2,2
'T17 culasse st chg', 1.72,   1.90,    -9.3, 2,2
'T17 dent st chg',    1.85,   1.68,    +9.9, 2,2
'T17 dent rot chg',   2.05,   1.87,    +9.6, 2,2
'T17 culasse rot chg',1.51,   1.58,    -4.0, 2,2
};

fprintf('=== T14 : audit d''intervalle des colonnes d''ecart (v2) ===\n');
fprintf('  Une cellule est CONFORME si l''ecart publie est atteignable par des\n');
fprintf('  valeurs pleines compatibles avec les valeurs affichees.\n\n');
fprintf('  %-22s %9s %9s %9s %9s  %s\n', ...
    'cellule','publie','min','max','naif','verdict');
nbad=0;
for k=1:size(T,1)
    a=T{k,2}; f=T{k,3}; e=T{k,4}; ha=0.5*10^-T{k,5}; hf=0.5*10^-T{k,6};
    lo=100*((a-ha)-(f+hf))/(f+hf);
    hi=100*((a+ha)-(f-hf))/(f-hf);
    naif=100*(a-f)/f;
    %  L'ECART PUBLIE EST LUI AUSSI ARRONDI (une decimale) : il represente
    %  [e-0.05, e+0.05]. Le test correct compare DEUX intervalles et ne
    %  signale que s'ils sont DISJOINTS. L'oubli de ce point produit 14 faux
    %  positifs sur 42 -- l'erreur a ete commise puis corrigee ici.
    he=0.05;
    ok = (e+he) >= lo-1e-9 && (e-he) <= hi+1e-9;
    if ~ok, nbad=nbad+1; end
    fprintf('  %-22s %+8.1f %+9.2f %+9.2f %+9.2f  %s\n', ...
        T{k,1},e,lo,hi,naif,char(9989*ok+10060*(~ok)));
end
fprintf('\n  %d cellule(s) HORS intervalle sur %d.\n',nbad,size(T,1));
if nbad==0
    fprintf(['  Aucune incompatibilite detectable. Cela NE prouve PAS que les\n' ...
             '  colonnes sont regenerees : une correction manuelle vers une\n' ...
             '  valeur plausible reste dans l''intervalle. Le test BORNE le\n' ...
             '  probleme ; il ne remplace pas la regeneration depuis les\n' ...
             '  sources numeriques.\n']);
end

%% ---- ecart entre l'ecart publie et l'ecart naif ----------------------
%  Un ecart publie tres proche de l'ecart NAIF signale au contraire une
%  cellule probablement calculee SUR LES VALEURS ARRONDIES -- le defaut
%  inverse, et lui aussi a corriger.
fprintf('\n  --- cellules ou publie ~= naif a moins de 0,01 point ---\n');
fprintf('  (indice d''un calcul fait sur les valeurs ARRONDIES)\n');
n2=0;
for k=1:size(T,1)
    naif=100*(T{k,3}*0+ (T{k,2}-T{k,3})/T{k,3});
    if abs(T{k,4}-naif)<0.01
        fprintf('    %-22s publie %+.1f  naif %+.4f\n',T{k,1},T{k,4},naif); n2=n2+1;
    end
end
fprintf('  %d cellule(s).\n',n2);
fprintf('\n=== T14 termine ===\n');
diary off;
