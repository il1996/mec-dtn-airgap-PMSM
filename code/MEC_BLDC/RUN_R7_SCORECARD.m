%% RUN_R7_SCORECARD - R7 : recompter les vingt-neuf grandeurs de la synthese
%
%  CONSTAT DU RAPPORTEUR. Le comptage visuel de la Fig. 10 donne au moins
%  SEPT barres hors +/-5 %, alors que le manuscrit annonce 23 sur 29 dans les
%  5 %, donc six depassements.
%
%  PRECEDENT. Le comptage voisin -- "each of the two reported being the better
%  on five of the seventeen" -- s'est revele FAUX au recomptage entree par
%  entree (douze / quatre avec une egalite). Le 23/29 est donc traite comme
%  non verifie tant qu'il n'a pas ete recompte DEPUIS LES DONNEES.
%
%  METHODE. On execute le programme maitre, qui construit lui-meme la liste
%  SC des 29 grandeurs, puis on recompte a PLEINE PRECISION. Aucune valeur
%  n'est relue d'un tableau ni transcrite.
%
%  CACHE. Le programme maitre reutilise mec_map.mat / mec_map_field.mat.
%  Ces fichiers datent du 28/07 alors que machine_bldc.m a change le 04/08,
%  mais R6 a verifie qu'ils redonnent la carte reconstruite a la meme grille
%  (104,318239 / 82,244042 / 62,309367 / 12,057481 mH a 0/5/10/25 A). Le
%  comptage porte donc bien sur la configuration PUBLIEE.
%
%  GARDE. Le recomptage doit reproduire EXACTEMENT le compteur nok que le
%  programme maitre calcule lui-meme (ligne 748-753). Sinon c'est le
%  recomptage qui est faux, pas le manuscrit.
clear; clc;
fprintf('=== R7 : execution du programme maitre (comptage a suivre) ===\n');
%  BLDC_MEC_COMPLET commence par "clear" : un tic range dans une variable
%  ordinaire ne survivrait pas. On le confie donc a l'appdata de la racine,
%  que clear ne touche pas.
setappdata(0,'r7tic',tic);
BLDC_MEC_COMPLET;                       %#ok<*NASGU>  % clear + tout le programme
tR7=toc(getappdata(0,'r7tic'));

diary('R7_scorecard_out.txt'); diary on;
TEX=fullfile('..','article','ArticleI_DtN_PMSM.tex');
fprintf('\n\n=== R7 : recomptage des vingt-neuf grandeurs ===\n\n');
fprintf('  CONFIGURATION\n');
fprintf('    programme    : BLDC_MEC_COMPLET (execute a l''instant, %.0f s)\n',tR7);
fprintf('    machine      : PMSM 15/14, 750 W\n');
fprintf('    cartes       : mec_map.mat / mec_map_field.mat, verifiees en R6\n');
fprintf('    seuil        : |ecart| < 5 %% (strict, comme le code : abs(e)<5)\n');
fprintf('    manuscrit    : %s\n',TEX);

E=arrayfun(@(s)(s.m-s.f)/s.f*100,SC);
N=numel(SC);

%% ---- 1. LA LISTE, A PLEINE PRECISION ---------------------------------
fprintf('\n  ---- 1. LES %d GRANDEURS, ECART SIGNE A PLEINE PRECISION ----\n',N);
fprintf('  %3s %-34s %16s %16s %14s %s\n','#','grandeur','MEC','FEA','ecart (%)','');
for k=1:N
    fl='';
    if abs(E(k))>=5, fl='  <-- HORS 5 %'; end
    fprintf('  %3d %-34s %16.6g %16.6g %+13.4f %s\n',k,SC(k).n,SC(k).m,SC(k).f,E(k),fl);
end

%% ---- 2. GARDE : le recomptage reproduit-il le compteur du programme ? --
myok=sum(abs(E)<5);
fprintf('\n  ---- 2. GARDE ----\n');
fprintf('    compteur du programme maitre (nok) : %d sur %d\n',nok,N);
fprintf('    recomptage independant             : %d sur %d\n',myok,N);
G=(myok==nok);
fprintf('    GARDE %s\n',tern(G,'PASSEE','ECHOUEE -- le recomptage est en cause'));

%% ---- 3. CE QUE LE MANUSCRIT ANNONCE ----------------------------------
fprintf('\n  ---- 3. CONTRE L''ANNONCE DU MANUSCRIT ----\n');
tx=fileread(TEX);
tk=regexp(tx,['twenty-nine quantities computed, of which\s*\\?\s*' ...
              'twenty-three fall within \\SI\{(\d+)\}\{\\percent\}'],'tokens','once');
tk2=regexp(tx,'(\w+)-(\w+) fall within \\SI\{5\}\{\\percent\}','tokens','once');
if ~isempty(tk2)
    fprintf('    texte relu du manuscrit : "%s-%s fall within 5 %%"\n',tk2{1},tk2{2});
end
fprintf('    annonce  : 23 sur 29 dans les 5 %%, donc 6 depassements\n');
fprintf('    mesure   : %d sur %d dans les 5 %%, donc %d depassements\n',myok,N,N-myok);
if myok==23
    fprintf('    => L''ANNONCE EST EXACTE sur le compte brut.\n');
else
    fprintf('    => L''ANNONCE EST FAUSSE sur le compte brut. A CORRIGER : %d/%d.\n',myok,N);
end

%% ---- 4. LES CAS LIMITES, QUI EXPLIQUENT LE COMPTAGE VISUEL -----------
fprintf('\n  ---- 4. CAS LIMITES (le comptage visuel de la figure) ----\n');
fprintf('    La figure trace des barres et le texte arrondit a 0,1 point.\n');
fprintf('    Grandeurs a moins de 1 point du seuil :\n');
bl=find(abs(abs(E)-5)<1);
if isempty(bl)
    fprintf('      aucune -- le comptage visuel ne peut pas etre ambigu de ce fait.\n');
else
    for k=bl(:).'
        fprintf('      %-34s %+8.4f %%  -> %s\n',SC(k).n,E(k), ...
            tern(abs(E(k))<5,'DANS les 5 %','HORS les 5 %'));
    end
end

%% ---- 5. EXCLUSION DES GRANDEURS DECLAREES NON VALIDEES ---------------
%  Chaque exclusion est justifiee par une phrase DU MANUSCRIT, citee.
EXC={ ...
 'Pertes aimant a vide',      '§6.6 : "The reference cannot resolve the quantity on which it is being compared."'; ...
 'Pertes aimant en charge',   '§6.6 + Table 19 : ecart -1,6 W sous le plancher de resolution de 52,7 W'; ...
 'Couple a l''arret (balayage)','§5 : "cannot be claimed as a validation ... excluded from the comparison"'};
fprintf('\n  ---- 5. HORS GRANDEURS DECLAREES NON VALIDEES PAR LE MANUSCRIT ----\n');
keep=true(1,N);
for j=1:size(EXC,1)
    q=find(strcmp({SC.n},EXC{j,1}));
    if isempty(q)
        fprintf('    [!] "%s" introuvable dans la liste -- exclusion NON appliquee\n',EXC{j,1});
        continue;
    end
    keep(q)=false;
    fprintf('    exclue #%d %-30s %+8.4f %%\n      motif : %s\n',q,SC(q).n,E(q),EXC{j,2});
end
Nk=sum(keep); okk=sum(abs(E(keep))<5);
fprintf('\n    compte hors grandeurs non validees : %d sur %d dans les 5 %%\n',okk,Nk);
fprintf('    soit %d depassement(s) sur %d.\n',Nk-okk,Nk);
fprintf('\n    les depassements qui SUBSISTENT apres exclusion :\n');
for k=find(keep & abs(E)>=5)
    fprintf('      %-34s %+8.4f %%\n',SC(k).n,E(k));
end

%% ---- 6. RESUME ------------------------------------------------------
fprintf('\n  ---- 6. RESUME ----\n');
fprintf('    annonce du manuscrit          : 23 / 29\n');
fprintf('    compte brut mesure            : %d / %d\n',myok,N);
fprintf('    compte hors non validees      : %d / %d\n',okk,Nk);
save('R7_scorecard.mat','SC','E','keep','myok','okk','N','Nk');
fprintf('\n  duree totale %.0f s\n=== R7 termine ===\n',tR7);
diary off;

% ======================================================================
function s=tern(c,a,b), if c, s=a; else, s=b; end, end
