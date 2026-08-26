%% RUN_T9_TEMPS - comparaison de temps MEC vs elements finis
%
%  T9. La §3.6 annonce une economie de calcul sans jamais la mesurer contre
%  la reference. On extrait ici les temps ANSYS des fichiers .profile des
%  deux projets, et on les confronte aux temps MATLAB mesures.
%
%  FORMAT. Chaque enregistrement 'p(t0, t1, ...)' d'un profil Maxwell porte
%  un instant de debut et un instant de fin. L'unite n'est pas documentee :
%  on la CALIBRE en exigeant que la duree totale soit physiquement plausible,
%  et on rapporte les deux hypotheses plutot que d'en imposer une.
clear; clc;
diary('T9_temps_out.txt'); diary on;

P = { ...
 'BLDC  opti0_2',  'C:\Users\hp\Desktop\ANSYS-\BLDC\BLDC.aedtresults\opti0_2.profile'
 'BLDC  opti44_0', 'C:\Users\hp\Desktop\ANSYS-\BLDC\BLDC.aedtresults\opti44_0.profile'
 'BLDC  opti53_0', 'C:\Users\hp\Desktop\ANSYS-\BLDC\BLDC.aedtresults\opti53_0.profile'
 'MAS   opti59_5', 'C:\Users\hp\Desktop\ANSYS-\moteur_18_5\IM_18kW_690V.aedtresults\opti59_5.profile'
};

fprintf('=== T9 : temps de calcul, elements finis vs MEC ===\n\n');
fprintf('  %-16s %8s %14s %14s %14s\n', ...
    'profil','records','somme brute','si 1e-6 s','si 1e-7 s');
tot=zeros(size(P,1),1);
for k=1:size(P,1)
    f=P{k,2};
    if ~isfile(f), fprintf('  %-16s  ABSENT : %s\n',P{k,1},f); continue; end
    txt=fileread(f);
    %  p(t0, t1, ...) : on ne retient que les deux premiers nombres
    m=regexp(txt,'\n\s*p\(([0-9.eE+-]+),\s*([0-9.eE+-]+),','tokens');
    if isempty(m), fprintf('  %-16s  aucun record p()\n',P{k,1}); continue; end
    t0=cellfun(@(c)str2double(c{1}),m);
    t1=cellfun(@(c)str2double(c{2}),m);
    d=t1-t0; d=d(d>0);
    tot(k)=sum(d);
    fprintf('  %-16s %8d %14.4e %11.1f s %11.1f s\n', ...
        P{k,1},numel(d),sum(d),sum(d)*1e-6,sum(d)*1e-7);
end

%% ---- calibration de l'unite -----------------------------------------
%  Un point de solution magnetostatique 2-D de cette taille ne peut pas
%  durer 0,2 s sur une machine de bureau de 2024 ; il ne peut pas non plus
%  durer 20 s pour les cas les plus legers. On retient l'hypothese qui
%  place la duree par variation dans une plage plausible et on le DIT.
fprintf('\n  --- calibration ---\n');
for k=1:size(P,1)
    if tot(k)==0, continue; end
    f=P{k,2}; txt=fileread(f);
    m=regexp(txt,'\n\s*p\(([0-9.eE+-]+),\s*([0-9.eE+-]+),','tokens');
    n=numel(m);
    fprintf('  %-16s %4d variations -> %6.2f s/var (1e-6) | %6.3f s/var (1e-7)\n', ...
        P{k,1},n,tot(k)*1e-6/n,tot(k)*1e-7/n);
end

%% ---- confrontation aux temps MATLAB ---------------------------------
tBLDC_mec = 133;      % s, mesure par l'utilisateur
tIM_mec   = 647;      % s, mesure par l'utilisateur
fprintf('\n  --- temps MEC mesures (MATLAB) ---\n');
fprintf('  BLDC 15/14  : %d s   (8 analyses, 29 grandeurs, 9 figures)\n',tBLDC_mec);
fprintf('  MAS 48/44   : %d s   (30 glissements, 3 conditions, cartes)\n',tIM_mec);

fprintf(['\n  --- RESERVE METHODOLOGIQUE, a lire avant tout ratio ---\n' ...
  '  T9 demande une comparaison CONTROLEE : une machine, un point de\n' ...
  '  fonctionnement, MEME MATERIEL, meme cible de precision. Les conditions\n' ...
  '  ci-dessus ne la satisfont pas :\n' ...
  '    (a) les runs EF datent de 2024 sur DESKTOP-1FRG5FP, les runs MEC\n' ...
  '        sont recents et sur une autre machine ;\n' ...
  '    (b) les profils couvrent des BALAYAGES optimetrics, le MEC couvre\n' ...
  '        une etude complete : les PERIMETRES different ;\n' ...
  '    (c) aucune cible de precision commune n''a ete fixee.\n' ...
  '  Un ratio tire de ces nombres serait une illustration, PAS une mesure.\n' ...
  '  Il peut etre publie comme ordre de grandeur explicitement date et\n' ...
  '  situe ; il ne peut pas etayer une affirmation d''economie.\n']);
fprintf('\n=== T9 termine ===\n');
diary off;
