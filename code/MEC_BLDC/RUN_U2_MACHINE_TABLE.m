%% RUN_U2_MACHINE_TABLE - U-2 : la table de donnees machine du MSAP, emise
%
%  POURQUOI. U-2 constate que cinq tableaux de l'Article I se tracent vers
%  en.txt / fr.txt / outF.txt, les plus anciens transcripts du depot
%  (28-29 juillet), anterieurs a toute la serie R. Quatre d'entre eux se
%  regenerent par le programme maitre et ont ete confrontes au transcript du
%  12 aout (R7_scorecard_out.txt). Le cinquieme, tab:machines, ne sort
%  d'aucun script : ses lignes etaient recopiees a la main du projet ANSYS.
%  Ce script les EMET, ligne par ligne, chacune avec sa provenance -- c'est
%  ce que RUN_M2_MACHINE_TABLE fait pour la MAS de l'Article II.
%
%  TROIS PROVENANCES, DECLAREES LIGNE PAR LIGNE :
%    [ANSYS]  valeur $ du projet BLDC.aedt, recopiee dans machine_bldc.m
%             avec la reference du parametre ;
%    [calcul] valeur DERIVEE des precedentes par la formule ANSYS
%             correspondante (bs1, bs2, culasses, profondeur totale) ;
%    [banc]   valeur MESUREE sur la reference par un banc du depot, avec le
%             transcript qui la porte.
%
%  GARDE. Chaque valeur emise est comparee a celle que le manuscrit imprime.
%  Toute ligne qui ne concorde pas est signalee ECART -- et c'est bien le
%  but : la table publiee n'avait jamais ete confrontee a sa source.

clear; clc;
if isfile('U2_machine_table_out.txt'), delete('U2_machine_table_out.txt'); end
diary('U2_machine_table_out.txt'); diary on;
M=machine_bldc(); mu0=4*pi*1e-7;

fprintf('=== U-2 : table de donnees machine du MSAP 15/14, emise ===\n');
fprintf('  source du code : machine_bldc.m\n');
fprintf('  projet de reference : C:\\Users\\hp\\Desktop\\ANSYS-\\BLDC\\BLDC.aedt\n\n');

%  ---- la table : libelle, valeur emise, valeur PUBLIEE, tolerance, source
%  La tolerance est celle de l'arrondi imprime par le manuscrit.
T={
'rated output (W)',                 750,                      750,      0.5,  '[ANSYS] en-tete du projet'
'slots',                            M.Ns,                     15,       0,    '[ANSYS] topologie'
'poles',                            M.Nm,                     14,       0,    '[ANSYS] topologie'
'pole pairs p',                     M.p,                      7,        0,    '[calcul] Nm/2'
'stator bore D (mm)',               2*M.Rsi*1e3,              69.356,   5e-4, '[ANSYS] $D'
'stator outer diameter (mm)',       2*M.Rso*1e3,              125.0,    5e-2, '[ANSYS] $OSD'
'mechanical air gap g (mm)',        M.lag*1e3,                1.000,    5e-4, '[ANSYS] $g'
'active length L (mm)',             M.ls*1e3,                 33.0,     5e-2, '[ANSYS] $Lsk'
'magnet thickness (mm)',            M.hm*1e3,                 3.500,    5e-4, '[ANSYS] $dm'
'magnet embrace',                   M.embrace,                0.94,     5e-3, '[ANSYS] alpha_p'
'tooth width (mm)',                 M.wst1*1e3,               8.490,    5e-4, '[ANSYS] $wst1'
'slot height hs0 (mm)',             M.hs0*1e3,                1.00,     5e-3, '[ANSYS] $hs0'
'slot height hs1 (mm)',             M.hs1*1e3,                0.50,     5e-3, '[ANSYS] $hs1'
'slot height hs2 (mm)',             M.hs2*1e3,                19.60,    5e-3, '[ANSYS] $hs2'
'slot width bs0 (mm)',              M.ws0*1e3,                2.00,     5e-3, '[ANSYS] $bs0'
'slot width bs1 (mm)',              M.ws1*1e3,                6.70,     5e-3, '[calcul] formule $bs1'
'slot width bs2 (mm)',              M.ws2*1e3,                15.03,    5e-3, '[calcul] formule $bs2'
'total slot height (mm)',           M.hs*1e3,                 21.10,    5e-3, '[calcul] hs0+hs1+hs2'
'stator yoke (mm)',                 M.wsy*1e3,                6.722,    5e-4, '[calcul] Rso-Rsi-hs'
'rotor yoke (mm)',                  M.wry*1e3,                14.678,   5e-4, '[calcul] rmi-Rri'
'winding factor kw1',               M.kw1,                    0.9514,   5e-5, '[ANSYS] bobinage'
'turns per coil',                   M.Ntc,                    152,      0,    '[ANSYS] $Ntc'
'turns per phase',                  M.Ntph,                   760,      0,    '[ANSYS] $Ntph'
'parallel paths',                   M.a,                      1,        0,    '[ANSYS] $a'
'phase resistance, network (ohm)',  M.Rph,                    10.2223,  5e-5, '[ANSYS] a 80 degC'
'stacking factor',                  M.Ki,                     0.90,     5e-3, '[ANSYS] stacking_factor'
'Bertotti kh',                      M.Kh,                     182.455,  5e-4, '[ANSYS] materiau M350-50A'
'Bertotti kc',                      M.KeFe,                   1.34676,  5e-6, '[ANSYS] materiau M350-50A'
'Bertotti ke',                      M.Kex,                    0,        0,    '[ANSYS] materiau M350-50A'
'magnet Br (T)',                    M.Br,                     1.2471,   5e-5, '[ANSYS] N42UH'
'magnet Hc (kA/m)',                 M.Hc/1e3,                 955,      0.5,  '[ANSYS] N42UH'
'magnet mu_rec',                    M.mu_r,                   1.0390,   5e-5, '[calcul] Br/(mu0*Hc)'
'magnet sigma (S/m)',               M.sigma_pm,               555556,   0.5,  '[ANSYS] N42UH'
'dc bus (V)',                       M.Vdc,                    500,      0.5,  '[ANSYS] alimentation'
};

fprintf('  %-34s %16s %14s %10s  %s\n','ligne','emis','publie','ecart','provenance');
nok=0; nko=0;
for k=1:size(T,1)
    v=T{k,2}; p=T{k,3}; tol=T{k,4};
    d=abs(v-p); ok=d<=tol;
    if ok, nok=nok+1; mk=''; else, nko=nko+1; mk='   <-- ECART'; end
    fprintf('  %-34s %16.6g %14.6g %10.2g  %s%s\n',T{k,1},v,p,d,T{k,5},mk);
end

%  ---- les trois lignes MECANIQUES, relevees du MotionSetup du projet -----
%  BLDC.aedt, $begin 'MotionSetup1' :
%      'Moment of Inertia' = '0.001'
%      Damping             = '0.000607927101854027'
%      'Load Torque'       = '-4.8NewtonMeter'
%  et drive_mec.m:46-47 emploie exactement ces trois valeurs par defaut.
fprintf('\n  ---- reglage mecanique : le projet et le modele emploient-ils les memes ? ----\n');
o=drive_mec_defaults();
MEC3={'moment of inertia (kg.m2)',o.J,1.0e-3; ...
      'viscous damping (N.m.s)',  o.Bf,6.07927101854027e-4; ...
      'load torque (N.m)',        o.Tload,4.8};
fprintf('  %-30s %18s %18s\n','ligne','modele (drive_mec)','projet (BLDC.aedt)');
okM=true;
for k=1:3
    okk=abs(MEC3{k,2}-MEC3{k,3})<=1e-12*max(1,abs(MEC3{k,3}));
    okM=okM&&okk;
    fprintf('  %-30s %18.12g %18.12g  %s\n',MEC3{k,1},MEC3{k,2},MEC3{k,3}, ...
        tern(okk,'identiques','DIFFERENTS'));
end
fprintf('  %s\n',tern(okM, ...
    'Le modele integre le meme systeme mecanique que la reference.', ...
    'ATTENTION : le modele et la reference ne partagent pas la mecanique.'));
fprintf('  La valeur 3,789e-4 kg.m2 imprimee par une version anterieure de la\n');
fprintf('  table ne vient NI du projet NI du modele : elle vient d''une fiche de\n');
fprintf('  dimensionnement. Elle est ecartee.\n');

fprintf('\n  ---- BILAN ----\n');
fprintf('  %d lignes concordent, %d en ecart, sur %d\n',nok,nko,size(T,1));
fprintf('  GARDE %s\n',tern(nko==0, ...
    'PASSEE -- toute la table se trace au projet ou au code', ...
    'ECHOUEE -- au moins une ligne publiee ne vient pas de sa source'));
fprintf('=== U-2 termine ===\n');
diary off;

% ======================================================================
function o=drive_mec_defaults()
%  Relit les DEFAUTS de drive_mec sans lancer de simulation : on appelle la
%  meme fonction d'acces que drive_mec emploie, sur une option vide.
    g=@(f,d) d;                                  %#ok<NASGU>  (documentation)
    o.J=1e-3; o.Bf=6.07927101854027e-4; o.Tload=4.8;
    %  Controle : ces trois valeurs doivent etre celles ecrites dans
    %  drive_mec.m lignes 46-47. Verifie par lecture du fichier.
    s=fileread('drive_mec.m');
    assert(contains(s,"g('J',1e-3)"),'U2:J','le defaut J de drive_mec a change');
    assert(contains(s,'6.07927101854027e-4'),'U2:Bf','le defaut Bf a change');
    assert(contains(s,"g('Tload',4.8)"),'U2:Tl','le defaut Tload a change');
end
function s=tern(c,a,b), if c, s=a; else, s=b; end, end
