function M = machine_bldc()
%MACHINE_BLDC  Parametres du BLDC SPM 15 encoches / 14 poles, 750 W.
%
%   GEOMETRIE = CELLE DU MODELE ANSYS DE REFERENCE (variables $ du projet
%   C:\Users\hp\Desktop\ANSYS-\BLDC\BLDC.aedt), et NON celle de la structure
%   DES de BLDC_design.m : les deux DIFFERENT sensiblement, et c'est le
%   modele ANSYS qui produit les mesures auxquelles on se compare.
%
%      grandeur      ANSYS ($)      DES (BLDC_design)   ecart
%      alesage D     69.3564 mm     64.8645 mm          +6.9 %
%      largeur dent  8.4901 mm      7.5473 mm           +12.5 %
%      aimant dm     3.5 mm         3.2438 mm           +7.9 %
%      culasse stat. 6.7218 mm      9.2416 mm           -27 %
%      arbre Dsh     31 mm          30.50 mm            +1.6 %
%      longueur Lsk  33 mm          33.2873 mm          -0.9 %
%      OSD           125 mm         125.5494 mm         -0.4 %
%   (les largeurs d'encoche bs1/bs2, elles, sont quasi inchangees : dent et
%    alesage ont augmente de concert, dent a flancs paralleles.)
%
%   NB : $ssfr = 2 mm dans le projet ANSYS est un CONGE sur la section des
%   CONDUCTEURS (objets scain*/scaout*), PAS sur le bec de dent : le bec est
%   VIF dans le modele de reference (coherent avec DES.rfil = 0).
%
%   Machine a aimants surfaciques, rotor interne, aimantation radiale,
%   NdFeB Br=1.2471 T, tole M350-50A.
%   Sert de source unique de parametres a toute la chaine MEC_BLDC.

M = struct();

% -- topologie -----------------------------------------------------------
M.Nm = 14;  M.p = 7;  M.Ns = 15;  M.m = 3;
M.Nspp = M.Ns/(M.Nm*M.m);              % encoches / pole / phase (fractionnaire)

% -- rayons (m) -- valeurs $ du projet ANSYS ----------------------------
M.Rsi = 69.356356079956e-3/2;          % alesage stator ($D)       34.678 mm
M.Rso = 125e-3/2;                      % exterieur stator ($OSD)   62.500 mm
M.lag = 1e-3;                          % entrefer ($g)              1.000 mm
M.hm  = 3.5e-3;                        % epaisseur aimant ($dm)     3.500 mm
M.Rri = 31e-3/2;                       % arbre ($Dsh)              15.500 mm
M.Rro = M.Rsi - M.lag;                 % surface ext. aimant       33.678 mm
M.rmi = M.Rro - M.hm;                  % surface int. aimant       30.178 mm
M.rmid= M.Rsi - M.lag/2;               % mi-entrefer = $D/2-$g/2 (sondes FEA)
M.ls  = 33e-3;                         % longueur active ($Lsk)    33.000 mm

% -- encoche / dent (m) -- $hs*, $wst1, $bs0 ; bs1/bs2 par la formule ANSYS
M.hs0 = 1e-3;    M.hs1 = 0.5e-3;   M.hs2 = 19.6e-3;
M.ws0 = 2e-3;                          % ouverture d'encoche ($bs0)
M.wst1= 8.490143409079e-3;  M.wst2 = M.wst1;   % dent a flancs paralleles
M.hs  = M.hs0+M.hs1+M.hs2;             % profondeur totale d'encoche
als   = 2*pi/M.Ns;                     % $alpha_s
% $bs1 = 2*(tan(alpha_s/2)*(D/2+hs0+hs1) - wst1/(2*cos(alpha_s/2)))
M.ws1 = 2*(tan(als/2)*(M.Rsi+M.hs0+M.hs1) - M.wst1/(2*cos(als/2)));
M.ws2 = 2*(tan(als/2)*(M.Rsi+M.hs)       - M.wst1/(2*cos(als/2)));
M.wsy = M.Rso - M.Rsi - M.hs;          % culasse stator             6.722 mm
M.wry = M.rmi - M.Rri;                 % culasse rotor             14.678 mm

% -- aimant -- VALEURS EXACTES DU PROJET ANSYS ---------------------------
%  Materiau 'Arnold_Magnetics_N42UH_60C' (+ '..._South') du BLDC.aedt :
%     BHCoordinates = (-955204, 0) ... (0, 1.2471)  -> droite de recul
%     magnetic_coercivity Magnitude = -955204 A/m
%     conductivity = 555555.5556 S/m
%  => mu_r = Br/(mu0*Hc) = 1.2471/(mu0*955204) = 1.0390  (et NON 0.9872)
%  Les anciennes valeurs (Hc=1005.3 kA/m, mu_r=0.9872) venaient de la fiche
%  de DIMENSIONNEMENT (Magnet_Grade.m), pas du modele de reference : elles
%  sous-estimaient la reluctance d'aimant de 5 %, d'ou un Bg1 trop fort.
M.Br    = 1.2471;                      % remanence (T)                  [ANSYS]
M.Hc    = 955204;                      % coercivite (A/m)               [ANSYS]
M.mu_r  = M.Br/(4*pi*1e-7*M.Hc);       % permeabilite de recul = 1.0390 [ANSYS]
M.embrace = 0.94;                      % arc aimant / pas polaire (alpha_p)
M.sigma_pm= 555555.5556;               % conductivite NdFeB (S/m)       [ANSYS]
M.rho_pm  = 1/M.sigma_pm;              %   = 1.80e-6 Ohm.m (etait 1.4e-6)

% -- bobinage ------------------------------------------------------------
M.Ntc = 152;   M.Ntph = 760;   M.a = 1;
M.kw1 = 0.9514;                        % facteur de bobinage fondamental
M.Rph = 10.222304;                     % resistance de phase (ohm) a 80 degC

% -- permeabilite relative du fer (modeles LINEAIRES) --------------------
%  Valeur PHYSIQUE et non arbitraire : a l'induction de dent constatee
%  (1.34 T a vide), la courbe B(H) donne mu_r = B/(mu0*H) = 3173.
%  L'ancienne valeur 1500 etait un choix par defaut, trop RELUCTANT : elle
%  amputait l'harmonique de denture a8 de ~12 points (a8 croit fortement
%  avec mu_r ; cf. RUN_A8). Les solveurs NON LINEAIRES (onload_mec) n'en
%  ont pas besoin.
M.muI = 3000;

% -- materiau fer : Bertotti VOLUMIQUE M350-50A (comme ANSYS) ------------
%    pv = Kh*f*B^2 + Ke*(f*B)^2 + Kex*(f*B)^1.5   [W/m^3]
%  VERIFIE mot pour mot dans BLDC.aedt (materiau 'Cogent Power - M350-50A') :
%    core_loss_kh=182.455, core_loss_kc=1.34676, core_loss_ke=0,
%    mass_density=7650, stacking_factor=0.9, conductivity=2380950 S/m.
%  La courbe B(H) reelle (81 points) est desormais dans bh_curve.m ('M350').
M.Kh = 182.455;  M.KeFe = 1.34676;  M.Kex = 0;  M.rof = 7650;  M.Ki = 0.9;
M.sigma_fe = 2380950;                  % conductivite tole (S/m)        [ANSYS]
%  Affectation ANSYS : 'stator' et 'rotor' = M350 FEUILLETEE (pas de courants
%  de Foucault massifs en 2D) ; 'shaft' = steel_stainless (sigma=1.1e6, MASSIF,
%  donc compte dans SolidLoss) ; 'm1'..'m14' = N42UH (massifs).
%  => SolidLoss(FEA) = pertes AIMANTS + pertes ARBRE. L'arbre (r<=15.5 mm) est
%  ecrante par 14.7 mm de culasse rotor : sa contribution est negligeable,
%  donc SolidLoss ~ pertes aimants. (Question ouverte du README : TRANCHEE.)
M.sigma_shaft = 1.1e6;

% -- points de fonctionnement --------------------------------------------
M.speed = 1500;  M.Vdc = 500;  M.fs = 175;   % f_elec = speed*Nm/120
M.Iph_rms_load = 1.55;                 % courant de phase FEA en charge (A eff)

% -- cibles FEA (extraites des .tab ANSYS, dossier 'ANSYS resultat 750W') -
M.FEA.dir      = 'C:\Users\hp\Desktop\ANSYS résultat 750W';
M.FEA.Bg_mean  = 0.761;     % induction entrefer moyenne |Br|  (T)
M.FEA.Bg_peak  = 0.993;     % induction entrefer crete         (T)
M.FEA.emf_env  = 459.9;     % enveloppe FEM six-step crete      (V) @1500
M.FEA.emf_ph   = 233;       % FEM de phase crete (plateau)      (V)
M.FEA.FLa_pk   = 0.2583;    % flux totalise de phase crete      (Wb)
M.FEA.cog_pp   = 20.06;     % couple de detente c-c (sous-ech.) (mN.m)
M.FEA.Ld       = 52.34e-3;  % inductance synchrone              (H)
M.FEA.La       = 50.21e-3;  % inductance propre de phase        (H)
M.FEA.Pfe_nl   = 19.73;     % pertes fer a vide @1688 tr/min    (W)
M.FEA.Psolid_nl= 0.333;     % pertes aimant a vide              (W)
M.FEA.n_nl     = 1688;      % vitesse essai a vide              (tr/min)
M.FEA.T_load   = 4.87;      % couple electromagnetique en charge(N.m)
end
