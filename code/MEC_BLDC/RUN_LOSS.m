%% RUN_LOSS  Pertes fer locales sur champ dente + CHAMP ROTATIONNEL (C2)
%  Les flux de dent et de culasse sont extraits du RESEAU couple sur une
%  periode ELECTRIQUE. Bertotti volumique M350-50A (formulation ANSYS) :
%      pv = Kh*f*Bpk^2 + Ke*<(dB/dt)^2>/(2*pi^2)      [W/m^3]   (Kex = 0)
%
%  CHAMP VECTORIEL DE CULASSE - decomposition a DEUX COMPOSANTES (le modele
%  employe par ANSYS Maxwell lui-meme) :   pv = pv(B_tang) + pv(B_rad)
%  La culasse porte simultanement :
%    * un flux TANGENTIEL (circulation vers les poles voisins), et
%    * le flux RADIAL que la dent y injecte (Phi_T sur un pas d'encoche).
%  Le locus mesure est QUASI RECTILIGNE a ~40 deg (composantes en phase) :
%  le champ de culasse est donc LINEAIREMENT POLARISE SELON UN AXE INCLINE,
%  et non tournant. Un modele 1D ne comptait qu'UNE composante (0.612 T)
%  alors que la RESULTANTE vaut sqrt(0.612^2+0.500^2) = 0.790 T.
%  Pour des composantes en phase la somme est EXACTE :
%      pv(B_t) + pv(B_r) = k*(B_t^2 + B_r^2) = k*|B|^2
%  Ce n'est donc PAS une majoration empirique mais la simple restitution du
%  module du vecteur B. (L'ouverture residuelle de l'ellipse traduit une
%  petite part reellement rotationnelle, du second ordre.)
%  Les DENTS restent quasi 1D (flux radial sur 19.6 des 21.1 mm de hauteur).
%
%  Reference FEA : essai transitoire A VIDE, CoreLoss = 19.73 W a 1688 tr/min
%  (grandeur INTEGRALE bien convergee ; rotor ~0 au synchronisme ; les
%  0.33 W de pertes aimant sont comptes separement par ANSYS en SolidLoss).
clear; clc;
M=machine_bldc(); Ns=M.Ns; p=M.p;
if isfile('kfringe_ident.mat'), S=load('kfringe_ident.mat'); kfr=S.kbest; else, kfr=0.75; end
f=M.FEA.n_nl*M.Nm/120;
fprintf('=== PERTES FER LOCALES + CHAMP ROTATIONNEL (C2) ===\n');
fprintf('  kfringe identifie = %.2f | f = %.1f Hz (%g tr/min)\n',kfr,f,M.FEA.n_nl);

Np=721;
R=cogging_mec(M,840,0,Np,M.muI,kfr,2*pi/p);

% ---- geometrie des elements ------------------------------------------
Vt1=M.wst1*(M.hs0+M.hs1+M.hs2)*M.ls*M.Ki;              % une dent
tau_y=2*pi*(M.Rso-M.wsy/2)/Ns;                          % pas de culasse (rayon moyen)
Vy1=tau_y*M.wsy*M.ls*M.Ki;                              % un element de culasse
Ryi=M.Rsi+M.hs0+M.hs1+M.hs2;                            % rayon interieur de culasse
tau_yi=2*pi*Ryi/Ns;                                     % pas d'encoche a l'entree de culasse
fprintf('  Element de culasse : pas %.2f mm, epaisseur %.2f mm\n',tau_y*1e3,M.wsy*1e3);
fprintf('  Surface d''entree RADIALE dans la culasse : %.2f mm x L (pas d''encoche)\n',tau_yi*1e3);

% ---- inductions locales ----------------------------------------------
Bt_rad = R.PhiT/(M.wst1*M.ls*M.Ki);                     % dent : radial
% Culasse : element i = l'arc entre la dent i et la dent i+1.
%   - tangentiel  : flux de BRANCHE Phi_Y,i (rigoureux, pas de moyenne !)
%   - radial      : flux Phi_T,i que la dent i injecte a sa frontiere gauche
%     (chaque flux de dent est ainsi affecte a UN SEUL element : pas de
%      double comptage).
%   NB : ne PAS moyenner les segments voisins. En bobinage dentaire 15/14 les
%   dents voisines portent des flux quasi OPPOSES (pas de 24 deg contre une
%   periode polaire de 51.4 deg) : la moyenne s'annule (0.05 T au lieu de
%   0.62 T) alors que |B| tangentiel est bien reel dans chaque arc.
By_tan = R.PhiY/(M.wsy*M.ls*M.Ki);                      % tangentiel
By_rad = R.PhiT/(tau_yi*M.ls*M.Ki);                     % radial (flux de dent)

tt=linspace(0,1/f,Np); dt=tt(2)-tt(1);
bert=@(B,V) deal( sum(M.Kh*f*((max(B,[],2)-min(B,[],2))/2).^2)*V, ...
                  sum((M.KeFe/(2*pi^2))*mean(gradient(B,dt).^2,2))*V );

[Pht,Pet]=bert(Bt_rad,Vt1);                             % dents (radial)
[Phy1,Pey1]=bert(By_tan,Vy1);                           % culasse tangentiel
[Phy2,Pey2]=bert(By_rad,Vy1);                           % culasse radial

P_teeth  = Pht+Pet;
P_yoke_1c= Phy1+Pey1;
P_yoke_2c= Phy1+Pey1+Phy2+Pey2;

fprintf('\n[1] Inductions cretes\n');
fprintf('    dent (radial)          : %.3f T\n',max((max(Bt_rad,[],2)-min(Bt_rad,[],2))/2));
fprintf('    culasse (tangentiel)   : %.3f T\n',max((max(By_tan,[],2)-min(By_tan,[],2))/2));
fprintf('    culasse (radial)       : %.3f T   <- ignoree par un modele 1D\n',...
        max((max(By_rad,[],2)-min(By_rad,[],2))/2));
bta=max((max(By_tan,[],2)-min(By_tan,[],2))/2);
bra=max((max(By_rad,[],2)-min(By_rad,[],2))/2);
fprintf('    RESULTANTE culasse     : %.3f T   (axe incline de %.0f deg)\n',...
        hypot(bta,bra),atan2d(bra,bta));
fprintf('    -> le modele 1D ne comptait que %.3f T sur une resultante de %.3f T :\n',bta,hypot(bta,bra));
fprintf('       il manquait un facteur %.2f sur B^2 dans la culasse.\n',(hypot(bta,bra)/bta)^2);

fprintf('\n[2] Bilan des pertes fer\n');
fprintf('    %-34s %8s %8s %8s\n','','hyst','Foucault','total');
fprintf('    %-34s %8.2f %8.2f %8.2f W\n','dents (radial)',Pht,Pet,P_teeth);
fprintf('    %-34s %8.2f %8.2f %8.2f W\n','culasse - tangentiel seul (1D)',Phy1,Pey1,P_yoke_1c);
fprintf('    %-34s %8.2f %8.2f %8.2f W\n','culasse - vecteur complet',Phy1+Phy2,Pey1+Pey2,P_yoke_2c);

fprintf('\n[3] CONFRONTATION FEA (CoreLoss = %.2f W)\n',M.FEA.Pfe_nl);
row=@(n,a) fprintf('    %-40s %8.2f W %+8.1f %%\n',n,a,(a-M.FEA.Pfe_nl)/M.FEA.Pfe_nl*100);
row('modele 1D (une seule composante)',P_teeth+P_yoke_1c);
row('modele VECTORIEL (2 composantes)',P_teeth+P_yoke_2c);
fprintf('    (pour memoire : modele a champ LISSE = 27.30 W, +38.4 %%)\n');
Ptot=P_teeth+P_yoke_2c;
fprintf('\n    Repartition finale : dents %.1f %% | culasse %.1f %%\n',...
        P_teeth/Ptot*100,P_yoke_2c/Ptot*100);
fprintf('    Residu %+.1f %% : pertes aimant (0.33 W, SolidLoss ANSYS separe),\n',...
        (Ptot-M.FEA.Pfe_nl)/M.FEA.Pfe_nl*100);
fprintf('    rotationnel du PIED DE DENT non traite (dent supposee 1D), et\n');
fprintf('    harmoniques de denture encore ~25 %% bas (identification A1).\n');

figure('Name','Pertes fer rotationnelles','Color','w','Position',[60 60 1100 420]);
subplot(1,3,1);
ang=R.phis*p*180/pi;
plot(ang,Bt_rad(1,:),'b','LineWidth',1.3); hold on;
plot(ang,By_tan(1,:),'r','LineWidth',1.3);
plot(ang,By_rad(1,:),'m--','LineWidth',1.3); grid on;
xlabel('angle elec (deg)'); ylabel('B (T)');
legend('dent (radial)','culasse (tang.)','culasse (radial)','Location','best');
title('Inductions locales');
subplot(1,3,2);
plot(By_tan(1,:),By_rad(1,:),'b','LineWidth',1.4); grid on; axis equal;
xlabel('B_{tangentiel} (T)'); ylabel('B_{radial} (T)');
title('Locus de B dans la culasse');
subplot(1,3,3);
bar([P_teeth P_yoke_1c 0; P_teeth P_yoke_2c-P_yoke_1c P_yoke_1c].','stacked'); hold on;
bb=bar([P_teeth+P_yoke_1c, P_teeth+P_yoke_2c, M.FEA.Pfe_nl]); grid on;
set(gca,'XTick',1:3,'XTickLabel',{'1D','2 comp.','FEA'}); ylabel('pertes fer (W)');
title(sprintf('%.1f -> %.1f W (FEA %.1f)',P_teeth+P_yoke_1c,Ptot,M.FEA.Pfe_nl));
saveas(gcf,fullfile(fileparts(mfilename('fullpath')),'FIG_loss.png'));
fprintf('\n  (figure FIG_loss.png sauvee)\n');
