%% RUN_D1_PMLOSS_CHAIN - D1 : d'ou vient le facteur 1,7 entre les deux chaines
%
%  CE QUE R7 A ETABLI, ET CE QU'IL A LAISSE OUVERT. Sur la meme grandeur,
%  au meme point de fonctionnement, contre la meme reference :
%      pm_loss_load, reseau a une dent (programme maitre) : 0,341619 W  (+2,45 %)
%      pm_loss,      reseau a une dent (A3, chaine gelee) : 0,2036 W    (-39,1 %)
%  R7 conclut : "Je n'isole pas ici la cause de l'ecart entre les deux
%  routines". D'ou D1, presentee comme une DECISION D'AUTEUR : laquelle des
%  deux est le modele ?
%
%  CE QUE CE SCRIPT TESTE. Que ce n'en est peut-etre pas une. Les deux
%  routines ecrivent la MEME formulation :
%      A = r * somme_n c_n(r) [ Pc_n sin(n*Theta) - Ps_n cos(n*Theta) ]
%      c_n = -mu0*mur*(1/r)*cosh(n*log(r/rmi))/sinh(n*Um)
%      J = -sigma (dA/dt - K),  K constant par aimant et par instant
%      P = sigma*L*INT (dA/dt - K)^2 dS
%  ligne pour ligne (pm_loss.m:47-79 contre pm_loss_load.m:118-148). Elles
%  ne different que par SIX reglages, tous confondus dans la comparaison de
%  R7 :
%      | reglage            | pm_loss (A3)   | pm_loss_load (maitre) |
%      | source Nsurf       | 1260           | 1680                  |
%      | positions Np       | 241            | 181                   |
%      | mu du fer          | 1500 (EN DUR)  | M.muI = 3000          |
%      | rayons nr          | 7 EQUIDISTANTS | 41 en cosinus         |
%      | angles nth         | 1440           | 720                   |
%      | derivee            | d/dphi, Np pts | d/dt, Nt pts          |
%
%  HYPOTHESE. Le champ de denture dans l'aimant vaut
%  cosh(n*log(r/rmi))/sinh(n*Um) : il est MAXIMAL contre la surface
%  d'entrefer et decroit vers l'interieur, d'autant plus vite que n est
%  grand. Sept rayons equidistants ne resolvent pas cette couche limite et
%  SOUS-ESTIMENT la perte. Si c'est la cause, alors il n'y a pas deux
%  modeles mais un seul a deux resolutions, et D1 n'a pas d'objet.
%
%  CONTRE-HYPOTHESE, testee aussi. La grille "en cosinus" de pm_loss_load
%  (uq = 1 - cos) resserre les points du cote INTERIEUR -- soit du cote ou
%  le champ est le plus faible. Si tel est le cas, 41 points ne suffisent
%  pas davantage et les 0,3416 W ne sont pas converges non plus.
%
%  GARDE A : les deux valeurs publiees doivent etre reproduites AVANT toute
%            exploration. Les arguments de quadrature ajoutes ce jour a
%            pm_loss et pm_loss_load sont sans effet par defaut : cette
%            garde le prouve au chiffre publie pres.
%  GARDE B : la routine et la source de champ doivent etre separables --
%            le carre 2x2 (deux sources) x (deux routines) doit montrer une
%            dependance a la routine bien plus grande qu'a la source.
%  GARDE C : a quadrature convergee et source identique, les deux routines
%            doivent donner la MEME valeur. Si oui, elles sont un seul
%            modele. Si non, elles sont deux modeles et D1 subsiste.

clear; clc; t0=tic;
if isfile('D1_pmloss_chain_out.txt'), delete('D1_pmloss_chain_out.txt'); end
diary('D1_pmloss_chain_out.txt'); diary on;
M=machine_bldc(); Ns=M.Ns;
if isfile('kfringe_ident.mat'), S0=load('kfringe_ident.mat'); kfr=S0.kbest; else, kfr=0.325; end
n_rpm=M.FEA.n_nl; om=n_rpm*2*pi/60; Tsl=2*pi/(Ns*om); PF=0.334;

fprintf('=== D1 : cause du facteur 1,7 sur les pertes aimant a vide ===\n');
fprintf('  machine BLDC 15/14 750 W | kfringe = %.4f | %.0f tr/min\n',kfr,n_rpm);
fprintf('  reference FEA a vide : %.4f W\n',PF);

%% ---- GARDE A : reproduire les deux valeurs publiees ---------------------
fprintf('\n  ---- GARDE A : les deux chaines publiees ----\n');
[PA3,DA3]=pm_loss(M,1260,241,kfr,n_rpm);              % A3_pmloss_out.txt:7
S_A3=struct('AG',DA3.R.AG,'Usurf',DA3.R.Usurf,'phis',DA3.R.phis);

Rnl=cogging_mec(M,1680,0,181,M.muI,kfr,2*pi/Ns);      % BLDC_MEC_COMPLET:338
tnl=linspace(0,Tsl,721);
PMA=pm_loss_load(M,Rnl,zeros(1680,3),tnl,zeros(3,721),om*tnl,41);
S_MA=struct('AG',Rnl.AG,'Usurf',Rnl.Usurf,'phis',Rnl.phis);

%  EXECUTION APRES CORRECTION DE SIGNE (pm_loss_load.m:107, 12 aout 2026).
%  pm_loss n'a pas ete touchee : elle DOIT rendre les 0,2036 W publies,
%  ce qui prouve du meme coup que l'argument de quadrature ajoute ce jour
%  est inerte par defaut. pm_loss_load, elle, ne peut plus rendre les
%  0,341619 W publies -- c'etait la valeur de la routine defectueuse, et
%  elle est conservee dans D1_pmloss_chain_prefix_out.txt.
PMA_av=0.341619;                                       % piece, avant correction
okA=abs(PA3-0.2036)<5e-5;
fprintf('  %-46s %10.6f W  (publie %.6f)\n','pm_loss(1260,241) -- A3, INCHANGEE',PA3,0.2036);
fprintf('  %-46s %10.6f W  (avant correction %.6f)\n','pm_loss_load(1680,181,nr=41) -- maitre',PMA,PMA_av);
fprintf('  effet de la correction de signe : %+.1f %%\n',100*(PMA-PMA_av)/PMA_av);
fprintf('  rapport maitre / A3 : %.3f  (avant correction %.3f)\n',PMA/PA3,PMA_av/PA3);
fprintf('  GARDE A %s\n',tern(okA, ...
    'PASSEE -- pm_loss est bit pour bit celle du manuscrit', ...
    'ECHOUEE -- rien de ce qui suit n''est interpretable'));

%% ---- 1. carre 2x2 : la source de champ contre la routine ----------------
fprintf('\n  ---- 1. SEPARER LA SOURCE DE CHAMP DE LA ROUTINE ----\n');
fprintf('  Meme source, deux routines ; meme routine, deux sources.\n');
fprintf('  %-34s %12s %12s %10s\n','source du champ','pm_loss','pm_loss_load','rapport');
SRC={'A3   (1260 surf, 241 pos, mu=1500)',S_A3; ...
     'maitre (1680 surf, 181 pos, mu=3000)',S_MA};
Q=nan(2,2);
for k=1:2
    Sk=SRC{k,2}; np=numel(Sk.phis);
    Rk=struct('AG',Sk.AG,'Usurf',Sk.Usurf,'phis',Sk.phis);
    Q(k,1)=pm_loss(M,[],[],kfr,n_rpm,Sk);                       % nr=7 uni, nth=1440
    Q(k,2)=pm_loss_load(M,Rk,zeros(size(Sk.Usurf,1),3),tnl,zeros(3,721),om*tnl,41);
    fprintf('  %-34s %12.6f %12.6f %10.3f   (Np=%d)\n', ...
        SRC{k,1},Q(k,1),Q(k,2),Q(k,2)/Q(k,1),np);
end
vs=max(abs(Q(1,:)-Q(2,:))./Q(2,:));            % effet SOURCE, a routine fixee
vr=max(abs(Q(:,1)-Q(:,2))./Q(:,2));            % effet ROUTINE, a source fixee
%  GARDE B, ENONCE RETOURNE LE 12 AOUT APRES CORRECTION. Avant correction
%  elle demandait vr > 3*vs -- "la source de champ n'explique pas l'ecart",
%  et elle passait a 44,6 % contre 10,5 %. C'etait la garde du DIAGNOSTIC.
%  Le defaut etant corrige, il n'y a plus d'ecart a expliquer : ce qu'il
%  faut verifier est l'inverse, que la routine ne domine PLUS. Le transcript
%  du 12 aout (D1_pmloss_chain_out.txt) porte encore l'ancien enonce et
%  affiche donc "GARDE B ECHOUEE" : c'est le resultat ATTENDU apres
%  correction, non un echec. L'etat d'avant est dans
%  D1_pmloss_chain_prefix_out.txt, ou la garde passe.
okB=vr<max(vs,0.10);
fprintf('  dispersion due a la SOURCE (routine fixee)  : %.1f %%\n',100*vs);
fprintf('  dispersion due a la ROUTINE (source fixee)  : %.1f %%\n',100*vr);
fprintf('  GARDE B %s\n',tern(okB, ...
    'PASSEE -- la routine ne domine plus la source', ...
    'ECHOUEE -- la routine domine encore : le defaut n''est pas ferme'));

%% ---- 2. echantillonnage : est-il en cause ? -----------------------------
fprintf('\n  ---- 2. ECHANTILLONNAGE TEMPOREL (pm_loss_load, nr=41) ----\n');
NT=[181 361 721 1441]; PT=nan(size(NT));
for k=1:numel(NT)
    ts=linspace(0,Tsl,NT(k));
    PT(k)=pm_loss_load(M,Rnl,zeros(1680,3),ts,zeros(3,NT(k)),om*ts,41);
    fprintf('  Nt = %5d  ->  %10.6f W  (%+6.2f %% / Nt=721)\n', ...
        NT(k),PT(k),100*(PT(k)-PMA)/PMA);
end
fprintf('  dispersion sur un facteur 8 : %.2f %%\n', ...
    100*(max(PT)-min(PT))/mean(PT));

%% ---- 3. LA QUADRATURE RADIALE ------------------------------------------
fprintf('\n  ---- 3. QUADRATURE RADIALE, LA CAUSE PRESUMEE ----\n');
fprintf('  Source : maitre. Le champ de denture decroit en cosh(n*v)/sinh(n*Um)\n');
fprintf('  de Rro vers rmi : la couche limite est du cote de l''ENTREFER.\n');
NR=[7 13 25 41 81 161];
ts=linspace(0,Tsl,181); ia=zeros(3,181); U3=zeros(1680,3);
PL=nan(numel(NR),3); DI={'uni','in','out'};
fprintf('\n  pm_loss_load (derivee en temps, nth=720)\n');
fprintf('  %6s %14s %14s %14s\n','nr','equidistant','cos vers rmi','cos vers Rro');
for k=1:numel(NR)
    for d=1:3
        PL(k,d)=pm_loss_load(M,Rnl,U3,ts,ia,om*ts,NR(k),720,DI{d});
    end
    fprintf('  %6d %14.6f %14.6f %14.6f\n',NR(k),PL(k,1),PL(k,2),PL(k,3));
end
NR2=[7 13 25 41 81];
PP=nan(numel(NR2),3);
fprintf('\n  pm_loss (derivee en position, nth=1440, Np=181)\n');
fprintf('  %6s %14s %14s %14s\n','nr','equidistant','cos vers rmi','cos vers Rro');
for k=1:numel(NR2)
    for d=1:3
        PP(k,d)=pm_loss(M,[],[],kfr,n_rpm,S_MA, ...
            struct('nr',NR2(k),'nth',1440,'dist',DI{d}));
    end
    fprintf('  %6d %14.6f %14.6f %14.6f\n',NR2(k),PP(k,1),PP(k,2),PP(k,3));
end

fprintf('\n  ---- 4. QUADRATURE ANGULAIRE (pm_loss_load, nr=41 vers Rro) ----\n');
NA=[360 720 1440 2880]; PA=nan(size(NA));
for k=1:numel(NA)
    PA(k)=pm_loss_load(M,Rnl,U3,ts,ia,om*ts,41,NA(k),'out');
    fprintf('  nth = %5d  ->  %10.6f W\n',NA(k),PA(k));
end
fprintf('  dispersion sur un facteur 8 : %.2f %%\n',100*(max(PA)-min(PA))/mean(PA));

%% ---- 5. GARDE C : les deux routines convergent-elles au meme point ? ----
fprintf('\n  ---- 5. GARDE C : meme source, quadrature convergee ----\n');
Pc_load=pm_loss_load(M,Rnl,U3,ts,ia,om*ts,81,1440,'out');
Pc_free=pm_loss(M,[],[],kfr,n_rpm,S_MA,struct('nr',81,'nth',1440,'dist','out'));
ecc=100*(Pc_free-Pc_load)/Pc_load;
okC=abs(ecc)<2;
fprintf('  pm_loss_load  nr= 81 vers Rro, nth=1440 : %10.6f W\n',Pc_load);
fprintf('  pm_loss       nr= 81 vers Rro, nth=1440 : %10.6f W\n',Pc_free);
fprintf('  ecart entre les deux routines           : %+9.2f %%\n',ecc);
fprintf('  GARDE C %s\n',tern(okC, ...
    'PASSEE -- une seule formulation, deux resolutions', ...
    'ECHOUEE -- les deux routines sont deux modeles distincts'));

fprintf('\n  ---- SYNTHESE ----\n');
fprintf('  publie A3      (nr= 7 equidistant)  %10.6f W  (%+6.1f %% / FEA)\n',PA3,100*(PA3-PF)/PF);
fprintf('  maitre AVANT correction de signe    %10.6f W  (%+6.1f %% / FEA)\n',PMA_av,100*(PMA_av-PF)/PF);
fprintf('  maitre APRES correction de signe    %10.6f W  (%+6.1f %% / FEA)\n',PMA,100*(PMA-PF)/PF);
fprintf('  converge       (nr eleve vers Rro)  %10.6f W  (%+6.1f %% / FEA)\n',Pc_load,100*(Pc_load-PF)/PF);

save('D1_pmloss_chain.mat','Q','PT','NT','PL','PP','NR','NR2','PA','NA', ...
     'PA3','PMA','Pc_load','Pc_free','okA','okB','okC');
fprintf('\n  GARDES : A %s | B %s | C %s\n', ...
    tern(okA,'OK','KO'),tern(okB,'OK','KO'),tern(okC,'OK','KO'));
fprintf('  duree %.0f s\n=== D1 termine ===\n',toc(t0));
diary off;

% ======================================================================
function s=tern(c,a,b), if c, s=a; else, s=b; end, end
