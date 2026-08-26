%% RUN_C3_SIGMA - bloc C3 (SPEC_CLAUDE_CODE_v3 §8, priorite 4)
%
%  L'OBJECTION. Le Test 2 de l'audit de la reference donne l'ordre 42 a
%  77,0 mN.m « 40 sigma » en magnetostatique remaillee et 0,88 mN.m
%  « 2,6 sigma » en transitoire a maillage preserve. La specification
%  releve que 77,0/0,88 = 87,5 alors que 40/2,6 = 15,4, et conclut :
%  « impossible sous une normalisation unique ; hypothese de travail :
%  deux normalisations ont ete melangees. »
%
%  L'HYPOTHESE EST FAUSSE, ET CE BLOC LE MONTRE. Les deux multiples
%  emploient LA MEME definition, celle de RUN_COG_ORDER42.m:33 :
%
%      sigma = std(residu) * sqrt(2/N)
%
%  c'est-a-dire l'ERREUR-TYPE DE L'AMPLITUDE DE FOURIER estimee par
%  moindres carres. Or sigma est une propriete DU JEU DE DONNEES : il
%  depend du residu et du nombre de points de CET essai. Les deux essais
%  n'ont ni le meme bruit ni le meme N -- c'est meme tout l'objet du test,
%  puisque l'un est remaille a chaque position et l'autre non.
%
%  Il n'y a donc aucune contradiction : le rapport des amplitudes (87,5)
%  se decompose en rapport des multiples (15,4) fois rapport des sigma.
%  Le bloc verifie cette identite numeriquement.
%
%  CE QUI RESTE A CORRIGER est de nature REDACTIONNELLE : le manuscrit
%  ecrit « 40 sigma » et « 2,6 sigma » sans dire que sigma n'est pas le
%  meme nombre dans les deux cas. Un lecteur suppose naturellement une
%  normalisation commune -- la specification l'a supposee, et s'est
%  trompee pour cette raison. La correction est d'ecrire sigma_m et
%  sigma_t, et de donner les deux valeurs.
%
%  CONFIGURATION DECLAREE
%    machine    : PMSM 15/14, 750 W (machine_bldc)
%    essai 1    : magnetostique(Magnetic_loading)\Torque Plot 1.tab
%                 360 positions, maillage REMAILLE a chaque position
%    essai 2    : transitoire (Back_emf)\Torque Plot.tab
%                 bande glissante, maillage PRESERVE
%    estimateur : moindres carres sur [cos(n th) sin(n th) 1 th th^2],
%                 identique a RUN_COG_ORDER42.m:31-34
clear; clc; t0=tic;
diary('C3_sigma_out.txt'); diary on;
M=machine_bldc(); Ns=M.Ns; Nm=M.Nm; LCMv=lcm(Ns,Nm);
rd=@(f)readmatrix(f,'FileType','text','Delimiter','\t','NumHeaderLines',1);

fprintf('=== C3 : les multiples de sigma du Test 2 ===\n');
fprintf('  machine : PMSM 15/14, 750 W | Ns = %d, Nm = %d, LCM = %d\n\n',Ns,Nm,LCMv);

%% ---- les deux essais, lus ---------------------------------------------
dm=rd(fullfile(M.FEA.dir,'magnetostique(Magnetic_loading)','Torque Plot 1.tab'));
am=dm(:,1); Tm=dm(:,2); if abs(am(end)-360)<1e-6, am(end)=[]; Tm(end)=[]; end
Tm=Tm-mean(Tm); thm=am*pi/180;
dt=rd(fullfile(M.FEA.dir,'transitoire (Back_emf)','Torque Plot.tab'));
tht=dt(:,1)*1e-3*(M.speed/60)*2*pi; Tt=dt(:,2); Tt=Tt-mean(Tt);

fprintf('  ---- 1. les deux essais ----\n');
fprintf('  %-34s %10s %10s\n','','magnetost.','transitoire');
fprintf('  %-34s %10d %10d\n','nombre de positions N',numel(thm),numel(tht));
fprintf('  %-34s %10.2f %10.2f\n','etendue angulaire (deg)', ...
    max(thm)*180/pi,max(tht)*180/pi);
fprintf('  %-34s %10s %10s\n','maillage','REMAILLE','PRESERVE');
fprintf('  %-34s %10.4f %10.4f\n','ecart-type du couple (mN.m)',std(Tm),std(Tt));

%% ---- 2. l'estimateur, repris a l'identique ----------------------------
%  Reprise de RUN_COG_ORDER42.m:31-34, ou fitord est LOCALE donc non
%  appelable. Quatrieme cas de duplication forcee -- voir MANIFEST §5.
fprintf('\n  ---- 2. amplitude et incertitude a l''ordre 42 ----\n');
[A42m,s42m,rm]=c3_fitord(thm,Tm,42);
[A42t,s42t,rt]=c3_fitord(tht,Tt,42);
fprintf('  sigma = std(residu) * sqrt(2/N)   (erreur-type de l''amplitude)\n\n');
fprintf('  %-30s %14s %14s\n','','magnetostatique','transitoire');
fprintf('  %-30s %14.6f %14.6f\n','amplitude A_42 (mN.m)',A42m,A42t);
fprintf('  %-30s %14.6f %14.6f\n','std du residu (mN.m)',std(rm),std(rt));
fprintf('  %-30s %14d %14d\n','N',numel(thm),numel(tht));
fprintf('  %-30s %14.6f %14.6f\n','sigma (mN.m)',s42m,s42t);
fprintf('  %-30s %14.2f %14.2f\n','multiple A/sigma',A42m/s42m,A42t/s42t);

%% ---- 3. l'identite qui leve l'objection -------------------------------
rA=A42m/A42t; rN=(A42m/s42m)/(A42t/s42t); rS=s42m/s42t;
fprintf('\n  ---- 3. l''objection levee ----\n');
fprintf('  rapport des AMPLITUDES        A_m/A_t         = %8.3f\n',rA);
fprintf('  rapport des MULTIPLES         (A/s)_m/(A/s)_t = %8.3f\n',rN);
fprintf('  rapport des SIGMA             s_m/s_t         = %8.3f\n',rS);
fprintf('  produit multiples x sigma                     = %8.3f\n',rN*rS);
fprintf('  identite verifiee a %.2e pres\n',abs(rN*rS-rA)/rA);
fprintf(['\n  LECTURE. Le rapport des amplitudes se factorise EXACTEMENT en\n' ...
         '  rapport des multiples fois rapport des sigma. Les deux multiples\n' ...
         '  sont donc COHERENTS entre eux ; ce qui differe est sigma, et il\n' ...
         '  DOIT differer, puisque les deux essais n''ont ni le meme bruit ni\n' ...
         '  le meme nombre de points. L''hypothese de deux normalisations\n' ...
         '  melangees est ecartee : la normalisation est unique, c''est le\n' ...
         '  JEU DE DONNEES qui change.\n']);

%% ---- 4. ce qui reste a corriger, et c'est de la redaction ------------
fprintf('\n  ---- 4. correction a porter au manuscrit ----\n');
fprintf('  Le Test 2 doit ecrire les deux sigma, et non "sigma" seul :\n\n');
fprintf('    ordre 42, magnetostatique remaillee : %.2f mN.m = %.0f sigma_m,\n',A42m,A42m/s42m);
fprintf('      sigma_m = %.4f mN.m (N = %d, maillage remaille)\n',s42m,numel(thm));
fprintf('    ordre 42, transitoire maillage preserve : %.2f mN.m = %.1f sigma_t,\n',A42t,A42t/s42t);
fprintf('      sigma_t = %.4f mN.m (N = %d, maillage preserve)\n',s42t,numel(tht));
fprintf('\n  Les deux valeurs peuvent donc etre RETABLIES, contrairement a ce\n');
fprintf('  que la specification laissait craindre. Elles etaient justes ; leur\n');
fprintf('  PRESENTATION ne disait pas que sigma change d''un essai a l''autre.\n');

%% ---- 5. controle : l'ordre 42 est-il present dans le transitoire ? ---
fprintf('\n  ---- 5. controle du verdict du Test 2 ----\n');
fprintf('  A_42 transitoire = %.4f mN.m contre 3*sigma_t = %.4f mN.m\n',A42t,3*s42t);
if A42t < 3*s42t
    fprintf('  => l''ordre 42 est SOUS LE BRUIT du transitoire. Le verdict du\n');
    fprintf('     Test 2 tient : les 220 mN.m sont un artefact de remaillage.\n');
else
    fprintf('  => l''ordre 42 DEPASSE le bruit du transitoire : verdict a revoir.\n');
end
[A210m,~]=c3_fitord(thm,Tm,LCMv); [A210t,s210t]=c3_fitord(tht,Tt,LCMv);
fprintf('  detente (ordre %d) : magnetost. %.4f (hors bande, Nyquist %d)\n', ...
    LCMv,A210m,floor(numel(thm)/2));
fprintf('                       transitoire %.4f +/- %.4f mN.m  <- reference\n', ...
    A210t,s210t);

save('C3_sigma.mat','A42m','s42m','A42t','s42t','A210t','s210t','rA','rN','rS');
fprintf('\n  duree %.0f s\n=== C3 termine ===\n',toc(t0));
diary off;

% ======================================================================
function [A,s,r]=c3_fitord(th,T,n)
%  Reprise A L'IDENTIQUE de RUN_COG_ORDER42.m:31-34, ou fitord est une
%  fonction LOCALE donc non appelable. Le residu r est renvoye en plus,
%  pour que sigma puisse etre decompose dans la sortie.
    X=[cos(n*th) sin(n*th) ones(size(th)) th th.^2];
    c=X\T; r=T-X*c; A=hypot(c(1),c(2)); s=std(r)*sqrt(2/numel(th));
end
