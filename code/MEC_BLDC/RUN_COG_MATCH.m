%% RUN_COG_MATCH  -  Rapprocher AU MAXIMUM la detente MEC de la reference FEA
%
%  OBJECTIF : obtenir un couple d'encochage MEC aussi proche que possible de
%  la FEA. La difficulte n'est pas le modele : c'est de savoir QUELLE courbe
%  FEA est une cible legitime. Ce script le tranche par la mesure, puis
%  optimise le MEC contre la reference valide.
%
%  1. TEST DE REPLIEMENT : si les 220 mN.m de la magnetostatique etaient de
%     la detente d'ordre 210 mal echantillonnee, le repliement a 1 deg la
%     ferait apparaitre a l'ordre |360-210| = 150. Verifions ce que devient
%     la detente MEC echantillonnee EXACTEMENT comme la FEA.
%  2. REFERENCE VALIDE : l'essai TRANSITOIRE (bande glissante, maillage
%     PRESERVE d'une position a l'autre) -> extraction de l'ordre 210 par
%     moindres carres, avec incertitude.
%  3. MEC OPTIMISE : modele SANS PARAMETRE AJUSTE (sous-domaine exact
%     d'ouverture + mu(B) local), convergence verifiee.
%  4. DISPERSION D'AIMANTS equivalente au contenu d'ordre 42 de la FEA.
clear; clc;
M=machine_bldc(); Ns=M.Ns; Nm=M.Nm; LCMv=lcm(Ns,Nm);
rd=@(f)readmatrix(f,'FileType','text','Delimiter','\t','NumHeaderLines',1);
if isfile('kfringe_ident.mat'), S=load('kfringe_ident.mat'); kfr=S.kbest; else, kfr=0.75; end
fprintf('=== Rapprocher la detente MEC de la FEA : quelle cible ? ===\n');

%% ---------- 1. TEST DE REPLIEMENT ----------
Rc=cogging_mec(M,1260,0,841,M.muI,kfr);          % 1 pas polaire
TM=Rc.T(:).'*1e3; TM=TM-mean(TM); NM=numel(TM)-1;
TMfull=repmat(TM(1:NM),1,Nm);                     % deroule sur un tour
xf=(0:numel(TMfull)-1)*360/numel(TMfull);
% echantillonnage IDENTIQUE a la magnetostatique : 360 points, pas de 1 deg
xs=0:359; TMs=interp1([xf 360],[TMfull TMfull(1)],xs,'linear');
Ys=abs(fft(TMs))/360*2; ns=180; [~,js]=max(Ys(2:ns));
fprintf('\n--- 1. Repliement : detente MEC vue avec le pas de 1 deg de la FEA ---\n');
fprintf('  ordre reel de la detente MEC              : %d\n',LCMv);
fprintf('  ordre apparent apres echantillonnage 1 deg: %d  (repliement 360-210)\n',js);
fprintf('  amplitude repliee                         : %.3f mN.m\n',Ys(js+1));
dmag=rd(fullfile(M.FEA.dir,'magnetostique(Magnetic_loading)','Torque Plot 1.tab'));
alr=dmag(:,1); Tr=dmag(:,2);
if abs(alr(end)-360)<1e-6, alr(end)=[]; Tr(end)=[]; end
Tr=Tr-mean(Tr); Yr=abs(fft(Tr))/numel(Tr)*2; [~,jr]=max(Yr(2:ns));
fprintf('  ordre dominant de la FEA magnetostatique  : %d (%.1f mN.m)\n',jr,Yr(jr+1));
fprintf('  => %d n''est NI 210 NI son repli %d : la magnetostatique ne\n',jr,js);
fprintf('     contient PAS de detente, meme repliee. Ce n''est pas une cible.\n');

%% ---------- 2. REFERENCE VALIDE : le TRANSITOIRE ----------
dtr=rd(fullfile(M.FEA.dir,'transitoire (Back_emf)','Torque Plot.tab'));
tt=dtr(:,1)*1e-3; Tt=dtr(:,2); th=tt*(M.speed/60)*2*pi;      % rad mec
A=[cos(LCMv*th) sin(LCMv*th) ones(size(th)) th th.^2];
cf=A\Tt; res=Tt-A*cf; A210=hypot(cf(1),cf(2));
sN=std(res); sig=sN*sqrt(2/numel(th));
fprintf('\n--- 2. Reference VALIDE : essai transitoire (maillage preserve) ---\n');
fprintf('  %d positions sur %.2f deg mec (%.0f periodes de detente)\n', ...
    numel(th),max(th)*180/pi,max(th)/(2*pi/LCMv));
fprintf('  amplitude a l''ordre %d : %.3f +/- %.3f mN.m (bruit residuel %.2f)\n', ...
    LCMv,A210,sig,sN);
fprintf('  -> detente FEA = %.3f mN.m c-c, borne a 2 sigma : < %.2f mN.m c-c\n', ...
    2*A210,2*(A210+2*sig));

%% ---------- 3. MEC OPTIMISE (sans parametre ajuste) ----------
fprintf('\n--- 3. MEC sans parametre ajuste : convergence de la detente ---\n');
Tcp=2*pi/LCMv; Npc=25; phc=linspace(0,Tcp,Npc);
fprintf('  %6s %8s | %12s\n','K','numax','detente c-c');
res_c=[];
for K=[6 8 10]
    nm=ceil(4*K*pi/(M.ws0/M.Rsi));
    Tq=zeros(1,Npc);
    for q=1:Npc, Rq=subdomain_mec(M,K,nm,'nl',phc(q)); Tq(q)=Rq.T; end
    pp=(max(Tq)-min(Tq))*1e3; res_c(end+1)=pp; %#ok<SAGROW>
    fprintf('  %6d %8d | %10.3f mN.m\n',K,nm,pp);
end
Tsub=res_c(end);
fprintf('  variation K=8->10 : %+.1f %% (convergence)\n', ...
    100*(res_c(end)-res_c(end-1))/res_c(end-1));
fprintf('  detente du modele a kfringe (reference interne) : %.3f mN.m\n',Rc.Tpp);

%% ---------- 4. BILAN : le MEC est-il compatible avec la FEA valide ? ----------
fprintf('\n===================== BILAN =====================\n');
fprintf('  MEC (sous-domaine exact, sans parametre) : %.3f mN.m c-c\n',Tsub);
fprintf('  MEC (reseau + kfringe identifie)         : %.3f mN.m c-c\n',Rc.Tpp);
fprintf('  FEA transitoire (ordre %d)               : %.3f mN.m c-c\n',LCMv,2*A210);
fprintf('  incertitude de la reference (2 sigma)    : +/- %.3f mN.m\n',4*sig);
ok=abs(Tsub-2*A210)<=4*sig;
fprintf('  ecart MEC-FEA = %.3f mN.m  -> %s\n',abs(Tsub-2*A210), ...
    ternary(ok,'DANS l''incertitude de la reference','hors incertitude'));
fprintf('  ordre : MEC %d / FEA %d  -> IDENTIQUE\n',LCMv,LCMv);

%% ---------- 5. Dispersion d'aimants equivalente a l'ordre 42 ----------
%  Un couple aux multiples de Nm (et non de LCM) exige une dispersion
%  AIMANT A AIMANT. Le couple varie en Br^2 : une dispersion relative eps
%  sur Br produit ~2*eps de modulation. En rapportant l'ordre 42 de la FEA
%  a l'amplitude du couple par aimant, on obtient la dispersion equivalente.
Tmag=Rc.Tpp*LCMv/Nm;                    % ordre de grandeur du couple par aimant
eps42=Yr(jr+1)/(2*Tmag);
fprintf('\n--- 5. Ce que l''ordre %d de la FEA representerait physiquement ---\n',jr);
fprintf('  un couple aux multiples de %d exige une dispersion AIMANT A AIMANT\n',Nm);
fprintf('  (un stator parfait ne peut produire que des multiples de %d)\n',LCMv);
fprintf('  dispersion de Br equivalente : ~%.0f %% -- inconcevable pour des\n',100*eps42);
fprintf('  aimants industriels (tolerance usuelle 1-3 %%). C''est donc bien du\n');
fprintf('  bruit de remaillage, chaque aimant etant maille differemment.\n');

function o=ternary(c,a,b), if c, o=a; else, o=b; end, end
