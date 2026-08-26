%% RUN_LOAD_MEC  -  Fonctionnement EN CHARGE entierement calcule par le MEC
%
%  Jusqu'ici la section 5 ne fournissait que le couple moyen et les pertes ;
%  toutes les grandeurs TEMPORELLES (vitesse, courants, rendement) venaient
%  de la FEA. Ce script construit le regime transitoire complet a partir du
%  seul reseau de reluctances :
%     lambda_pm(theta), L_a, M  -> reseau MEC + DtN etendu a l'aimant
%     T_detente(theta)          -> tenseur de Maxwell sur le meme reseau
%     onduleur + mecanique      -> netlist 'vf-15-14.ckt' et MotionSetup1
%  puis compare CHAQUE grandeur a l'essai ANSYS 'transitoire (en charge)'.
clear; clc; close all;
t0=tic;
M=machine_bldc(); Ns=M.Ns; Nm=M.Nm; p=M.p; Ntc=M.Ntc; fea=M.FEA.dir;
PA=[1 -2 -15 3 14]; PB=[6 -7 -5 8 4]; PC=[11 -12 -10 13 9];
if isfile('kfringe_ident.mat'), S=load('kfringe_ident.mat'); kfr=S.kbest; else, kfr=0.75; end
rd=@(f)readmatrix(f,'FileType','text','Delimiter','\t','NumHeaderLines',1);
ol=fullfile(fea,'transitoire (en charge)');
w4=@(x)x(round(numel(x)*0.75):end);

%% ---------- 1. parametres electromagnetiques : MEC seul -----------------
fprintf('=== MEC : parametres de l''entrainement ===\n');
R=cogging_mec(M,1260,0,721,M.muI,kfr,2*pi/p);
phis=R.phis; the=phis*p;
lamf=@(P) Ntc*sum(sign(P(:)).*R.PhiT(abs(P(:)),:),1);
lam=[lamf(PA); lamf(PB); lamf(PC)];
RI=inductance_mec(M,1260,5000,0); Leff=RI.Ld;
Rc=cogging_mec(M,1260,0,841,M.muI,kfr);
fprintf('  lambda_pm crete = %.4f Wb | L_d = L_a - M = %.3f mH | R = %.3f ohm\n', ...
    max(abs(lam(1,:))),Leff*1e3,M.Rph);
fprintf('  constante electrique L_d/R = %.2f ms\n',Leff/M.Rph*1e3);

%% ---------- 2. calage angulaire etabli sur l'essai A VIDE ---------------
%  Convention de repere entre le MEC et ANSYS ($theta0 = 43 deg de rotation
%  des aimants). Elle est etablie UNE FOIS sur l'essai a vide (aucun courant,
%  donc aucune ambiguite) puis appliquee TELLE QUELLE a l'essai en charge.
dFL=rd(fullfile(fea,'transitoire (Back_emf)','Output Variables Plot 4.tab'));
posFL=dFL(:,2)*p; FLaF=dFL(:,3);
cLM=2*mean(lam(1,:).*exp(-1i*the));
cLF=2*mean(FLaF.'.*exp(-1i*posFL.'*pi/180));
shL=angle(cLM/cLF); dphi=-shL;
fprintf('  calage MEC/ANSYS (essai A VIDE) : %+.2f deg elec\n',shL*180/pi);
%  controle : le meme calage doit ressortir de l'essai EN CHARGE
dFLo=rd(fullfile(ol,'Output Variables Plot 1.tab'));
dI  =rd(fullfile(ol,'BranchCurrent Plot 1.tab'));
dS  =rd(fullfile(ol,'Speed Plot 1.tab'));
tI=dI(:,1); iF=dI(:,2:4);
FLa_ol=interp1(dFLo(:,1),dFLo(:,2),tI,'linear','extrap');
omF=interp1(dS(:,1),dS(:,2),tI,'linear','extrap')*2*pi/60;
thF=cumtrapz(tI*1e-3,omF)*p;                      % angle elec ANSYS
mS=tI>0.6*tI(end);
c1F=2*mean((FLa_ol(mS)-Leff*iF(mS,1)).'.*exp(-1i*thF(mS).'));
shC=angle(cLM/c1F);
fprintf('  meme calage retrouve EN CHARGE : %+.2f deg elec (ecart %+.2f deg)\n', ...
    shC*180/pi,mod((shC-shL)*180/pi+180,360)-180);

%% ---------- 3. simulation transitoire MEC -------------------------------
opt=struct('dt',1e-6,'tend',dS(end,1)*1e-3,'J',1e-3, ...
    'Bf',6.07927101854027e-4,'Tload',4.8,'Vdc',M.Vdc,'Rph',M.Rph, ...
    'dphi',dphi,'Tcog',Rc.T,'thc',Rc.phis*p);
Dl=drive_mec(M,the,lam,Leff,opt);                 % modele LINEAIRE
if isfile('mec_map.mat')
    Q=load('mec_map.mat'); opt.map=Q.S;
    D=drive_mec(M,the,lam,Leff,opt);              % modele SATURABLE
else
    D=Dl; warning('mec_map.mat absent : modele lineaire seul');
end
fprintf('  transitoire MEC : %d pas de %.1f us sur %.2f ms (%.1f s)\n', ...
    numel(D.t),opt.dt*1e6,opt.tend*1e3,toc(t0));

%% ---------- 4. comparaison grandeur par grandeur ------------------------
tms=D.t*1e3; ws=tms>0.75*tms(end);                % fenetre de regime
dIt=rd(fullfile(ol,'BranchCurrent Plot 2.tab'));
dL =rd(fullfile(ol,'Plot 1_loss.tab'));
dE =rd(fullfile(ol,'Output Variables Plot 3.tab'));
wf=@(x)mean(w4(x));
nF=wf(dS(:,2)); PemF=wf(dL(:,7)); PcuF=wf(dL(:,4));
PfeF=wf(dL(:,3))/1e3; PpmF=wf(dL(:,2))/1e3; IbF=wf(dIt(:,2)); etaF=wf(dE(:,2));
IrmsF=sqrt(mean(w4(dI(:,2)).^2));

nM=mean(D.n(ws)); TM=mean(D.T(ws)); PcuM=mean(D.Pcu(ws));
IrmsM=sqrt(mean(D.i(1,ws).^2)); IbM=mean(D.idc(ws));
omM=mean(D.om(ws)); PemM=TM*omM;

fprintf('\n===== REGIME ETABLI : MEC vs FEA =====\n');
cmp=@(n,a,b,u)fprintf('  %-26s %11.4f %11.4f %8.1f %%   %s\n',n,a,b,100*(a-b)/b,u);
cmp('vitesse',nM,nF,'tr/min');
cmp('couple electromagnetique',TM,PemF/(nF*2*pi/60),'N.m');
cmp('puissance Pem',PemM,PemF,'W');
cmp('courant de phase rms',IrmsM,IrmsF,'A');
cmp('courant de bus moyen',IbM,IbF,'A');
cmp('pertes Joule',PcuM,PcuF,'W');

%% ---------- 5. dynamique de demarrage -----------------------------------
tr=@(t,n,nf)t(find(n>=0.95*nf,1));
fprintf('\n===== DEMARRAGE =====\n');
fprintf('  %-26s %11s %11s\n','','MEC','FEA');
fprintf('  %-26s %11.2f %11.2f  ms\n','temps de montee a 95 %', ...
    tr(tms,D.n,nM),tr(dS(:,1),dS(:,2),nF));
fprintf('  %-26s %11.2f %11.2f  A\n','crete de courant', ...
    max(abs(D.i(1,:))),max(abs(dI(:,2))));
omF2=interp1(dS(:,1),dS(:,2),dL(:,1)*1e-6,'linear','extrap')*2*pi/60;
fprintf('  %-26s %11.2f %11.2f  N.m\n','couple crete',max(D.T), ...
    max(dL(omF2>5,7)./omF2(omF2>5)));
fprintf('  %-26s %11.2f %11.2f  tr/min\n','vitesse min (recul)',min(D.n),min(dS(:,2)));
fprintf('\n  --- apport du modele SATURABLE ---\n');
fprintf('  %-26s %11s %11s %11s\n','','MEC lineaire','MEC sature','FEA');
fprintf('  %-26s %11.2f %11.2f %11.2f  A\n','crete de courant', ...
    max(abs(Dl.i(1,:))),max(abs(D.i(1,:))),max(abs(dI(:,2))));
fprintf('  %-26s %11.2f %11.2f %11.2f  ms\n','temps de montee a 95 %', ...
    tr(tms,Dl.n,mean(Dl.n(ws))),tr(tms,D.n,nM),tr(dS(:,1),dS(:,2),nF));
fprintf('  %-26s %11.0f %11.0f %11.0f  tr/min\n','vitesse finale', ...
    mean(Dl.n(ws)),nM,nF);

%% ---------- 6. figure de controle ---------------------------------------
f=figure('Color','w','Position',[30 30 1400 820],'Visible','off');
subplot(2,3,1);
plot(dS(:,1),dS(:,2),'r.','MarkerSize',7); hold on; plot(tms,D.n,'b','LineWidth',1);
grid on; xlabel('t (ms)'); ylabel('n (tr/min)'); legend('FEA','MEC','Location','se');
title(sprintf('(a) Vitesse : %.0f vs %.0f tr/min',nM,nF));
subplot(2,3,2);
plot(dI(:,1),dI(:,2),'r.','MarkerSize',6); hold on; plot(tms,D.i(1,:),'b','LineWidth',0.8);
grid on; xlabel('t (ms)'); ylabel('i_a (A)'); legend('FEA','MEC');
title('(b) Courant de phase a');
subplot(2,3,3);
plot(dI(:,1),dI(:,2),'r.-','MarkerSize',7); hold on; plot(tms,D.i(1,:),'b','LineWidth',1);
grid on; xlim([tms(end)-8 tms(end)]); xlabel('t (ms)'); ylabel('i_a (A)');
title('(c) Zoom regime etabli');
subplot(2,3,4);
plot(dIt(:,1),dIt(:,2),'r.','MarkerSize',6); hold on; plot(tms,D.idc,'b','LineWidth',0.8);
grid on; xlabel('t (ms)'); ylabel('i_{bus} (A)'); legend('FEA','MEC');
title(sprintf('(d) Courant de bus : %.3f vs %.3f A',IbM,IbF));
subplot(2,3,5);
plot(dL(:,1)*1e-6,dL(:,7),'r.','MarkerSize',6); hold on;
plot(tms,D.T.*D.om,'b','LineWidth',0.8);
grid on; xlabel('t (ms)'); ylabel('P_{em} (W)'); legend('FEA','MEC');
title(sprintf('(e) Puissance : %.0f vs %.0f W',PemM,PemF));
subplot(2,3,6);
plot(dL(:,1)*1e-6,dL(:,4),'r.','MarkerSize',6); hold on; plot(tms,D.Pcu,'b','LineWidth',0.8);
grid on; xlabel('t (ms)'); ylabel('p_{cu} (W)'); legend('FEA','MEC');
title(sprintf('(f) Pertes Joule : %.1f vs %.1f W',PcuM,PcuF));
sgtitle('Fonctionnement en charge : transitoire calcule par le MEC vs ANSYS 2D');
exportgraphics(f,'RUN_LOAD_MEC.png','Resolution',120);
fprintf('\n-> RUN_LOAD_MEC.png  (%.1f s)\n',toc(t0));
