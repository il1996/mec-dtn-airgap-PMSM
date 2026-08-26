%% ======================================================================
%  RUN_COMPARE_ALL  -  Reproduction MEC de TOUTES les etudes FEA + comparaison
%  ----------------------------------------------------------------------
%  Les 6 etudes du dossier "ANSYS resultat 750W" sont reproduites par le
%  modele MEC et confrontees une a une :
%     1. magnetostique(Magnetic_loading) : Br, Bt, spectre, lignes de flux
%     2. CIRCUIT MAGNETIQUE : FMM par aimant / par entrefer, k_l, k_r
%        (grandeurs INTERNES du circuit -- cibles ideales pour un MEC,
%         jamais exploitees jusqu'ici)
%     3. magnetostique(Armature-Field) : inductances propre et mutuelle
%     4. transitoire (Back_emf) : FEM de phase, enveloppe, flux, detente
%     5. transitoire (a vide) + (en charge) : pertes, couple, courants
%     6. balayage : couple / puissance / courant / rendement / FEM vs vitesse
%
%  Toutes les grandeurs materiau proviennent du PROJET ANSYS lui-meme
%  (BLDC.aedt) : courbe M350-50A reelle, aimant N42UH (mu_r=1.039,
%  Hc=955204 A/m, sigma=555556 S/m).
%% ======================================================================
clear; clc; close all;
M=machine_bldc(); mu0=4*pi*1e-7;
Rs=M.Rsi; L=M.ls; Ns=M.Ns; Nm=M.Nm; p=M.p; Ntc=M.Ntc; fea=M.FEA.dir;
PA=[1 -2 -15 3 14]; PB=[6 -7 -5 8 4]; PC=[11 -12 -10 13 9];
if isfile('kfringe_ident.mat'), S=load('kfringe_ident.mat'); kfr=S.kbest; else, kfr=0.75; end
rd=@(f)readmatrix(f,'FileType','text','Delimiter','\t','NumHeaderLines',1);
amp=@(B,th,k) 2*abs(mean(B(:).*exp(-1i*k*th(:))));
SC=struct('n',{},'m',{},'f',{},'u',{});
add=@(S,n,m,f,u)[S,struct('n',n,'m',m,'f',f,'u',u)];
fprintf('=== MEC BLDC 15/14 750 W : reproduction de TOUTES les etudes FEA ===\n');
fprintf('    materiaux = ceux du projet ANSYS (M350 reelle, N42UH mu_r=%.4f)\n\n',M.mu_r);

%% ---------- resolution MEC (1 periode electrique) ----------
Nsurf=1260; Np=721;
R=cogging_mec(M,Nsurf,0,Np,M.muI,kfr,2*pi/p);
phis=R.phis; om=M.speed*2*pi/60;

%% ================= 1. CHAMP D'ENTREFER =================
d4=rd(fullfile(fea,'magnetostique(Magnetic_loading)','Calculator Expressions Plot 4.tab'));
d2=rd(fullfile(fea,'magnetostique(Magnetic_loading)','Calculator Expressions Plot 2.tab'));
dA=rd(fullfile(fea,'magnetostique(Magnetic_loading)','Plot 1.tab'));
angF=d4(:,2)*pi/180; BrF=d4(:,3); BtF=d2(:,2);
AF_fea=dA(:,2);                                   % lignes de flux (potentiel vecteur)
thu=linspace(0,2*pi,3601); thu(end)=[];
BrFu=interp1([angF;2*pi],[BrF;BrF(1)],thu,'linear','extrap');
BtFu=interp1([angF;2*pi],[BtF;BtF(1)],thu,'linear','extrap');
thm=R.thq(1:end-1); Brm=R.Br(1:end-1); Btm=R.Bt(1:end-1);
% recalage de phase sur le fondamental (origine angulaire FEA arbitraire)
sh=angle( (2*mean(Brm.*exp(-1i*p*thm))) / (2*mean(BrFu.*exp(-1i*p*thu))) )/p;
thmr=mod(thm-sh,2*pi); [thmr,is]=sort(thmr); Brmr=Brm(is); Btmr=Btm(is);
ords=[p 8 22 23 21 35 37];
a_m=arrayfun(@(k)amp(Brm,thm,k),ords); a_f=arrayfun(@(k)amp(BrFu,thu,k),ords);
SC=add(SC,'Bg1 fondamental',a_m(1),a_f(1),'T');
SC=add(SC,'B entrefer moyen |Br|',mean(abs(Brm)),mean(abs(BrFu)),'T');
SC=add(SC,'B entrefer crete',max(Brm),max(BrFu),'T');
SC=add(SC,'Bt tangentiel RMS',sqrt(mean(Btm.^2)),sqrt(mean(BtFu.^2)),'T');
% potentiel vecteur A (lignes de flux) : A(th) = Rs*int(Br dth)
Amec=cumtrapz(thmr,Brmr)*M.rmid; Amec=Amec-mean(Amec);
Afea=AF_fea-mean(AF_fea);
SC=add(SC,'Flux lines A crete',max(Amec),max(Afea),'Wb/m');

%% ================= 2. CIRCUIT MAGNETIQUE (FMM, k_l, k_r) =================
%  FEA : mmf_m1..14 (chute dans chaque aimant), mmf_g1..14 (chaque entrefer),
%  k_l (facteur de fuite), k_r (facteur de reluctance).
dm=rd(fullfile(fea,'magnetostique(Magnetic_loading)','Calculator Expressions Table 1.tab'));
mmf_m_F=dm(1,2:15); mmf_g_F=dm(1,16:29);
kl_F=rd(fullfile(fea,'magnetostique(Magnetic_loading)','Output Variables Table 1.tab')); kl_F=kl_F(1,2);
kr_F=rd(fullfile(fea,'magnetostique(Magnetic_loading)','Output Variables Table 2.tab')); kr_F=kr_F(1,2);
% MEC : potentiel a l'interface aimant/entrefer (culasse rotor = 0)
AG=R.AG; nu=AG.nu; Us=R.U(1:Nsurf,1);
Usc=AG.Wc*Us; Uss=AG.Ws*Us;
pic=AG.alphau.*Usc + AG.betasrc;           % phi_i (cos), rotor a phi=0
pis=AG.alphau.*Uss;
thq=linspace(0,2*pi,2001); thq(end)=[];
phi_i =(pic.'*cos(nu*thq) + pis.'*sin(nu*thq));      % potentiel en Rro
phi_b =(Usc.'*cos(nu*thq) + Uss.'*sin(nu*thq));      % potentiel au bore
%  IMPORTANT : ANSYS sonde une LIGNE RADIALE au centre de chaque aimant
%  (mmf_mk = int H.dl le long de cette ligne). Moyenner sur l'arc d'aimant
%  (24 deg) lisserait la denture et ecraserait artificiellement la
%  dispersion -> on evalue donc AU CENTRE de chaque aimant.
mmf_m_M=zeros(1,Nm); mmf_g_M=zeros(1,Nm);
for k=1:Nm
    c=(k-1)*2*pi/Nm;
    [~,j]=min(abs(angle(exp(1i*(thq-c)))));
    mmf_m_M(k)=abs(phi_i(j));                       % chute dans l'aimant
    mmf_g_M(k)=abs(phi_i(j)-phi_b(j));              % chute dans l'entrefer
end
% recalage : la position rotorique de l'essai FEA est arbitraire -> on
% compare la STRUCTURE (moyenne + dispersion), invariante par rotation.
SC=add(SC,'FMM aimant (moyenne)',mean(mmf_m_M),mean(mmf_m_F),'A');
SC=add(SC,'FMM entrefer (moyenne)',mean(mmf_g_M),mean(mmf_g_F(2:end)),'A');
SC=add(SC,'FMM aimant : dispersion',std(mmf_m_M)/mean(mmf_m_M)*100, ...
        std(mmf_m_F)/mean(mmf_m_F)*100,'%');
SC=add(SC,'Facteur de reluctance k_r',1+ (mean(mmf_m_M)-mean(mmf_g_M))/mean(mmf_g_M),kr_F,'-');

%% ================= 3. INDUCTANCES =================
dL=rd(fullfile(fea,'magnetostique(Armature-Field)','Output Variables Table 1.tab'));
LaF=dL(1,2); MF=(dL(1,3)+dL(1,4))/2; LdF=LaF-MF;
RI=inductance_mec(M,1260,5000,0);
SC=add(SC,'Inductance propre La',RI.La*1e3,LaF*1e3,'mH');
SC=add(SC,'Mutuelle M',RI.M*1e3,MF*1e3,'mH');
SC=add(SC,'Inductance synchrone Ld',RI.Ld*1e3,LdF*1e3,'mH');

%% ================= 4. BACK-EMF ET DETENTE =================
lamf=@(P) Ntc*sum(sign(P(:)).*R.PhiT(abs(P(:)),:),1);
lamA=lamf(PA); lamB=lamf(PB); lamC=lamf(PC);
eA=-gradient(lamA,phis)*om; eB=-gradient(lamB,phis)*om; eC=-gradient(lamC,phis)*om;
env=max([eA;eB;eC],[],1)-min([eA;eB;eC],[],1);
dE=rd(fullfile(fea,'transitoire (Back_emf)','Output Variables Plot.tab'));
dV=rd(fullfile(fea,'transitoire (Back_emf)','NodeVoltage Plot.tab'));
dFL=rd(fullfile(fea,'transitoire (Back_emf)','Output Variables Plot 4.tab'));
posE=dE(:,2)*p; envF=dE(:,3); posV=dV(:,2)*p; EaF=dV(:,3);
posFL=dFL(:,2)*p; FLaF=dFL(:,3);
SC=add(SC,'FEM de phase crete',max(abs(eA)),max(abs(EaF)),'V');
SC=add(SC,'Enveloppe six-step crete',max(env),max(envF),'V');
SC=add(SC,'Flux totalise crete',max(abs(lamA)),max(abs(FLaF)),'Wb');
Rc=cogging_mec(M,1260,0,841,M.muI,kfr);
dC=rd(fullfile(fea,'transitoire (Back_emf)','Torque Plot.tab'));
cogF_pp=max(dC(:,2))-min(dC(:,2));

%% ================= 5. CHARGE ET PERTES =================
Iflat=M.Iph_rms_load*sqrt(3/2);
ss=@(e)sixstep(e,Iflat);
T_ideal=mean((eA.*ss(eA)+eB.*ss(eB)+eC.*ss(eC))/om);   % creneau idealise
dOL=rd(fullfile(fea,'transitoire (en charge)','Plot 1_loss.tab'));
dOLs=rd(fullfile(fea,'transitoire (en charge)','Speed Plot 1.tab'));
dOLe=rd(fullfile(fea,'transitoire (en charge)','Output Variables Plot 3.tab'));
dOLi=rd(fullfile(fea,'transitoire (en charge)','BranchCurrent Plot 1.tab'));
w4=@(x)x(round(numel(x)*0.75):end);
nOL=mean(w4(dOLs(:,2))); PemOL=mean(w4(dOL(:,7))); TF_load=PemOL/(nOL*2*pi/60);
% --- couple RIGOUREUX : courant REEL de la FEA + alignement corrige -------
%  En charge, FL_a = lambda_aimant + Ld*i_a : caler la phase sur FL_a BRUT
%  biaise l'alignement de ~20 deg elec (Ld*i_a = 0.13 Wb pour 0.26 Wb d'aimant)
%  et le couple varie de ~1 %/deg. On retranche donc Ld*i_a avant recalage.
dFLo=rd(fullfile(fea,'transitoire (en charge)','Output Variables Plot 1.tab'));
tI=dOLi(:,1); iabc=dOLi(:,2:4);
tF=dFLo(:,1); FLa_ol=interp1(tF,dFLo(:,2),tI,'linear','extrap');
lam_pm_fea=FLa_ol - LdF*iabc(:,1);            % flux d'aimant seul (FEA)
omOL=nOL*2*pi/60; theOL=omOL*p*tI*1e-3;       % angle electrique de l'essai
% fondamental de lambda_pm(FEA) et de lambda_a(MEC) -> dephasage
c1F=2*mean(lam_pm_fea.'.*exp(-1i*theOL.'));
c1M=2*mean(lamA.*exp(-1i*phis*p));
dps=angle(c1M/c1F);
% back-EMF MEC evaluee aux angles de l'essai, alignee
eAi=interp1(phis*p,eA,mod(theOL+dps,2*pi),'linear','extrap');
eBi=interp1(phis*p,eB,mod(theOL+dps,2*pi),'linear','extrap');
eCi=interp1(phis*p,eC,mod(theOL+dps,2*pi),'linear','extrap');
kOL=omOL/om;                                   % FEM proportionnelle a la vitesse
T_load=mean(kOL*(eAi.*iabc(:,1)+eBi.*iabc(:,2)+eCi.*iabc(:,3)))/omOL;
SC=add(SC,'Couple en charge (courant reel)',T_load,TF_load,'N.m');
SC=add(SC,'Couple (creneau idealise)',T_ideal,TF_load,'N.m');
SC=add(SC,'Courant de phase RMS',M.Iph_rms_load,sqrt(mean(w4(dOLi(:,2)).^2)),'A');
% pertes a vide
f=M.FEA.n_nl*Nm/120;
Vt1=M.wst1*(M.hs0+M.hs1+M.hs2)*L*M.Ki; Vy1=(2*pi*(M.Rso-M.wsy/2)/Ns)*M.wsy*L*M.Ki;
tau_yi=2*pi*(M.Rsi+M.hs0+M.hs1+M.hs2)/Ns;
Bt_=R.PhiT/(M.wst1*L*M.Ki); By_=R.PhiY/(M.wsy*L*M.Ki); Byr_=R.PhiT/(tau_yi*L*M.Ki);
tt=linspace(0,1/f,Np);
bl=@(B,V)deal(sum(M.Kh*f*((max(B,[],2)-min(B,[],2))/2).^2)*V, ...
              sum((M.KeFe/(2*pi^2))*mean(gradient(B,tt(2)-tt(1)).^2,2))*V);
[Ph_t,Pe_t]=bl(Bt_,Vt1); [Ph_y,Pe_y]=bl(By_,Vy1); [Ph_r,Pe_r]=bl(Byr_,Vy1);
Pfe=Ph_t+Pe_t+Ph_y+Pe_y+Ph_r+Pe_r;
[Ppm,~]=pm_loss(M,840,91,kfr);
dNL=rd(fullfile(fea,'transitoire (a vide)','Plot 1_loss.tab'));
SC=add(SC,'Pertes fer a vide',Pfe,mean(w4(dNL(:,3)))/1000,'W');
SC=add(SC,'Pertes aimant a vide',Ppm,mean(w4(dNL(:,2)))/1000,'W');

%% ================= 6. BALAYAGE EN VITESSE =================
nS =rd(fullfile(fea,'balayage','Torque Plot 2.tab'));
TS =nS(:,2); nS=nS(:,1);
PS =rd(fullfile(fea,'balayage','Output Variables Plot 3.tab')); PS=PS(:,2)*1000;
ES =rd(fullfile(fea,'balayage','Output Variables Plot 2.tab')); ES=ES(:,2);
FS =rd(fullfile(fea,'balayage','Output Variables Plot 6.tab')); FS=FS(:,2);
IS =rd(fullfile(fea,'balayage','BranchCurrent Plot 3.tab')); IS=IS(:,2);
etaS=rd(fullfile(fea,'balayage','Output Variables Plot 1.tab')); etaS=etaS(:,2);
[nS,o]=sort(nS); TS=TS(o); PS=PS(o); ES=ES(o); FS=FS(o); IS=IS(o); etaS=etaS(o);
% MEC : FEM proportionnelle a la vitesse ; flux totalise constant
kE=max(env)/M.speed;                       % V par tr/min (enveloppe)
E_mec=kE*nS;
SC=add(SC,'FEM enveloppe @1500 (balayage)',kE*1500,interp1(nS,ES,1500),'V');
SC=add(SC,'Constante de FEM k_E',kE*1000,interp1(nS,ES,1500)/1500*1000,'mV/tr/min');
%  NB : rms(FL_a) du balayage vaut 0.268 Wb > la CRETE du flux d'aimant
%  (0.258 Wb) : c'est le flux EN CHARGE (lambda_aimant + L*i), pas le flux
%  d'aimant seul -> non comparable au flux a vide. On ne le confronte donc pas.
FL_mec=[];

%% ====================== SCORECARD ======================
fprintf('  %-34s %11s %11s %9s\n','Grandeur','MEC','FEA','Ecart');
for k=1:numel(SC)
    e=(SC(k).m-SC(k).f)/SC(k).f*100;
    fprintf('  %-34s %11.4f %11.4f %+8.1f%%  %s\n',SC(k).n,SC(k).m,SC(k).f,e,SC(k).u);
end
fprintf('  %-34s %11.3f %11s %9s  mN.m\n','Couple de detente c-c',Rc.Tpp,'< 1.4','(borne)');
fprintf('  %-34s %11d %11d %9s\n','Ordre de detente',Rc.order,lcm(Ns,Nm),'exact');

%% ====================== FIGURES ======================
% --- FIG 1 : champ d'entrefer ---
f1=figure('Color','w','Position',[40 40 1200 780],'Visible','off');
subplot(2,2,1); plot(thu*180/pi,BrFu,'r','LineWidth',1.3); hold on;
plot(thmr*180/pi,Brmr,'b--','LineWidth',1.2); grid on; xlim([0 360]);
xlabel('\theta mec (deg)'); ylabel('B_r (T)'); legend('FEA','MEC','Location','best');
title('(a) Induction radiale d''entrefer');
subplot(2,2,2); plot(thu*180/pi,BrFu,'r','LineWidth',1.5); hold on;
plot(thmr*180/pi,Brmr,'b--','LineWidth',1.4); grid on; xlim([0 72]);
xlabel('\theta mec (deg)'); ylabel('B_r (T)'); title('(b) zoom : 3 pas d''encoche');
subplot(2,2,3); plot(thu*180/pi,BtFu,'r','LineWidth',1.3); hold on;
plot(thmr*180/pi,Btmr,'b--','LineWidth',1.2); grid on; xlim([0 72]);
xlabel('\theta mec (deg)'); ylabel('B_t (T)'); legend('FEA','MEC'); title('(c) Induction tangentielle');
subplot(2,2,4); kk=1:40;
aFk=arrayfun(@(k)amp(BrFu,thu,k),kk); aMk=arrayfun(@(k)amp(Brm,thm,k),kk);
bar(kk-0.2,aFk,0.4,'r'); hold on; bar(kk+0.2,aMk,0.4,'b'); set(gca,'YScale','log');
ylim([1e-4 2]); grid on; xlabel('ordre spatial'); ylabel('|B_{r,n}| (T)');
legend('FEA','MEC'); title('(d) Spectre spatial');
sgtitle('Etude 1 - Magnetic loading : champ d''entrefer MEC vs FEA');
exportgraphics(f1,'CMP_1_champ.png','Resolution',120);

% --- FIG 2 : circuit magnetique ---
f2=figure('Color','w','Position',[40 40 1200 420],'Visible','off');
subplot(1,3,1); bar([mmf_m_F(:) mmf_m_M(:)]); grid on;
xlabel('aimant #'); ylabel('FMM (A)'); legend('FEA','MEC','Location','south');
title('(a) FMM par aimant');
subplot(1,3,2); bar([mmf_g_F(:) mmf_g_M(:)]); grid on;
xlabel('entrefer #'); ylabel('FMM (A)'); legend('FEA','MEC','Location','south');
title('(b) FMM par entrefer');
subplot(1,3,3);
bar([kl_F kr_F; NaN 1+(mean(mmf_m_M)-mean(mmf_g_M))/mean(mmf_g_M)].'); grid on;
set(gca,'XTickLabel',{'k_l','k_r'}); ylabel('-'); legend('FEA','MEC');
title('(c) Facteurs de fuite / reluctance');
sgtitle('Etude 2 - Circuit magnetique interne (FMM par element) MEC vs FEA');
exportgraphics(f2,'CMP_2_circuit.png','Resolution',120);

% --- FIG 3 : back-EMF + detente ---
f3=figure('Color','w','Position',[40 40 1200 780],'Visible','off');
ae=phis*p*180/pi;
shE=angle((2*mean(eA.*exp(-1i*ae*pi/180)))/(2*mean(EaF.'.*exp(-1i*posV.'*pi/180))))*180/pi;
subplot(2,2,1); plot(mod(posV,360),EaF,'r.','MarkerSize',7); hold on;
plot(mod(ae-shE,360),eA,'b.','MarkerSize',4); grid on; xlim([0 360]);
xlabel('angle elec (deg)'); ylabel('e_a (V)'); legend('FEA','MEC'); title('(a) FEM de phase');
subplot(2,2,2); plot(mod(posE,360),envF,'r.','MarkerSize',7); hold on;
plot(mod(ae-shE,360),env,'b.','MarkerSize',4); grid on; xlim([0 360]);
xlabel('angle elec (deg)'); ylabel('e_{LL} (V)'); legend('FEA','MEC'); title('(b) Enveloppe six-step');
subplot(2,2,3); plot(mod(posFL,360),FLaF,'r.','MarkerSize',7); hold on;
plot(mod(ae-shE,360),lamA,'b.','MarkerSize',4); grid on; xlim([0 360]);
xlabel('angle elec (deg)'); ylabel('\lambda_a (Wb)'); legend('FEA','MEC'); title('(c) Flux totalise');
subplot(2,2,4); plot(Rc.phis*180/pi,Rc.T*1e3,'b','LineWidth',1.3); hold on;
plot(dC(:,1)/max(dC(:,1))*max(Rc.phis)*180/pi,dC(:,2)-mean(dC(:,2)),'r','LineWidth',0.8);
grid on; xlabel('position rotor (deg mec)'); ylabel('T (mN.m)');
legend('MEC','FEA (bruit)'); title(sprintf('(d) Detente : MEC %.2f mN.m, ordre %d',Rc.Tpp,Rc.order));
sgtitle('Etude 3-4 - Back-EMF, flux totalise et couple de detente');
exportgraphics(f3,'CMP_3_bemf.png','Resolution',120);

% --- FIG 4 : charge et pertes ---
f4=figure('Color','w','Position',[40 40 1200 420],'Visible','off');
subplot(1,3,1);
plot(dOLs(:,1),dOLs(:,2),'r','LineWidth',1.1); grid on;
xlabel('t (ms)'); ylabel('n (tr/min)'); title(sprintf('(a) Vitesse FEA (%.0f tr/min)',nOL));
subplot(1,3,2);
tI=dOLi(:,1); plot(tI,dOLi(:,2:4),'LineWidth',0.9); grid on;
xlabel('t (ms)'); ylabel('i (A)'); title('(b) Courants de phase FEA (six-step)');
subplot(1,3,3);
bar([Pfe mean(w4(dNL(:,3)))/1000; Ppm*10 mean(w4(dNL(:,2)))/1000*10; ...
     T_load TF_load].'); grid on;
ylabel('W  /  10xW  /  N.m');
set(gca,'XTickLabel',{'MEC','FEA'}); legend('P_{fe}(W)','P_{aim}(W)','T(N.m)');
title('(c) Pertes a vide et couple');
sgtitle('Etude 5 - Fonctionnement en charge et pertes');
exportgraphics(f4,'CMP_4_charge.png','Resolution',120);

% --- FIG 5 : balayage vitesse ---
f5=figure('Color','w','Position',[40 40 1200 780],'Visible','off');
subplot(2,2,1); plot(nS,TS,'r-o','LineWidth',1.2,'MarkerSize',3); grid on;
xlabel('n (tr/min)'); ylabel('T (N.m)'); title('(a) Couple - vitesse (FEA)');
subplot(2,2,2); plot(nS,PS,'r-o','LineWidth',1.2,'MarkerSize',3); grid on;
xlabel('n (tr/min)'); ylabel('P_{em} (W)'); title('(b) Puissance - vitesse (FEA)');
subplot(2,2,3); plot(nS,ES,'r-o','LineWidth',1.2,'MarkerSize',3); hold on;
plot(nS,E_mec,'b--','LineWidth',1.4); grid on;
xlabel('n (tr/min)'); ylabel('BEMF (V)'); legend('FEA','MEC (k_E\cdotn)','Location','se');
title('(c) FEM - vitesse : MEC vs FEA');
subplot(2,2,4); plot(nS,IS,'r-o','LineWidth',1.2,'MarkerSize',3); grid on;
xlabel('n (tr/min)'); ylabel('rms(i_a) (A)');
title('(d) Courant de phase - vitesse (FEA)');
sgtitle('Etude 6 - Balayage en vitesse');
exportgraphics(f5,'CMP_5_balayage.png','Resolution',120);

% --- FIG 6 : scorecard ---
f6=figure('Color','w','Position',[60 60 900 520],'Visible','off');
er=arrayfun(@(s)(s.m-s.f)/s.f*100,SC); nmv={SC.n};
b=barh(er,'FaceColor','flat'); grid on; set(gca,'YTick',1:numel(er),'YTickLabel',nmv);
xlabel('ecart MEC - FEA (%)'); title('Scorecard MEC vs FEA (ANSYS Maxwell 2D)');
for k=1:numel(er)
    if abs(er(k))<5, b.CData(k,:)=[0.15 0.65 0.15];
    elseif abs(er(k))<15, b.CData(k,:)=[0.9 0.7 0.1];
    else, b.CData(k,:)=[0.85 0.33 0.1]; end
end
xline(0,'k'); exportgraphics(f6,'CMP_6_scorecard.png','Resolution',120);

fprintf('\n>> Figures : CMP_1_champ / CMP_2_circuit / CMP_3_bemf / CMP_4_charge /\n');
fprintf('             CMP_5_balayage / CMP_6_scorecard  (.png)\n');

%% ---------- utilitaires ----------
function i=sixstep(e,Iflat)
    n=numel(e); k=(0:n-1)/n*2*pi;
    psi=atan2(2*mean(e.*sin(k)),2*mean(e.*cos(k)));
    ke=mod(k-psi,2*pi); i=zeros(size(e));
    i(ke<pi/3|ke>5*pi/3)=Iflat; i(ke>2*pi/3&ke<4*pi/3)=-Iflat;
end
function y=Brmi(thq,th,B), y=interp1(th,B,thq,'linear','extrap'); end
function g=gradPhi(~,~), g=0; end
