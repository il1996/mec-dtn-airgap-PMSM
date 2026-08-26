%% ======================================================================
%  MEC_BLDC_MASTER  -  Circuit magnetique equivalent (MEC) du BLDC 15/14 750W
%  ----------------------------------------------------------------------
%  Applique la methode MEC (reseau de reluctances, Silva 2023) developpee et
%  validee sur le moteur asynchrone 18,5 kW (MEC_IM) a une 2e machine :
%  BLDC SPM 15 encoches / 14 poles, 750 W, 1500 tr/min, 500 V.
%
%  APPORT ORIGINAL : operateur d'entrefer DtN ETENDU A LA COURONNE D'AIMANT
%  (sous-domaines analytiques aimant + entrefer, magnetisation radiale =
%  source), COUPLE au reseau de reluctances statorique par les potentiels de
%  surface -> la DENTURE, le champ TANGENTIEL et le couple de DETENTE
%  emergent, au lieu d'etre absents d'un modele lisse.
%
%  Modele d'ouverture d'encoche IDENTIFIE sur les bandes laterales de denture
%  du champ FEA (kfringe = 0.75, RUN_SLOT_IDENT) - jamais sur les grandeurs
%  a predire. Pertes fer LOCALES (C2) sur les flux de dent/culasse du reseau.
%
%  Modules : machine_bldc, magnet_subdomain, airgap_magnet, cogging_mec.
%  Etudes  : RUN_COG_FEA (critique de la reference), RUN_SLOT_IDENT (A1),
%            RUN_COG_CONV (convergence), RUN_LOSS (C2), RUN_FIELD/BEMF/PERF.
%% ======================================================================
clear; clc; close all;
M=machine_bldc(); mu0=4*pi*1e-7;
Rs=M.Rsi; L=M.ls; Ns=M.Ns; p=M.p; Ntc=M.Ntc; taus=2*pi/Ns;
PA=[1 -2 -15 3 14]; PB=[6 -7 -5 8 4]; PC=[11 -12 -10 13 9];
if isfile('kfringe_ident.mat'), S=load('kfringe_ident.mat'); kfr=S.kbest; else, kfr=0.75; end
fprintf('=== MEC BLDC 15/14 750W : DtN etendu COUPLE au reseau stator ===\n');
fprintf('    (modele d''encoche identifie : kfringe = %.2f)\n\n',kfr);

%% ---- resolution couplee sur 1 periode ELECTRIQUE ----------------------
Nsurf=1260; Np=721;
R=cogging_mec(M,Nsurf,0,Np,M.muI,kfr,2*pi/p);
phis=R.phis;  om=M.speed*2*pi/60;

%% ---- 1. champ a vide (dente) ------------------------------------------
fea=M.FEA.dir;
d4=readmatrix(fullfile(fea,'magnetostique(Magnetic_loading)','Calculator Expressions Plot 4.tab'),'FileType','text','NumHeaderLines',1);
d2=readmatrix(fullfile(fea,'magnetostique(Magnetic_loading)','Calculator Expressions Plot 2.tab'),'FileType','text','NumHeaderLines',1);
angF=d4(:,2)*pi/180; BrF=d4(:,3); BtF=d2(:,2);
thu=linspace(0,2*pi,3601); thu(end)=[];
BrFu=interp1([angF;2*pi],[BrF;BrF(1)],thu,'linear','extrap');
amp=@(B,th,k) 2*abs(mean(B.*exp(-1i*k*th)));
thm=R.thq(1:end-1); Brm=R.Br(1:end-1);
a_mec=arrayfun(@(k) amp(Brm,thm,k),[p 8 22 23]);
a_fea=arrayfun(@(k) amp(BrFu,thu,k),[p 8 22 23]);

%% ---- 2. back-EMF a partir des FLUX DE DENT du reseau -------------------
lam=@(P) Ntc*sum(sign(P(:)).*R.PhiT(abs(P(:)),:),1);
lamA=lam(PA); lamB=lam(PB); lamC=lam(PC);
eA=-gradient(lamA,phis)*om; eB=-gradient(lamB,phis)*om; eC=-gradient(lamC,phis)*om;
env=max([eA;eB;eC],[],1)-min([eA;eB;eC],[],1);
Eph_pk=max(abs(eA)); env_pk=max(env); FLa_pk=max(abs(lamA));

%% ---- 3. couple en charge (six-step) ------------------------------------
Iflat=M.Iph_rms_load*sqrt(3/2);
ss=@(e) sixstep(e,Iflat);
T_load=mean((eA.*ss(eA)+eB.*ss(eB)+eC.*ss(eC))/om);

%% ---- 4. inductances : MEC couple + FMM distribuee (amelioration 3) -----
%  Trapeze d'encoche integre exactement, ampere-tours comptes PAR ENCOCHE,
%  part d'entrefer issue du reseau couple. kfringe = 0 ici : la frange sert
%  a l'entree RADIALE du flux a vide, pas au pont TANGENTIEL dent-a-dent
%  (deja dans lambda_tip) - cf. en-tete de inductance_mec.
RI=inductance_mec(M,1260,5000,0);
La=RI.La; Mmut=RI.M; Ld=RI.Ld;

%% ---- 5. pertes fer LOCALES sur le champ dente (C2) ---------------------
%  Champ ROTATIONNEL par decomposition a DEUX COMPOSANTES (modele d'ANSYS
%  lui-meme) : la culasse porte un flux tangentiel ET le flux radial injecte
%  par la dent ; le locus de B y est quasi circulaire (rapport 0.82).
f=M.FEA.n_nl*M.Nm/120;
Vt1=M.wst1*(M.hs0+M.hs1+M.hs2)*L*M.Ki; Vy1=(2*pi*(M.Rso-M.wsy/2)/Ns)*M.wsy*L*M.Ki;
tau_yi=2*pi*(M.Rsi+M.hs0+M.hs1+M.hs2)/Ns;
Bt_=R.PhiT/(M.wst1*L*M.Ki);                  % dent : radial
By_=R.PhiY/(M.wsy*L*M.Ki);                   % culasse : tangentiel (branche)
Byr_=R.PhiT/(tau_yi*L*M.Ki);                 % culasse : radial
tt=linspace(0,1/f,Np);
bl=@(B,V) deal(sum(M.Kh*f*((max(B,[],2)-min(B,[],2))/2).^2)*V, ...
                sum((M.KeFe/(2*pi^2))*mean(gradient(B,tt(2)-tt(1)).^2,2))*V);
[Ph_t,Pe_t]=bl(Bt_,Vt1); [Ph_y,Pe_y]=bl(By_,Vy1); [Ph_r,Pe_r]=bl(Byr_,Vy1);
Ph_y=Ph_y+Ph_r; Pe_y=Pe_y+Pe_r;
Pfe=Ph_t+Pe_t+Ph_y+Pe_y;

%% ---- 6. pertes AIMANT par courants de Foucault (C3) -------------------
[Ppm,Dpm]=pm_loss(M,840,91,kfr);

%% ---- 7. couple de detente (balayage 1 pas polaire) --------------------
Rc=cogging_mec(M,1260,0,841,M.muI,kfr);

%% ====================== SCORECARD =====================================
fprintf('  %-34s %10s %10s %9s\n','Grandeur','MEC','FEA','Ecart');
row=@(n,a,ff,u) fprintf('  %-34s %10.3f %10.3f %+8.1f%%  %s\n',n,a,ff,(a-ff)/ff*100,u);
fprintf(' [Champ a vide - magnetostatique]\n');
row('Bg1 fondamental (T)',a_mec(1),a_fea(1),'');
row('B entrefer moyen |Br| (T)',mean(abs(R.Br)),mean(abs(BrFu)),'');
row('B entrefer crete (T)',max(R.Br),max(BrFu),'');
row('Bt tangentiel RMS (T)',R.Btrms,sqrt(mean(BtF.^2)),'');
fprintf(' [Denture - bandes laterales (identification A1)]\n');
row('harmonique 8 = |Ns-p| (T)',a_mec(2),a_fea(2),'');
row('harmonique 22 = Ns+p (T)',a_mec(3),a_fea(3),'');
row('harmonique 23 = |2Ns-p| (T)',a_mec(4),a_fea(4),'');
fprintf(' [Back-EMF - 1500 tr/min, flux de dent du reseau]\n');
row('FEM phase crete (V)',Eph_pk,M.FEA.emf_ph,'');
row('Enveloppe six-step crete (V)',env_pk,M.FEA.emf_env,'');
row('Flux totalise crete (Wb)',FLa_pk,M.FEA.FLa_pk,'');
fprintf(' [Charge]\n');
row('Couple electromagnetique (N.m)',T_load,M.FEA.T_load,'(puissance)');
fprintf(' [Inductances - MEC couple + FMM distribuee]\n');
row('Inductance propre La (mH)',La*1e3,M.FEA.La*1e3,'');
row('Mutuelle M (mH)',Mmut*1e3,-2.132,'(cible independante)');
row('Inductance synchrone Ld (mH)',Ld*1e3,M.FEA.Ld*1e3,'');
fprintf(' [Pertes a vide - %g tr/min, LOCALES sur champ dente]\n',M.FEA.n_nl);
row('Pertes fer (W)',Pfe,M.FEA.Pfe_nl,'(vectoriel dents+culasse)');
row('Pertes aimant (W)',Ppm,M.FEA.Psolid_nl,'(Foucault, denture 422 Hz)');
fprintf(' [Couple de detente - ordre %d = LCM]\n',lcm(Ns,M.Nm));
fprintf('  %-34s %10.3f %10s %9s  %s\n','Detente c-c (mN.m)',Rc.Tpp,'< 1.4','(borne)', ...
        '(FEA : bruit >> signal)');
fprintf('  %-34s %10d %10d %9s\n','Ordre de detente (/tour)',Rc.order,lcm(Ns,M.Nm),'exact');

fprintf('\n  Pertes fer : %.2f W = hyst %.2f + Foucault %.2f (dents %.0f %%)\n',...
        Pfe,Ph_t+Ph_y,Pe_t+Pe_y,(Ph_t+Pe_t)/Pfe*100);
fprintf('  B dent crete %.3f T | B culasse crete %.3f T\n',...
        max((max(Bt_,[],2)-min(Bt_,[],2))/2),max((max(By_,[],2)-min(By_,[],2))/2));

%% ====================== FIGURES =======================================
figure('Name','MEC BLDC couple vs FEA','Color','w','Position',[60 60 1150 700]);
subplot(2,3,1);
plot(thu*180/pi,BrFu,'r--','LineWidth',1.0); hold on;
plot(R.thq*180/pi,R.Br,'b','LineWidth',1.3); grid on; xlim([0 75]);
xlabel('\theta mec (deg)'); ylabel('B_r (T)'); legend('FEA','MEC','Location','south');
title('Champ a vide (dente)');
subplot(2,3,2);
plot(angF*180/pi,BtF,'r--','LineWidth',1.0); hold on;
plot(R.thq*180/pi,R.Bt,'b','LineWidth',1.3); grid on; xlim([0 75]);
xlabel('\theta mec (deg)'); ylabel('B_t (T)'); title('Champ tangentiel');
subplot(2,3,3);
plot(phis*p*180/pi,eA,phis*p*180/pi,eB,phis*p*180/pi,eC,'LineWidth',1.1); hold on;
plot(phis*p*180/pi,env,'k:','LineWidth',1.2); grid on;
xlabel('angle elec (deg)'); ylabel('FEM (V)');
title(sprintf('Back-EMF (phase %.0f V)',Eph_pk));
subplot(2,3,4);
plot(phis*p*180/pi,Bt_(1,:),'b','LineWidth',1.3); hold on;
plot(phis*p*180/pi,By_(1,:),'r','LineWidth',1.3); grid on;
xlabel('angle elec (deg)'); ylabel('B (T)'); legend('dent','culasse');
title(sprintf('Inductions locales -> P_{fe} = %.1f W',Pfe));
subplot(2,3,5);
plot(Rc.phis*180/pi,Rc.T*1e3,'b','LineWidth',1.3); grid on;
xlabel('position rotor (deg mec)'); ylabel('T (mN.m)');
title(sprintf('Detente : %.2f mN.m c-c, ordre %d',Rc.Tpp,Rc.order));
subplot(2,3,6);
nmv={'Bg1','Bmoy','Bcrete','Bt','E_{ph}','E_{env}','\lambda','T','L_a','M','L_d','Pfe'};
er=[(a_mec(1)-a_fea(1))/a_fea(1), (mean(abs(R.Br))-mean(abs(BrFu)))/mean(abs(BrFu)), ...
    (max(R.Br)-max(BrFu))/max(BrFu), (R.Btrms-sqrt(mean(BtF.^2)))/sqrt(mean(BtF.^2)), ...
    (Eph_pk-M.FEA.emf_ph)/M.FEA.emf_ph, (env_pk-M.FEA.emf_env)/M.FEA.emf_env, ...
    (FLa_pk-M.FEA.FLa_pk)/M.FEA.FLa_pk, (T_load-M.FEA.T_load)/M.FEA.T_load, ...
    (La-M.FEA.La)/M.FEA.La, (Mmut-(-2.132e-3))/(-2.132e-3), (Ld-M.FEA.Ld)/M.FEA.Ld, ...
    (Pfe-M.FEA.Pfe_nl)/M.FEA.Pfe_nl, (Ppm-M.FEA.Psolid_nl)/M.FEA.Psolid_nl]*100;
nmv{end+1}='P_{aim}';
b=bar(er,'FaceColor','flat'); grid on; ylabel('ecart MEC-FEA (%)');
set(gca,'XTickLabel',nmv); title('Scorecard');
for k=1:numel(er)
    if abs(er(k))<10, b.CData(k,:)=[0.2 0.7 0.2]; else, b.CData(k,:)=[0.85 0.5 0.1]; end
end
sgtitle('MEC BLDC 15/14 750W : DtN etendu couple au reseau stator vs FEA ANSYS');
saveas(gcf,fullfile(fileparts(mfilename('fullpath')),'FIG_master.png'));
fprintf('\n>> FIG_master.png sauvee.\n');

function i=sixstep(e,Iflat)
    n=numel(e); k=(0:n-1)/n*2*pi;
    psi=atan2(2*mean(e.*sin(k)),2*mean(e.*cos(k)));
    ke=mod(k-psi,2*pi); i=zeros(size(e));
    i(ke<pi/3|ke>5*pi/3)=Iflat; i(ke>2*pi/3&ke<4*pi/3)=-Iflat;
end
