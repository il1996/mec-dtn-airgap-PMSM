%% RUN_SLOT  Couplage DtN etendu <-> reseau statorique : denture + Bt (ameli. 1)
%  Surface stator = 15 faces de dent (fer) + 15 ouvertures d'encoche (air).
%  Reseau : dent (radiale) -> culasse (anneau) ; ouverture -> becs vers dents
%  adjacentes. Systeme (K_reseau - Y_crown) U = Isrc.  Champ resolu -> denture.
clear; clc;
M = machine_bldc(); mu0=4*pi*1e-7;
Rs=M.Rsi; L=M.ls; Ns=M.Ns; p=M.p; taus=2*pi/Ns;

% --- 0. sanity : le champ aimant seul (Us=0) doit reproduire magnet_subdomain
wop=M.ws0/Rs; wtf=taus-wop;
ths=zeros(1,2*Ns); dths=ths;
for i=1:Ns
    ths(2*i-1)=(i-1)*taus;        dths(2*i-1)=wtf;    % face de dent
    ths(2*i)  =(i-1)*taus+taus/2; dths(2*i)  =wop;    % ouverture d'encoche
end
AG = airgap_magnet(M, ths, dths, 120);
thq = linspace(0,2*pi,3601);
[Br0,~] = AG.field(zeros(2*Ns,1), 0, thq);
fprintf('Sanity (Us=0) : Br moy %.3f T, Bg1 %.3f T  (attendu 0.77 / 1.09)\n', ...
        mean(abs(Br0)), hypot(2*mean(Br0.*cos(p*thq)),2*mean(Br0.*sin(p*thq))));

% --- 1. reseau de reluctances statorique (conductances = permeances) ---
muI = 1500;                                   % permeabilite relative fer (no-load)
Gt   = mu0*muI*M.wst1*L/(M.hs0+M.hs1+M.hs2);  % dent (radiale)
Gy   = mu0*muI*M.wsy*L/(taus*(M.Rso-M.wsy/2));% segment de culasse
Gtip = mu0*(M.hs0+M.hs1)*L/(M.ws0/2);         % bec (air, tangentiel)
Nn=3*Ns;                                      % 15 T + 15 O + 15 Y
Ti=@(i)2*i-1; Oi=@(i)2*i; Yi=@(i)2*Ns+i;
K=zeros(Nn);
addb=@(K,a,b,G) K + sparse([a a b b],[a b a b],[G -G -G G],Nn,Nn);
K=sparse(K);
for i=1:Ns
    j=mod(i,Ns)+1;                             % dent suivante
    K=addb(K,Ti(i),Yi(i),Gt);                 % dent -> culasse
    K=addb(K,Yi(i),Yi(j),Gy);                 % culasse -> culasse
    K=addb(K,Oi(i),Ti(i),Gtip);               % ouverture -> dent i
    K=addb(K,Oi(i),Ti(j),Gtip);               % ouverture -> dent i+1
end

% --- 2. crown : Y (surface) + Isrc ; assemblage plein ---
Yfull=zeros(Nn); Yfull(1:2*Ns,1:2*Ns)=AG.Y;
Isrc=[AG.Isrc(0); zeros(Ns,1)];
Asys = K - Yfull;                             % (flux reseau) = (flux crown)
% reference : masse sur un noeud de culasse (Isrc de moyenne nulle)
ref=Yi(1); keep=true(Nn,1); keep(ref)=false;
U=zeros(Nn,1);
U(keep) = Asys(keep,keep)\Isrc(keep);

% --- 3. champ resolu -> denture ---
[Br,Bt] = AG.field(U(1:2*Ns), 0, thq);

Bmean=mean(abs(Br)); Bpeak=max(Br);
Bg1=hypot(2*mean(Br.*cos(p*thq)),2*mean(Br.*sin(p*thq)));
Bt_rms=sqrt(mean(Bt.^2));

% --- FEA ---
fea=M.FEA.dir;
d4=readmatrix(fullfile(fea,'magnetostique(Magnetic_loading)','Calculator Expressions Plot 4.tab'),'FileType','text','NumHeaderLines',1);
d2=readmatrix(fullfile(fea,'magnetostique(Magnetic_loading)','Calculator Expressions Plot 2.tab'),'FileType','text','NumHeaderLines',1);
BrF=d4(:,3); angF=d4(:,2)*pi/180; BtF=d2(:,2);
Bt_rms_F=sqrt(mean(BtF.^2));

fprintf('\n=== DtN COUPLE AU RESEAU STATOR : denture emergente vs FEA ===\n');
r=@(n,a,f) fprintf('  %-26s %9.3f %9.3f %+8.1f%%\n',n,a,f,(a-f)/f*100);
r('B entrefer moyen |Br| (T)', Bmean, mean(abs(BrF)));
r('B entrefer crete (T)',      Bpeak, max(BrF));
r('Fondamental Bg1 (T)',       Bg1,   1.0749);
r('Bt tangentiel RMS (T)',     Bt_rms,Bt_rms_F);
r('Bt tangentiel crete (T)',   max(abs(Bt)), max(abs(BtF)));

figure('Name','Denture emergente','Color','w','Position',[80 80 1000 400]);
subplot(1,2,1);
plot(thq*180/pi,Br,'b','LineWidth',1.3); hold on;
plot(angF*180/pi,BrF,'r--','LineWidth',1.0); grid on; xlim([0 120]);
xlabel('\theta mec (deg)'); ylabel('B_r (T)'); title('B_r : denture MEC vs FEA (zoom)');
legend('MEC couple','FEA','Location','south');
subplot(1,2,2);
plot(thq*180/pi,Bt,'b','LineWidth',1.3); hold on;
plot(angF(1:numel(BtF))*180/pi,BtF,'r--','LineWidth',1.0); grid on; xlim([0 120]);
xlabel('\theta mec (deg)'); ylabel('B_t (T)'); title('B_t tangentiel MEC vs FEA (zoom)');
saveas(gcf,fullfile(fileparts(mfilename('fullpath')),'FIG_slot.png'));
fprintf('  (figure FIG_slot.png sauvee)\n');
