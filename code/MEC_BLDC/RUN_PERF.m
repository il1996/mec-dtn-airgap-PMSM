%% RUN_PERF  Inductance, couple en charge, pertes fer : MEC vs FEA.
clear; clc;
M = machine_bldc(); mu0=4*pi*1e-7;
Rs=M.Rsi; L=M.ls; Ns=M.Ns; p=M.p; Ntc=M.Ntc; taus=2*pi/Ns;
PA=[1 -2 -15 3 14]; PB=[6 -7 -5 8 4]; PC=[11 -12 -10 13 9];

% ===================== champ & harmoniques ============================
th = linspace(0,2*pi,7201); th(end)=[];
F  = magnet_subdomain(M, th, 15);
Brb = F.Br_bore(:).';
nord=1:2:29; nu=nord*p;
Brn = arrayfun(@(v) 2*mean(Brb.*cos(v*th)), nu);
apert=(2./nu).*sin(nu*taus/2);

% ===================== 1. INDUCTANCE ==================================
%  L_airgap par fonction de bobinage (capte magnetisant + sous-harmonique
%  ET differentielle des 15/14) ; + fuite d'encoche Pyrhonen ; entrefer
%  magnetique effectif = kc*(lag + hm/mu_r) (aimant vu comme de l'air).
gmag = M.lag + M.hm/M.mu_r;
b0=M.ws0; ub=b0/gmag; gammaC=ub^2/(5+ub); kc=taus*Rs/(taus*Rs-gammaC*gmag);
geff = kc*gmag;
% winding function de la phase A (dents-bobines)
thg = linspace(0,2*pi,3600); dth=thg(2)-thg(1);
NA = zeros(size(thg));
for c=PA
    i=abs(c); s=sign(c); c0=(i-1)*taus;
    win = abs(angle(exp(1j*(thg-c0)))) <= taus/2;   % pas d'encoche de la dent
    NA = NA + s*Ntc*win;
end
WA = NA - mean(NA);                                  % fonction de bobinage
L_airgap = mu0*Rs*L/geff * sum(WA.^2)*dth;           % self airgap [H]
% fuite d'encoche (Pyrhonen, double couche) + bec + tetes de bobines
bs_avg=(M.ws1+M.ws2)/2;
lam_slot = M.hs2/(3*bs_avg) + M.hs1/(M.ws1) + M.hs0/M.ws0;   % corps+biseau+ouverture
lam_tip  = (5*(M.lag/b0))/(5+4*(M.lag/b0));
Ncp = numel(PA);  Nph = Ntc*Ncp;                     % bobines/phase, spires serie/phase
L_slot = (4*M.m/Ns)*mu0*L*Nph^2*(lam_slot+lam_tip);  % Pyrhonen (phase, 4m/Q)
La_MEC = L_airgap + L_slot;
Ld_MEC = La_MEC;                                     % mutuelle ~faible (FSCW)

fprintf('\n=== INDUCTANCE : MEC vs FEA ===\n');
fprintf('  L_airgap (fct bobinage)      = %6.2f mH\n', L_airgap*1e3);
fprintf('  L_fuite encoche+bec          = %6.2f mH  (lam_slot=%.3f)\n', L_slot*1e3, lam_slot);
fprintf('  L_propre La : MEC %6.2f mH | FEA %6.2f mH  (%+.1f %%)\n', ...
        La_MEC*1e3, M.FEA.La*1e3, (La_MEC-M.FEA.La)/M.FEA.La*100);

% ===================== 2. COUPLE EN CHARGE (six-step) =================
phi = linspace(0,2*pi/p,3601);
lamA=flux_phase(PA,phi,Brn,apert,nu,Ntc,Rs,L,Ns);
lamB=flux_phase(PB,phi,Brn,apert,nu,Ntc,Rs,L,Ns);
lamC=flux_phase(PC,phi,Brn,apert,nu,Ntc,Rs,L,Ns);
om=M.speed*2*pi/60;
eA=-gradient(lamA,phi)*om; eB=-gradient(lamB,phi)*om; eC=-gradient(lamC,phi)*om;
Iflat = M.Iph_rms_load*sqrt(3/2);                    % rms 120deg -> plateau
% courant six-step aligne sur le fondamental de chaque back-EMF
iA=sixstep(eA,Iflat); iB=sixstep(eB,Iflat); iC=sixstep(eC,Iflat);
Tinst=(eA.*iA+eB.*iB+eC.*iC)/om;
T_pow=mean(Tinst);
fprintf('\n=== COUPLE EN CHARGE (I_ph=%.2f A eff, plateau %.2f A) : MEC vs FEA ===\n',...
        M.Iph_rms_load,Iflat);
fprintf('  Couple (puissance sum e.i/w) : MEC %6.3f N.m | FEA %6.3f N.m  (%+.1f %%)\n',...
        T_pow, M.FEA.T_load, (T_pow-M.FEA.T_load)/M.FEA.T_load*100);

% ===================== 3. PERTES FER A VIDE (Bertotti B(t)) ===========
%  B(t) de dent et de culasse sur 1 periode elec ; Bertotti volumique M350.
f = M.FEA.n_nl*M.Nm/120;                              % f a la vitesse FEA a vide
% flux de dent i sur une periode elec -> B_dent(t)
phie = linspace(0,2*pi,2001); phr=phie/p;             % 1 periode elec = 2pi/p mec
PhiT = Rs*L*sum((Brn.*apert).' .* cos(nu.'*(0 - phr)),1);  % dent 1
Bt_t = PhiT/(M.wst1*L*M.Ki);
% culasse stator : porte ~ Phi_pole/2 ; Phi_pole a partir du fondamental
Phi_pole = 2*Rs*L*Brn(1)/nu(1)*2;                     % integrale demi-onde ~ (2/nu)Brn1 * Rs L *2
Bys_t = 0.5*Phi_pole/(M.wsy*L*M.Ki) * cos(p*phr);     % alternance a f
Vt = M.Ns*M.wst1*M.hs*L*M.Ki;                         % volume dents
Rys = M.Rso - M.wsy/2; Vy = 2*pi*Rys*M.wsy*L*M.Ki;    % volume culasse
Pfe = bertotti_loss(Bt_t,f,M)*Vt + bertotti_loss(Bys_t,f,M)*Vy;
fprintf('\n=== PERTES FER A VIDE @ %g tr/min (f=%.0f Hz) : MEC vs FEA ===\n',M.FEA.n_nl,f);
fprintf('  Pertes fer : MEC %6.2f W | FEA %6.2f W  (%+.1f %%)\n',...
        Pfe, M.FEA.Pfe_nl, (Pfe-M.FEA.Pfe_nl)/M.FEA.Pfe_nl*100);
fprintf('  (B dent crete %.3f T, B culasse crete %.3f T)\n',max(abs(Bt_t)),max(abs(Bys_t)));

% ===================== fonctions locales =============================
function lam=flux_phase(P,phi,Brn,apert,nu,Ntc,Rs,L,Ns)
    taus=2*pi/Ns; lam=zeros(size(phi));
    for c=P, i=abs(c); s=sign(c); thc=(i-1)*taus;
        lam=lam+s*Ntc*Rs*L*sum((Brn.*apert).'.*cos(nu.'*(thc-phi)),1); end
end
function i=sixstep(e,Iflat)
    % fondamental de e -> phase ; courant flat 120deg en phase avec e
    n=numel(e); k=(0:n-1)/n*2*pi;                    % angle elec sur la periode
    c=2*mean(e.*cos(k)); s=2*mean(e.*sin(k)); psi=atan2(s,c);
    ke=mod(k-psi,2*pi);
    i=zeros(size(e));
    i(ke<pi/3 | ke>5*pi/3)=Iflat;                    % +120deg centre sur le pic
    i(ke>2*pi/3 & ke<4*pi/3)=-Iflat;
end
function pv=bertotti_loss(Bt,f,M)
    % Bertotti volumique : hyst (crete) + Foucault (dB/dt), Kex=0
    dt=1/(f*(numel(Bt)-1)); % pas de temps sur 1 periode
    tt=linspace(0,1/f,numel(Bt)); dBdt=gradient(Bt,tt);
    Bpk=max(abs(Bt));
    pv = M.Kh*f*Bpk^2 + (M.KeFe/(2*pi^2))*mean(dBdt.^2);
end
