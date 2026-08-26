%% RUN_BEMF  Back-EMF a vide par fonction de bobinage (flux de dent) vs FEA.
%  Flux de dent i = integrale de Br(bore) sur le pas d'encoche, calcule
%  analytiquement a partir des harmoniques du champ de sous-domaine ; flux
%  de phase = Ntc * somme signee des dents de la phase (etoile d'encoches) ;
%  back-EMF = -dlambda/dphi * omega.
clear; clc;
M = machine_bldc();
Rs = M.Rsi; L = M.ls; Ns = M.Ns; p = M.p; Ntc = M.Ntc;
taus = 2*pi/Ns;                                  % pas d'encoche mec

% --- champ radial au bore, harmoniques (champ pair -> cosinus) ---
th = linspace(0,2*pi,7201); th(end)=[];
F  = magnet_subdomain(M, th, 15);
Brb = F.Br_bore(:).';
nord = 1:2:29;  nu = nord*p;                     % ordres spatiaux presents
Brn = arrayfun(@(v) 2*mean(Brb.*cos(v*th)), nu); % amplitudes au bore
apert = (2./nu).*sin(nu*taus/2);                 % ouverture de dent (sinc)

% --- assignation des dents par phase (winding.m, 15/14) ---
PA=[1 -2 -15 3 14]; PB=[6 -7 -5 8 4]; PC=[11 -12 -10 13 9];

% --- balayage sur 1 periode electrique ---
phi = linspace(0, 2*pi/p, 2001);
lamA = flux_phase(PA, phi, Brn, apert, nu, Ntc, Rs, L, Ns);
lamB = flux_phase(PB, phi, Brn, apert, nu, Ntc, Rs, L, Ns);
lamC = flux_phase(PC, phi, Brn, apert, nu, Ntc, Rs, L, Ns);

om = M.speed*2*pi/60;                            % vitesse mecanique rad/s
eA = -gradient(lamA, phi)*om;
eB = -gradient(lamB, phi)*om;
eC = -gradient(lamC, phi)*om;
env = max([eA;eB;eC],[],1) - min([eA;eB;eC],[],1);   % enveloppe DC-link six-step

Eph_pk  = max(abs(eA));
FLa_pk  = max(abs(lamA));
env_pk  = max(env);  env_mean = mean(env);

fprintf('\n=== BACK-EMF A VIDE @ %g tr/min : MEC vs FEA ===\n', M.speed);
r=@(n,a,f,u) fprintf('  %-30s %9.3f %9.3f %+8.1f%%  %s\n',n,a,f,(a-f)/f*100,u);
r('FEM de phase crete (V)',      Eph_pk,   M.FEA.emf_ph, '');
r('Enveloppe six-step crete (V)',env_pk,   M.FEA.emf_env,'');
r('Enveloppe six-step moy (V)',  env_mean, 444.6,        '(FEA moy)');
r('Flux totalise phase crete (Wb)',FLa_pk, M.FEA.FLa_pk, '');

figure('Name','Back-EMF MEC','Color','w');
subplot(2,1,1);
plot(phi*p*180/pi,eA,phi*p*180/pi,eB,phi*p*180/pi,eC,'LineWidth',1.2); hold on;
plot(phi*p*180/pi,env,'k:','LineWidth',1.3); grid on;
xlabel('angle elec (deg)'); ylabel('FEM (V)');
legend('e_A','e_B','e_C','enveloppe'); title(sprintf('Back-EMF MEC (phase crete %.1f V, env %.1f V)',Eph_pk,env_pk));
subplot(2,1,2);
plot(phi*p*180/pi,lamA,phi*p*180/pi,lamB,phi*p*180/pi,lamC,'LineWidth',1.2); grid on;
xlabel('angle elec (deg)'); ylabel('\lambda (Wb)'); title('Flux totalise de phase');
saveas(gcf, fullfile(fileparts(mfilename('fullpath')),'FIG_bemf.png'));

% ===================== fonction locale =================================
function lam = flux_phase(P, phi, Brn, apert, nu, Ntc, Rs, L, Ns)
    taus = 2*pi/Ns;
    lam = zeros(size(phi));
    for c = P
        i = abs(c); s = sign(c);
        thc = (i-1)*taus;
        lam = lam + s*Ntc*Rs*L*sum((Brn.*apert).' .* cos(nu.'*(thc - phi)), 1);
    end
end
