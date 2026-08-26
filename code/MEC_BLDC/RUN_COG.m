%% RUN_COG  Couple de detente par DtN couple au reseau statorique (ameli. 1)
%  Grille de surface FINE (Nsurf noeuds) : faces de dent (fer -> culasse) et
%  ouvertures d'encoche (air, pas de puits de flux -> creux). A chaque
%  position rotor : source aimant Isrc(phi) -> reseau couple -> couple MST.
%  Couple MST calcule par HARMONIQUES (Parseval) : rapide, sans grille theta.
clear; clc;
M = machine_bldc(); mu0=4*pi*1e-7;
Rs=M.Rsi; L=M.ls; Ns=M.Ns; Nm=M.Nm; p=M.p; taus=2*pi/Ns;
NL=lcm(Ns,Nm); numax=250;

% --- grille de surface fine + classification fer/air ---
Nsurf=420; dth=2*pi/Nsurf;
ths=(0:Nsurf-1)*dth; dths=dth*ones(1,Nsurf);
wtf=taus-M.ws0/Rs;                                   % arc de face de dent
tooth=zeros(1,Nsurf); isFe=false(1,Nsurf);
for j=1:Nsurf
    i0=mod(round(ths(j)/taus),Ns);                   % dent la plus proche (0..Ns-1)
    d=angle(exp(1j*(ths(j)-i0*taus)));
    tooth(j)=i0+1; isFe(j)=abs(d)<wtf/2;             % fer si sous la face de dent
end
AG=airgap_magnet(M, ths, dths, numax);

% --- reseau statorique (fer) : face -> corps de dent -> culasse ---
muI=1500;
g_face=mu0*muI*(dth*Rs)*L/(M.hs0+M.hs1);            % shoe (par noeud de face)
Gt=mu0*muI*M.wst1*L/M.hs2;                           % corps de dent
Gy=mu0*muI*M.wsy*L/(taus*(M.Rso-M.wsy/2));           % culasse
Ntot=Nsurf+2*Ns; TB=@(i)Nsurf+i; YY=@(i)Nsurf+Ns+i;
I=[];J=[];V=[];
addb=@(a,b,G) deal([a a b b],[a b a b],[G -G -G G]);
for j=1:Nsurf
    if isFe(j)                                       % face de dent -> corps
        [ii,jj,vv]=addb(j,TB(tooth(j)),g_face); I=[I ii];J=[J jj];V=[V vv];
    end                                              % ouverture : aucun puits fer
end
for i=1:Ns
    j=mod(i,Ns)+1;
    [ii,jj,vv]=addb(TB(i),YY(i),Gt); I=[I ii];J=[J jj];V=[V vv];
    [ii,jj,vv]=addb(YY(i),YY(j),Gy); I=[I ii];J=[J jj];V=[V vv];
end
K=sparse(I,J,V,Ntot,Ntot);
Yfull=sparse(1:Nsurf,1:Nsurf,0,Ntot,Ntot); Yfull(1:Nsurf,1:Nsurf)=AG.Y;
Asys=K-Yfull; ref=YY(1); keep=true(Ntot,1); keep(ref)=false;

% --- source aimant pour toutes les positions (vectorise) ---
Np=421; phis=linspace(0,2*pi/Nm,Np);                 % 1 pas polaire
nu=AG.nu; brm=AG.brm;
Isrc_all=zeros(Ntot,Np);
Isrc_all(1:Nsurf,:)=L*Rs*pi*( AG.Wc.'*(brm.*cos(nu*phis)) + AG.Ws.'*(brm.*sin(nu*phis)) );
U=zeros(Ntot,Np);
U(keep,:)=Asys(keep,keep)\Isrc_all(keep,:);

% --- couple MST par harmoniques (Parseval) : T=(L r^2/mu0)*pi*sum(Brc Btc+Brs Bts)
Us=U(1:Nsurf,:);                                     % potentiels de surface
Usc=AG.Wc*Us; Uss=AG.Ws*Us;                          % (numax x Np)
bru=AG.bru; btu=AG.btu; brmq=AG.brmq; btmq=AG.btmq;
cphi=cos(nu*phis); sphi=sin(nu*phis);
Brc=bru.*Usc + brmq.*cphi;   Brs=bru.*Uss + brmq.*sphi;
Btc=-btu.*Uss - btmq.*sphi;  Bts=btu.*Usc + btmq.*cphi;
Tcog=(L*M.rmid^2/mu0)*pi*sum(Brc.*Btc + Brs.*Bts, 1);
Tcog=Tcog-mean(Tcog);
Tpp=(max(Tcog)-min(Tcog))*1e3;

% ordre dominant (cycles par pas polaire ; attendu NL/Nm=15)
Yf=abs(fft(Tcog)); Yf=Yf(2:floor(Np/2)); [~,km]=max(Yf);
fprintf('\n=== COUPLE DE DETENTE : DtN couple (grille fine %d) vs FEA ===\n',Nsurf);
fprintf('  Amplitude c-c    : MEC %6.2f mN.m | FEA %6.2f mN.m (sous-echantillonnee)\n',Tpp,M.FEA.cog_pp);
fprintf('  Cycles/pas polaire : %d  (attendu %d = LCM/Nm)  -> ordre %d/tour (LCM=%d)\n',...
        km, NL/Nm, km*Nm, NL);

figure('Name','Couple de detente (grille fine)','Color','w');
plot(phis*180/pi,Tcog*1e3,'b','LineWidth',1.3); grid on;
xlabel('position rotor (deg mec)'); ylabel('T_{detente} (mN.m)');
title(sprintf('Couple de detente MEC couple (c-c %.2f mN.m, %d cycles/pas polaire, ordre %d)',Tpp,km,km*Nm));
saveas(gcf,fullfile(fileparts(mfilename('fullpath')),'FIG_cog.png'));
fprintf('  (figure FIG_cog.png sauvee)\n');
