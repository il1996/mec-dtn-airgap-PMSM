%% RUN_PMLOSS_LOAD - optimisation des pertes aimant EN CHARGE
%  Cible : 3.183 W (FEA). Controle prealable fait sur BLDC.aedt :
%  'Eddy Effect'=true porte sur 14 objets EXACTEMENT, les 14 aimants ->
%  SolidLoss est bien la perte aimant seule, la comparaison est homogene.
%  Le premier calcul donnait 1.366 W (-57 %). On teste une par une les
%  causes possibles de sous-estimation :
%    (1) resolution temporelle (les fronts de commutation durent ~50 us) ;
%    (2) derivee de ROTATION par difference finie au lieu d'analytique ;
%    (3) frange kfringe des potentiels de surface du champ d'INDUIT ;
%    (4) deficit connu de l'harmonique de denture 8 (deja documente a vide).
clear; clc; t0=tic;
M=machine_bldc(); p=M.p; Ntc=M.Ntc; fea=M.FEA.dir;
PA=[1 -2 -15 3 14]; PB=[6 -7 -5 8 4]; PC=[11 -12 -10 13 9];
if isfile('kfringe_ident.mat'), S0=load('kfringe_ident.mat'); kfr=S0.kbest; else, kfr=0.75; end
rd=@(f)readmatrix(f,'FileType','text','Delimiter','\t','NumHeaderLines',1);
ol=fullfile(fea,'transitoire (en charge)'); w4=@(x)x(round(numel(x)*0.75):end);
dLo=rd(fullfile(ol,'Plot 1_loss.tab')); PpmF=mean(w4(dLo(:,2)))/1e3;

%% ---- transitoire MEC ----
R=cogging_mec(M,1260,0,721,M.muI,kfr,2*pi/p);
the=R.phis*p; lamf=@(P) Ntc*sum(sign(P(:)).*R.PhiT(abs(P(:)),:),1);
lam=[lamf(PA); lamf(PB); lamf(PC)];
RI=inductance_mec(M,1260,M.muI,0);
d4=rd(fullfile(fea,'transitoire (Back_emf)','Output Variables Plot 4.tab'));
cLM=2*mean(lam(1,:).*exp(-1i*the));
cLF=2*mean(d4(:,3).'.*exp(-1i*d4(:,2).'*p*pi/180));
dphi=-angle(cLM/cLF);
Q=load('mec_map.mat'); fn=fieldnames(Q); MAP=Q.(fn{1});
opt=struct('dt',1e-6,'tend',34.26e-3,'J',1e-3,'Bf',6.07927101854027e-4, ...
    'Tload',4.8,'Vdc',M.Vdc,'Rph',M.Rph,'dphi',dphi,'map',MAP);
D=drive_mec(M,the,lam,RI.Ld,opt);
ws=D.t>0.75*D.t(end); nM=mean(D.n(ws)); Te=60/(nM*M.Nm/2);
mw=D.t>D.t(end)-Te; tw=D.t(mw); iw=D.i(:,mw); thw=D.th(mw);
fprintf('=== Pertes aimant EN CHARGE : optimisation ===\n');
fprintf('  reference FEA (SolidLoss, 14 aimants) : %.3f W\n',PpmF);
fprintf('  regime MEC %.0f tr/min, periode electrique %.3f ms, %d pas dispo\n', ...
    nM,Te*1e3,numel(tw));

%% ---- (1) resolution temporelle --------------------------------------
fprintf('\n  --- (1) resolution temporelle (kfringe induit = 0) ---\n');
fprintf('  %8s %12s %12s %10s\n','Nt','pas (us)','P (W)','ecart');
for Nt=[121 241 481 961 1921]
    js=round(linspace(1,numel(tw),Nt));
    Pm=pm_loss_load(M,R,RI.Usurf,tw(js),iw(:,js),thw(js));
    fprintf('  %8d %12.2f %12.4f %9.1f %%\n',Nt,(tw(js(2))-tw(js(1)))*1e6, ...
        Pm,100*(Pm-PpmF)/PpmF);
end

%% ---- (2)-(3) frange du champ d'induit -------------------------------
%  kfringe = 0 est le bon choix pour l'INDUCTANCE (il evite de compter deux
%  fois la fuite de bec) mais pas pour le CHAMP qui traverse l'entrefer :
%  la, la valeur identifiee sur les bandes de denture FEA est 0.75.
js=round(linspace(1,numel(tw),1921));
fprintf('\n  --- (2) frange du champ d''induit (Nt = 1921) ---\n');
fprintf('  %14s %12s %10s\n','kfringe induit','P (W)','ecart');
Pk=zeros(1,3); kv=[0 0.375 0.75];
for q=1:3
    RIk=inductance_mec(M,1260,M.muI,kv(q));
    Pk(q)=pm_loss_load(M,R,RIk.Usurf,tw(js),iw(:,js),thw(js));
    fprintf('  %14.3f %12.4f %9.1f %%\n',kv(q),Pk(q),100*(Pk(q)-PpmF)/PpmF);
end

%% ---- (3) decomposition ----------------------------------------------
RIk=inductance_mec(M,1260,M.muI,kfr);
Ptot=Pk(3);
Pden=pm_loss_load(M,R,RIk.Usurf,tw(js),iw(:,js)*0,thw(js));   % denture seule
Zs=R; Zs.Usurf=R.Usurf*0;
Pind=pm_loss_load(M,Zs,RIk.Usurf,tw(js),iw(:,js),thw(js));    % induit seul
fprintf('\n  --- (3) decomposition de la perte (kfringe = %.2f) ---\n',kfr);
fprintf('  denture seule (courants annules)     : %.4f W\n',Pden);
fprintf('  induit seul (aimants eteints)        : %.4f W\n',Pind);
fprintf('  total (avec interference)            : %.4f W   (FEA %.3f)\n',Ptot,PpmF);
fprintf('  controle : pm_loss a vide a %.0f tr/min = %.4f W\n', ...
    nM,pm_loss(M,840,91,kfr,nM));

%% ---- (4) deficit connu de l'harmonique de denture --------------------
%  La perte varie en B^2 et l'harmonique dominant vu par l'aimant est
%  nu = |Ns-p| = 8, le moins attenue radialement. Son deficit est deja
%  documente sur le champ d'entrefer a vide.
dm4=rd(fullfile(fea,'magnetostique(Magnetic_loading)','Calculator Expressions Plot 4.tab'));
thu=linspace(0,2*pi,3601); thu(end)=[];
BrFu=interp1([dm4(:,2)*pi/180;2*pi],[dm4(:,3);dm4(1,3)],thu,'linear','extrap');
ampf=@(B,th,k) 2*abs(mean(B(:).'.*exp(-1i*k*th(:).')));
a8m=ampf(R.Br(1:end-1),R.thq(1:end-1),8); a8f=ampf(BrFu,thu,8);
fprintf('\n  --- (4) deficit de l''harmonique 8 (deja documente a vide) ---\n');
fprintf('  harmonique 8 du champ d''entrefer : MEC %.5f T | FEA %.5f T (%+.1f %%)\n', ...
    a8m,a8f,100*(a8m-a8f)/a8f);
fprintf('  la part DENTURE de la perte varie en B8^2 -> facteur %.2f\n',(a8f/a8m)^2);
Pcor=Pden*(a8f/a8m)^2+(Ptot-Pden);
fprintf('  perte corrigee de ce seul deficit : %.4f W  (%+.1f %% vs FEA)\n', ...
    Pcor,100*(Pcor-PpmF)/PpmF);
fprintf('\n  (%.0f s)\n',toc(t0));
