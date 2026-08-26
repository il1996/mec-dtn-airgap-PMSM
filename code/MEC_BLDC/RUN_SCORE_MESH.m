%% RUN_SCORE_MESH - le maillage raffine sur les grandeurs du scorecard
%
%  Le maillage ne servait qu'aux cartes. On le confronte ici, grandeur par
%  grandeur, au reseau a une dent ET a la FEA, sur le bloc A VIDE du
%  scorecard. Tant que ce bloc n'est pas au niveau, il n'y a pas lieu de
%  basculer le programme maitre : une carte fine mais un flux faux serait
%  un recul.
%
%  Balayage rotor sur UNE periode electrique, Newton demarre a chaud sur la
%  position precedente (2 a 3 iterations au lieu de 9).
clear; clc; t0=tic;
M=machine_bldc(); Ns=M.Ns; Nm=M.Nm; p=M.p; L=M.ls; Ki=M.Ki; Ntc=M.Ntc;
PA=[1 -2 -15 3 14]; PB=[6 -7 -5 8 4]; PC=[11 -12 -10 13 9];
if isfile('kfringe_ident.mat'), S0=load('kfringe_ident.mat'); kfr=S0.kbest; else, kfr=0.325; end
rd=@(f)readmatrix(f,'FileType','text','Delimiter','\t','NumHeaderLines',1);
fea=M.FEA.dir;
amp=@(B,th,k) 2*abs(mean(B(:).'.*exp(-1i*k*th(:).')));

%% ---------- references FEA --------------------------------------------
d4=rd(fullfile(fea,'magnetostique(Magnetic_loading)','Calculator Expressions Plot 4.tab'));
d2=rd(fullfile(fea,'magnetostique(Magnetic_loading)','Calculator Expressions Plot 2.tab'));
dA=rd(fullfile(fea,'magnetostique(Magnetic_loading)','Plot 1.tab'));
thu=linspace(0,2*pi,3601); thu(end)=[];
BrF=interp1([d4(:,2)*pi/180;2*pi],[d4(:,3);d4(1,3)],thu,'linear','extrap');
BtF=interp1([linspace(0,2*pi,numel(d2(:,2))).'],[d2(:,2)],thu,'linear','extrap');
AF =max(abs(dA(:,2)-mean(dA(:,2))));
dFL=rd(fullfile(fea,'transitoire (Back_emf)','Output Variables Plot 4.tab'));
dV =rd(fullfile(fea,'transitoire (Back_emf)','NodeVoltage Plot.tab'));
dE =rd(fullfile(fea,'transitoire (Back_emf)','Output Variables Plot.tab'));
FLaF=max(abs(dFL(:,3))); EaF=max(abs(dV(:,3))); envF=max(dE(:,3));

%% ---------- balayage rotor sur le maillage raffine ---------------------
Ms=540; Np=61; phis=linspace(0,2*pi/p,Np);
PhiT=zeros(Ns,Np); PhiY=zeros(Ns,Np); Us=zeros(Ms,Np); Tq=zeros(1,Np);
U0=[];
for q=1:Np
    ME=mesh_bldc(M,Ms,4,3,[],phis(q),kfr,'nl');
    S=solve_bldc_mesh(ME,1e-9,30,U0); U0=S.U;
    Us(:,q)=S.Usurf;
    %  flux de DENT : branches radiales de fer a mi-corps d'encoche
    mT=ME.type==1 & ME.lay==ME.lamT & ME.iron;
    PhiT(:,q)=accumarray(ME.tooth(ME.col(mT)).',S.Phi(mT),[Ns 1]);
    %  flux de CULASSE : branches tangentielles des couches de culasse,
    %  regroupees par secteur d'encoche
    mY=ME.type==2 & ismember(ME.lay,ME.layY);
    sec=mod(floor(ME.ths(ME.col(mY))/ME.taus+0.5),Ns)+1;
    PhiY(:,q)=accumarray(sec.',S.Phi(mY),[Ns 1]);
    %  couple : tenseur de Maxwell harmonique (meme forme que cogging_mec)
    AG=ME.AG; nu=AG.nu;
    Usc=AG.Wc*S.Usurf; Uss=AG.Ws*S.Usurf;
    cq=cos(nu*phis(q)); sq=sin(nu*phis(q));
    Brc=AG.bru.*Usc+AG.brmq.*cq;  Brs=AG.bru.*Uss+AG.brmq.*sq;
    Btc=-AG.btu.*Uss-AG.btmq.*sq; Bts=AG.btu.*Usc+AG.btmq.*cq;
    Tq(q)=(L*M.rmid^2/(4*pi*1e-7))*pi*sum(Brc.*Btc+Brs.*Bts);
end
fprintf('=== Maillage raffine sur le scorecard (bloc a vide) ===\n');
fprintf('  Ms=%d, %d positions, Newton a chaud (%.0f s)\n',Ms,Np,toc(t0));

%% ---------- grandeurs ---------------------------------------------------
[Br,Bt]=ME.AG.field(Us(:,1),phis(1),thu);
lamf=@(P) Ntc*sum(sign(P(:)).*PhiT(abs(P(:)),:),1);
lam=[lamf(PA); lamf(PB); lamf(PC)]; the=phis*p; om=M.speed*2*pi/60;
e=zeros(3,Np); for k=1:3, e(k,:)=-gradient(lam(k,:),phis)*om; end
env=max(e,[],1)-min(e,[],1);
Acum=cumtrapz(thu,Br)*M.rmid; Acum=Acum-mean(Acum);
Tc=Tq-mean(Tq);

%% ---------- reseau a UNE DENT (reference interne) ----------------------
R1=cogging_mec(M,1260,0,Np,M.muI,kfr,2*pi/p);
lam1=[Ntc*sum(sign(PA(:)).*R1.PhiT(abs(PA(:)),:),1); ...
      Ntc*sum(sign(PB(:)).*R1.PhiT(abs(PB(:)),:),1); ...
      Ntc*sum(sign(PC(:)).*R1.PhiT(abs(PC(:)),:),1)];
e1=zeros(3,Np); for k=1:3, e1(k,:)=-gradient(lam1(k,:),R1.phis)*om; end
env1=max(e1,[],1)-min(e1,[],1);
[Br1,Bt1]=R1.AG.field(R1.U(1:1260,1),0,thu);
A1=cumtrapz(thu,Br1)*M.rmid; A1=A1-mean(A1);

%% ---------- tableau -----------------------------------------------------
fprintf('\n  %-26s %10s %10s %10s %9s %9s\n', ...
    'GRANDEUR','1 dent','maillage','FEA','1 dent','maillage');
pr=@(n,a1,a2,f,u)fprintf('  %-26s %10.4f %10.4f %10.4f %8.1f %% %8.1f %%  %s\n', ...
    n,a1,a2,f,100*(a1-f)/f,100*(a2-f)/f,u);
pr('Bg1 fondamental',amp(Br1,thu,p),amp(Br,thu,p),amp(BrF,thu,p),'T');
pr('B entrefer moyen |Br|',mean(abs(Br1)),mean(abs(Br)),mean(abs(BrF)),'T');
pr('B entrefer crete',max(Br1),max(Br),max(BrF),'T');
pr('Bt tangentiel RMS',sqrt(mean(Bt1.^2)),sqrt(mean(Bt.^2)),sqrt(mean(BtF.^2)),'T');
pr('bande de denture (rang 8)',amp(Br1,thu,8),amp(Br,thu,8),amp(BrF,thu,8),'T');
pr('Lignes de flux A crete',max(abs(A1)),max(abs(Acum)),AF,'Wb/m');
pr('Flux totalise crete',max(abs(lam1(1,:))),max(abs(lam(1,:))),FLaF,'Wb');
pr('FEM de phase crete',max(abs(e1(1,:))),max(abs(e(1,:))),EaF,'V');
pr('Enveloppe six-step',max(env1),max(env),envF,'V');
%  PAS DE DETENTE ICI. Elle est d'ordre LCM(15,14) = 210 : avec 61 positions
%  sur une periode electrique la frequence de Nyquist est l'ordre 30, et les
%  deux modeles renvoient 0.000 par sous-echantillonnage. La comparer
%  demande un balayage dedie (>= 841 positions), ce que fait deja
%  cogging_mec pour le reseau a une dent ; le porter au maillage coute
%  841 resolutions de Newton, soit ~10 min.
fprintf('  detente : NON evaluee ici (ordre 210, Nyquist de ce balayage = 30)\n');
fprintf('\n  (%.0f s)\n',toc(t0));
save('score_mesh.mat','PhiT','PhiY','Us','Tq','phis','lam','e','-v7.3');
