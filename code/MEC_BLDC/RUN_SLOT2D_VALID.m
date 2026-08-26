%% RUN_SLOT2D_VALID  Bouche d'encoche : scalaire identifie vs sous-modele 2D
%  Le juge n'est pas la metrique d'harmoniques mais les GRANDEURS PREDITES.
%  On confronte les deux representations de l'ouverture d'encoche sur tout ce
%  que la FEA fournit : champ, back-EMF, couple, pertes fer, pertes aimant.
clear; clc;
M=machine_bldc(); Ns=M.Ns; p=M.p; Nm=M.Nm;
if isfile('kfringe_ident.mat'), S=load('kfringe_ident.mat'); kfr=S.kbest; else, kfr=0.75; end
PA=[1 -2 -15 3 14]; PB=[6 -7 -5 8 4]; PC=[11 -12 -10 13 9];
f=M.FEA.n_nl*Nm/120; Np=721;

% --- reference FEA ---
fea=M.FEA.dir;
d4=readmatrix(fullfile(fea,'magnetostique(Magnetic_loading)','Calculator Expressions Plot 4.tab'),'FileType','text','NumHeaderLines',1);
d2=readmatrix(fullfile(fea,'magnetostique(Magnetic_loading)','Calculator Expressions Plot 2.tab'),'FileType','text','NumHeaderLines',1);
angF=d4(:,2)*pi/180; BrF=d4(:,3); BtF=d2(:,2);
thu=linspace(0,2*pi,3601); thu(end)=[];
BrFu=interp1([angF;2*pi],[BrF;BrF(1)],thu,'linear','extrap');

% --- geometrie de pertes ---
Vt1=M.wst1*(M.hs0+M.hs1+M.hs2)*M.ls*M.Ki;
Vy1=(2*pi*(M.Rso-M.wsy/2)/Ns)*M.wsy*M.ls*M.Ki;
tau_yi=2*pi*(M.Rsi+M.hs0+M.hs1+M.hs2)/Ns;
tt=linspace(0,1/f,Np); dt=tt(2)-tt(1);
bl=@(B,V) deal(sum(M.Kh*f*((max(B,[],2)-min(B,[],2))/2).^2)*V, ...
               sum((M.KeFe/(2*pi^2))*mean(gradient(B,dt).^2,2))*V);

fprintf('=== BOUCHE D''ENCOCHE : scalaire identifie vs sous-modele 2D ===\n');
res=struct();
for c=1:2
    if c==1
        nom='(a) kfringe = 0.75 identifie';
        Re=cogging_mec (M,840,0,Np,1500,kfr,2*pi/p);     % periode electrique
        Rs_=cogging_mec(M,840,0,181,1500,kfr,2*pi/Ns);   % pas d'encoche
    else
        nom='(b) bouche 2D maillee (nm=8)';
        Re=cogging_mec2(M,840,8,Np,1500,2*pi/p);
        Rs_=cogging_mec2(M,840,8,181,1500,2*pi/Ns);
    end
    % back-EMF depuis les flux de dent
    lam=@(P) M.Ntc*sum(sign(P(:)).*Re.PhiT(abs(P(:)),:),1);
    om=M.speed*2*pi/60;
    eA=-gradient(lam(PA),Re.phis)*om; eB=-gradient(lam(PB),Re.phis)*om; eC=-gradient(lam(PC),Re.phis)*om;
    env=max([eA;eB;eC],[],1)-min([eA;eB;eC],[],1);
    Iflat=M.Iph_rms_load*sqrt(3/2);
    sixs=@(e) sixstep(e,Iflat);
    T=mean((eA.*sixs(eA)+eB.*sixs(eB)+eC.*sixs(eC))/om);
    % pertes fer vectorielles
    Bt_=Re.PhiT/(M.wst1*M.ls*M.Ki); By_=Re.PhiY/(M.wsy*M.ls*M.Ki); Byr_=Re.PhiT/(tau_yi*M.ls*M.Ki);
    [h1,e1]=bl(Bt_,Vt1); [h2,e2]=bl(By_,Vy1); [h3,e3]=bl(Byr_,Vy1);
    Pfe=h1+e1+h2+e2+h3+e3;
    % pertes aimant
    Ppm=pm_loss_R(M,Rs_);
    res(c).nom=nom; res(c).Bg1=Re.Bg1; res(c).Bmean=mean(abs(Re.Br)); res(c).Bpk=max(Re.Br);
    res(c).Btrms=Re.Btrms; res(c).Eph=max(abs(eA)); res(c).env=max(env);
    res(c).T=T; res(c).Pfe=Pfe; res(c).Ppm=Ppm;
end

FEA=[1.0746, mean(abs(BrFu)), max(BrFu), sqrt(mean(BtF.^2)), M.FEA.emf_ph, ...
     M.FEA.emf_env, M.FEA.T_load, M.FEA.Pfe_nl, M.FEA.Psolid_nl];
nom={'Bg1 (T)','B moyen (T)','B crete (T)','Bt RMS (T)','FEM phase (V)', ...
     'Enveloppe (V)','Couple (N.m)','Pertes fer (W)','Pertes aimant (W)'};

fprintf('\n  %-18s %10s %10s %10s %10s %10s\n','Grandeur','FEA','(a) scal.','ecart','(b) 2D','ecart');
for k=1:9
    va=[res(1).Bg1 res(1).Bmean res(1).Bpk res(1).Btrms res(1).Eph res(1).env res(1).T res(1).Pfe res(1).Ppm];
    vb=[res(2).Bg1 res(2).Bmean res(2).Bpk res(2).Btrms res(2).Eph res(2).env res(2).T res(2).Pfe res(2).Ppm];
    fprintf('  %-18s %10.4f %10.4f %+9.1f%% %10.4f %+9.1f%%\n',nom{k},FEA(k), ...
            va(k),(va(k)-FEA(k))/FEA(k)*100, vb(k),(vb(k)-FEA(k))/FEA(k)*100);
end
ea=abs(([res(1).Bg1 res(1).Bmean res(1).Bpk res(1).Btrms res(1).Eph res(1).env res(1).T res(1).Pfe res(1).Ppm]-FEA)./FEA)*100;
eb=abs(([res(2).Bg1 res(2).Bmean res(2).Bpk res(2).Btrms res(2).Eph res(2).env res(2).T res(2).Pfe res(2).Ppm]-FEA)./FEA)*100;
fprintf('\n  |ecart| moyen : (a) %.1f %%   |   (b) %.1f %%\n',mean(ea),mean(eb));
fprintf('  |ecart| moyen hors pertes aimant : (a) %.1f %%   |   (b) %.1f %%\n',mean(ea(1:8)),mean(eb(1:8)));

function i=sixstep(e,Iflat)
    n=numel(e); k=(0:n-1)/n*2*pi;
    psi=atan2(2*mean(e.*sin(k)),2*mean(e.*cos(k)));
    ke=mod(k-psi,2*pi); i=zeros(size(e));
    i(ke<pi/3|ke>5*pi/3)=Iflat; i(ke>2*pi/3&ke<4*pi/3)=-Iflat;
end
