%% RUN_DIAG_LOAD - validation du modele d'entrainement A VITESSE IMPOSEE
%  On isole l'electromagnetisme de la mecanique : rotor entraine a la
%  vitesse mesuree par la FEA (1339.8 tr/min), meme onduleur, meme calage.
%  Si le MEC est bon, il doit retrouver le courant, le couple et la FEM.
clear; clc;
M=machine_bldc(); p=M.p; Ntc=M.Ntc; fea=M.FEA.dir;
PA=[1 -2 -15 3 14]; PB=[6 -7 -5 8 4]; PC=[11 -12 -10 13 9];
if isfile('kfringe_ident.mat'), S0=load('kfringe_ident.mat'); kfr=S0.kbest; else, kfr=0.75; end
rd=@(f)readmatrix(f,'FileType','text','Delimiter','\t','NumHeaderLines',1);
ol=fullfile(fea,'transitoire (en charge)'); w4=@(x)x(round(numel(x)*0.75):end);

R=cogging_mec(M,1260,0,721,M.muI,kfr,2*pi/p);
the=R.phis*p; lamf=@(P) Ntc*sum(sign(P(:)).*R.PhiT(abs(P(:)),:),1);
lam=[lamf(PA); lamf(PB); lamf(PC)];
RI=inductance_mec(M,1260,5000,0); Leff=RI.Ld;
dFL=rd(fullfile(fea,'transitoire (Back_emf)','Output Variables Plot 4.tab'));
cLM=2*mean(lam(1,:).*exp(-1i*the));
cLF=2*mean(dFL(:,3).'.*exp(-1i*dFL(:,2).'*p*pi/180));
dphi0=-angle(cLM/cLF);
Q=load('mec_map.mat'); MAP=Q.S;

dS=rd(fullfile(ol,'Speed Plot 1.tab')); dI=rd(fullfile(ol,'BranchCurrent Plot 1.tab'));
dL=rd(fullfile(ol,'Plot 1_loss.tab')); dIt=rd(fullfile(ol,'BranchCurrent Plot 2.tab'));
nF=mean(w4(dS(:,2))); omF=nF*2*pi/60; TF=mean(w4(dL(:,7)))/omF;
IF=sqrt(mean(w4(dI(:,2)).^2)); IbF=mean(w4(dIt(:,2)));
fprintf('=== Reference FEA en regime : %.1f tr/min, T = %.4f N.m, I = %.4f A ===\n',nF,TF,IF);

%% ---- balayage du calage : le couple doit etre maximal au bon angle ----
fprintf('\n--- sensibilite au calage (vitesse imposee %.0f tr/min) ---\n',nF);
fprintf('  %10s %10s %10s %10s %10s\n','dphi(deg)','T (N.m)','I rms (A)','I bus(A)','ecart T');
dd=-180:15:165; Tm=zeros(size(dd)); Im=Tm; Bm=Tm;
for k=1:numel(dd)
    o=struct('dt',2e-6,'tend',12e-3,'J',1e9,'Bf',0,'Tload',0,'Vdc',M.Vdc, ...
        'Rph',M.Rph,'dphi',dphi0+dd(k)*pi/180,'om0',omF,'map',MAP);
    D=drive_mec(M,the,lam,Leff,o); m=D.t>D.t(end)*0.5;
    Tm(k)=mean(D.T(m)); Im(k)=sqrt(mean(D.i(1,m).^2)); Bm(k)=mean(D.idc(m));
end
[~,kb]=min(abs(Tm-TF));
for k=1:numel(dd)
    mk=''; if k==kb, mk='  <-- meilleur'; end
    if mod(dd(k),45)==0 || k==kb
        fprintf('  %10.0f %10.3f %10.3f %10.3f %9.1f %%%s\n', ...
            dd(k),Tm(k),Im(k),Bm(k),100*(Tm(k)-TF)/TF,mk);
    end
end
[Tmx,kx]=max(Tm);
fprintf('  couple MAXIMAL %.3f N.m au calage %+.0f deg ; calage a vide -> %+.0f deg\n', ...
    Tmx,dd(kx),0);

%% ---- au calage etabli a vide : comparaison complete ----
o=struct('dt',1e-6,'tend',20e-3,'J',1e9,'Bf',0,'Tload',0,'Vdc',M.Vdc, ...
    'Rph',M.Rph,'dphi',dphi0,'om0',omF,'map',MAP);
D=drive_mec(M,the,lam,Leff,o); m=D.t>D.t(end)*0.5;
o.map=[]; Dl=drive_mec(M,the,lam,Leff,o);
fprintf('\n--- vitesse imposee %.0f tr/min, calage etabli a vide ---\n',nF);
fprintf('  %-24s %11s %11s %11s\n','','MEC lineaire','MEC sature','FEA');
pr=@(n,a,b,c,u)fprintf('  %-24s %11.4f %11.4f %11.4f  %s\n',n,a,b,c,u);
pr('couple',mean(Dl.T(m)),mean(D.T(m)),TF,'N.m');
pr('courant de phase rms',sqrt(mean(Dl.i(1,m).^2)),sqrt(mean(D.i(1,m).^2)),IF,'A');
pr('courant de bus moyen',mean(Dl.idc(m)),mean(D.idc(m)),IbF,'A');
pr('pertes Joule',mean(Dl.Pcu(m)),mean(D.Pcu(m)),mean(w4(dL(:,4))),'W');
pr('puissance Pem',mean(Dl.T(m))*omF,mean(D.T(m))*omF,mean(w4(dL(:,7))),'W');
