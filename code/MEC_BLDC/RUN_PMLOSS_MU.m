%% RUN_PMLOSS_MU - sensibilite des pertes aimant a la permeabilite d'induit
%  Le champ d'induit vu par l'aimant est calcule avec un fer de permeabilite
%  mu_I. Or les dents sont DEJA polarisees a ~1.5 T par les aimants : la
%  permeabilite INCREMENTALE offerte aux harmoniques d'induit y est bien
%  plus faible que la valeur a faible champ. On mesure l'effet.
clear; clc;
M=machine_bldc(); p=M.p; Ntc=M.Ntc; fea=M.FEA.dir; kfr=0.75;
PA=[1 -2 -15 3 14]; PB=[6 -7 -5 8 4]; PC=[11 -12 -10 13 9];
rd=@(f)readmatrix(f,'FileType','text','Delimiter','\t','NumHeaderLines',1);
ol=fullfile(fea,'transitoire (en charge)'); w4=@(x)x(round(numel(x)*0.75):end);
dLo=rd(fullfile(ol,'Plot 1_loss.tab')); PpmF=mean(w4(dLo(:,2)))/1e3;
R=cogging_mec(M,1260,0,721,M.muI,kfr,2*pi/p);
the=R.phis*p; lamf=@(P) Ntc*sum(sign(P(:)).*R.PhiT(abs(P(:)),:),1);
lam=[lamf(PA); lamf(PB); lamf(PC)];
RI=inductance_mec(M,1260,M.muI,0);
d4=rd(fullfile(fea,'transitoire (Back_emf)','Output Variables Plot 4.tab'));
dphi=-angle(2*mean(lam(1,:).*exp(-1i*the)) / ...
            (2*mean(d4(:,3).'.*exp(-1i*d4(:,2).'*p*pi/180))));
Q=load('mec_map.mat'); fn=fieldnames(Q); MAP=Q.(fn{1});
o=struct('dt',1e-6,'tend',34.26e-3,'J',1e-3,'Bf',6.07927101854027e-4, ...
    'Tload',4.8,'Vdc',500,'Rph',M.Rph,'dphi',dphi,'map',MAP);
D=drive_mec(M,the,lam,RI.Ld,o);
nM=mean(D.n(D.t>0.75*D.t(end))); Te=60/(nM*p);
mw=D.t>D.t(end)-Te; tw=D.t(mw); iw=D.i(:,mw); thw=D.th(mw);
js=round(linspace(1,numel(tw),961));
%  permeabilite incrementale reelle des dents, lue sur la courbe B(H) au
%  niveau de polarisation calcule par le MEC a vide
BH=bh_curve(); Btnl=max(abs(R.PhiT(:)))/(M.wst1*M.ls*M.Ki);
fprintf('=== Pertes aimant : sensibilite a mu_I du champ d''induit ===\n');
fprintf('  reference FEA : %.3f W\n',PpmF);
fprintf('  polarisation des dents par les aimants : B_t = %.2f T -> mu_r = %.0f\n', ...
    Btnl,BH.mur(Btnl));
fprintf('  %10s %10s %10s\n','mu_I','P (W)','ecart');
for mu=[3000 1000 400 150 50]
    RIm=inductance_mec(M,1260,mu,kfr);
    Pm=pm_loss_load(M,R,RIm.Usurf,tw(js),iw(:,js),thw(js),25);
    fprintf('  %10d %10.4f %8.1f %%\n',mu,Pm,100*(Pm-PpmF)/PpmF);
end
muE=BH.mur(Btnl);
RIe=inductance_mec(M,1260,muE,kfr);
Pe=pm_loss_load(M,R,RIe.Usurf,tw(js),iw(:,js),thw(js),25);
fprintf('\n  a la permeabilite INCREMENTALE physique (mu_r = %.0f) : %.4f W (%+.1f %%)\n', ...
    muE,Pe,100*(Pe-PpmF)/PpmF);
