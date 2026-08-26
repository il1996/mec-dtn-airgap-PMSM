%% RUN_ARMSAT - effet de la saturation TOURNANTE sur les pertes aimant
clear; clc; t0=tic;
M=machine_bldc(); p=M.p; Ntc=M.Ntc; fea=M.FEA.dir;
PA=[1 -2 -15 3 14]; PB=[6 -7 -5 8 4]; PC=[11 -12 -10 13 9];
if isfile('kfringe_ident.mat'), S0=load('kfringe_ident.mat'); kfr=S0.kbest; else, kfr=0.325; end
rd=@(f)readmatrix(f,'FileType','text','Delimiter','\t','NumHeaderLines',1);
ol=fullfile(fea,'transitoire (en charge)'); w4=@(x)x(round(numel(x)*0.75):end);
dLo=rd(fullfile(ol,'Plot 1_loss.tab')); PpmF=mean(w4(dLo(:,2)))/1e3;

R=cogging_mec(M,1260,0,721,M.muI,kfr,2*pi/p);
the=R.phis*p; lamf=@(P) Ntc*sum(sign(P(:)).*R.PhiT(abs(P(:)),:),1);
lam=[lamf(PA); lamf(PB); lamf(PC)];
RI=inductance_mec(M,1260,M.muI,0);
RIu=inductance_mec(M,1260,M.muI,kfr);
d4=rd(fullfile(fea,'transitoire (Back_emf)','Output Variables Plot 4.tab'));
dphi=-angle(2*mean(lam(1,:).*exp(-1i*the)) / ...
            (2*mean(d4(:,3).'.*exp(-1i*d4(:,2).'*p*pi/180))));
Q=load('mec_map.mat'); fn=fieldnames(Q); MAP=Q.(fn{1});
o=struct('dt',1e-6,'tend',34.26e-3,'J',1e-3,'Bf',6.07927101854027e-4, ...
    'Tload',4.8,'Vdc',500,'Rph',M.Rph,'dphi',dphi,'map',MAP);
D=drive_mec(M,the,lam,RI.Ld,o);
nM=mean(D.n(D.t>0.75*D.t(end))); Te=60/(nM*p);
Rfull=cogging_mec(M,1260,0,2161,M.muI,kfr,2*pi);
mw=D.t>D.t(end)-Te; tw=D.t(mw); iw=D.i(:,mw); thw=D.th(mw);
js=round(linspace(1,numel(tw),961));
fprintf('=== Pertes aimant en charge : saturation tournante ===\n');
fprintf('  reference FEA : %.3f W\n',PpmF);

%% ---- champ d'induit module par la saturation tournante ---------------
if isfile('arm_sat.mat')
    Qa=load('arm_sat.mat'); AF=Qa.AF;
else
    fprintf('  calcul du champ d''induit sature ... ');
    AF=arm_field_sat(M,1260,kfr,361); save('arm_sat.mat','AF','-v7.3');
    fprintf('%.0f s\n',toc(t0));
end
fprintf('  permeabilite de dent sur un tour : %.0f a %.0f (rapport %.1f)\n', ...
    AF.mu_min,AF.mu_max,AF.mu_max/AF.mu_min);

P0=pm_loss_load(M,Rfull,RIu.Usurf,tw(js),iw(:,js),thw(js),41);
P1=pm_loss_load(M,Rfull,AF,          tw(js),iw(:,js),thw(js),41);
Pd=pm_loss_load(M,Rfull,RIu.Usurf,tw(js),iw(:,js)*0,thw(js),41);
fprintf('\n  %-42s %9s %9s\n','','P (W)','ecart');
fprintf('  %-42s %9.4f %8.1f %%\n','induit LINEAIRE (etat precedent)',P0,100*(P0-PpmF)/PpmF);
fprintf('  %-42s %9.4f %8.1f %%\n','induit MODULE par la saturation',P1,100*(P1-PpmF)/PpmF);
fprintf('  %-42s %9.4f\n','denture seule (courants nuls)',Pd);
fprintf('\n  contribution d''induit : %.3f -> %.3f W (x%.2f)\n',P0-Pd,P1-Pd,(P1-Pd)/(P0-Pd));
fprintf('  facteur d''amplitude gagne sur le champ : %.2f\n',sqrt((P1-Pd)/(P0-Pd)));
fprintf('\n  (%.0f s)\n',toc(t0));
