%% RUN_PMLOSS_SPLIT - d'ou viennent les -50 % sur les pertes aimant en charge ?
%  A vide la chaine est validee (+2.4 %) : la formulation des courants
%  induits et le DtN etendu sont donc justes. Le deficit est entierement
%  dans la contribution d'INDUIT (decomposition MEC : 0.18 W denture +
%  1.41 W induit). Deux causes possibles, qu'on separe ici :
%    (1) le COURANT predit par le MEC est 6 % trop faible (P varie en i^2) ;
%    (2) le CHAMP d'induit par ampere vu par l'aimant est sous-estime.
%  On recalcule donc la perte en injectant les courants MESURES par la FEA,
%  toutes choses egales par ailleurs. Ce qui reste apres est imputable au
%  seul champ.
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

%% ---- courants FEA ramenes sur la meme grille --------------------------
dI=rd(fullfile(ol,'BranchCurrent Plot 1.tab'));
dS=rd(fullfile(ol,'Speed Plot 1.tab'));
tF=dI(:,1)*1e-3; iF=dI(:,2:4).';
omF=interp1(dS(:,1)*1e-3,dS(:,2),tF,'linear','extrap')*2*pi/60;
thF=cumtrapz(tF,omF)*p;                              % angle elec FEA
%  on rejoue les courants FEA en fonction de l'ANGLE (pas du temps), pour
%  qu'ils tombent au bon endroit du motif de denture
theM=mod(thw(js)+dphi,2*pi);
thFm=mod(thF+dphi,2*pi); mF=tF>0.75*tF(end);
[u1,ku]=unique(thFm(mF)); iFu=iF(:,mF); iFu=iFu(:,ku);
iFq=interp1(u1,iFu.',mod(theM,2*pi),'linear','extrap').';

rmsM=sqrt(mean(iw(1,js).^2)); rmsF=sqrt(mean(iFq(1,:).^2));
fprintf('=== Pertes aimant en charge : courant ou champ ? ===\n');
fprintf('  reference FEA : %.3f W\n',PpmF);
fprintf('  courant de phase rms : MEC %.4f A | FEA rejoue %.4f A (%+.1f %%)\n', ...
    rmsM,rmsF,100*(rmsM-rmsF)/rmsF);

P_mm=pm_loss_load(M,Rfull,RIu.Usurf,tw(js),iw(:,js),thw(js),41);
P_fm=pm_loss_load(M,Rfull,RIu.Usurf,tw(js),iFq,thw(js),41);
P_dn=pm_loss_load(M,Rfull,RIu.Usurf,tw(js),iw(:,js)*0,thw(js),41);
fprintf('\n  %-38s %9s %9s\n','','P (W)','ecart');
fprintf('  %-38s %9.4f %8.1f %%\n','courants MEC (etat actuel)',P_mm,100*(P_mm-PpmF)/PpmF);
fprintf('  %-38s %9.4f %8.1f %%\n','courants FEA, meme champ',P_fm,100*(P_fm-PpmF)/PpmF);
fprintf('  %-38s %9.4f\n','denture seule (courants nuls)',P_dn);
fprintf('\n  part imputable au COURANT : %.3f W (%.0f %% du deficit)\n', ...
    P_fm-P_mm,100*(P_fm-P_mm)/(PpmF-P_mm));
fprintf('  part imputable au CHAMP   : %.3f W (%.0f %% du deficit)\n', ...
    PpmF-P_fm,100*(PpmF-P_fm)/(PpmF-P_mm));
fprintf('  facteur d''amplitude manquant sur le champ d''induit : %.2f\n', ...
    sqrt((PpmF-P_dn)/max(P_fm-P_dn,eps)));
fprintf('\n  (%.0f s)\n',toc(t0));
