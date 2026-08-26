%% RUN_DIAG_CYCLE - lecture interne du modele d'entrainement en regime
clear; clc;
M=machine_bldc(); p=M.p; Ntc=M.Ntc; fea=M.FEA.dir;
PA=[1 -2 -15 3 14]; PB=[6 -7 -5 8 4]; PC=[11 -12 -10 13 9];
if isfile('kfringe_ident.mat'), S0=load('kfringe_ident.mat'); kfr=S0.kbest; else, kfr=0.75; end
rd=@(f)readmatrix(f,'FileType','text','Delimiter','\t','NumHeaderLines',1);
R=cogging_mec(M,1260,0,721,M.muI,kfr,2*pi/p);
the=R.phis*p; lamf=@(P) Ntc*sum(sign(P(:)).*R.PhiT(abs(P(:)),:),1);
lam=[lamf(PA); lamf(PB); lamf(PC)];
RI=inductance_mec(M,1260,5000,0);
d=rd(fullfile(fea,'transitoire (Back_emf)','Output Variables Plot 4.tab'));
cLM=2*mean(lam(1,:).*exp(-1i*the));
cLF=2*mean(d(:,3).'.*exp(-1i*d(:,2).'*p*pi/180));
dphi=-angle(cLM/cLF);
Q=load('mec_map.mat');
om=1339.8*2*pi/60;
o=struct('dt',1e-6,'tend',20e-3,'J',1e9,'Bf',0,'Tload',0,'Vdc',500, ...
    'Rph',M.Rph,'dphi',dphi,'om0',om,'map',Q.S);
D=drive_mec(M,the,lam,RI.Ld,o);
n0=numel(D.t); k=n0-6400:150:n0;
fprintf('  t(ms) the_e     ia     ib     ic     e_ab     v_ab  Leq(mH)      T\n');
for j=k
    fprintf('%7.2f %6.1f %6.2f %6.2f %6.2f %8.1f %8.0f %8.1f %7.3f\n', ...
        D.t(j)*1e3,D.the(j)*180/pi,D.i(1,j),D.i(2,j),D.i(3,j), ...
        D.e(1,j)-D.e(2,j),D.v(1,j)-D.v(2,j),D.Leq(j)*1e3,D.T(j));
end
m=D.t>0.5*D.t(end);
fprintf('\n  moyennes : T = %.3f N.m | i rms = %.3f A | max|e_ab| = %.0f V\n', ...
    mean(D.T(m)),sqrt(mean(D.i(1,m).^2)),max(abs(D.e(1,m)-D.e(2,m))));
fprintf('  FEA      : T = 4.867 N.m | i rms = 1.552 A | e_ab mesure = 456 V\n');
ea=D.e(1,m); fprintf('  max|e_a| MEC = %.0f V   (FEA Ie_a crete = 233 V)\n',max(abs(ea)));
