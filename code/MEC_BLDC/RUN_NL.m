%% RUN_NL  -  Effet de la SATURATION LOCALE (mu par dent) sur la denture
%  Les dents reelles travaillent a mu_r ~ 1150 (B_fer 1.47 T) et la culasse a
%  ~6700 : un mu UNIFORME ne peut pas representer cette repartition. Ici on
%  resout avec mu(B) local PAR DENT et PAR SEGMENT DE CULASSE (courbe M350
%  reelle du projet ANSYS), sans aucun parametre ajuste.
clear; clc;
M=machine_bldc();
D=readmatrix(fullfile(M.FEA.dir,'magnetostique(Magnetic_loading)', ...
    'Calculator Expressions Plot 4.tab'),'FileType','text','Delimiter','\t','NumHeaderLines',1);
thF=D(:,2)*pi/180; BrF=D(:,3); [thF,ix]=sort(thF); BrF=BrF(ix);
D2=readmatrix(fullfile(M.FEA.dir,'magnetostique(Magnetic_loading)', ...
    'Calculator Expressions Plot 2.tab'),'FileType','text','Delimiter','\t','NumHeaderLines',1);
fu=@(y,t,k)hypot(2*mean(y(:).*cos(k*t(:))),2*mean(y(:).*sin(k*t(:))));
F.Bg1=fu(BrF,thF,M.p); F.a8=fu(BrF,thF,8); F.a22=fu(BrF,thF,22); F.a23=fu(BrF,thF,23);
F.Bmean=mean(abs(BrF)); F.Bpeak=max(BrF); F.Btrms=sqrt(mean(D2(:,end).^2));

K=8; nm=ceil(4*K*pi/(M.ws0/M.Rsi));
fprintf('=== Saturation locale : mu UNIFORME vs mu(B) PAR DENT ===\n');
fprintf('%-22s %8s %8s %8s %8s %8s\n','modele','Bg1','a8','a22','a23','a22/a8');
for mu={3000,1e7,'nl'}
    R=subdomain_mec(M,K,nm,mu{1},0);
    if ischar(mu{1}), lab='mu(B) local M350';
    else, lab=sprintf('mu uniforme = %g',mu{1}); end
    fprintf('%-22s %8.4f %8.5f %8.5f %8.5f %8.3f\n',lab,R.Bg1,R.a8,R.a22,R.a23,R.a22/R.a8);
end
fprintf('%-22s %8.4f %8.5f %8.5f %8.5f %8.3f\n','FEA (ANSYS)',F.Bg1,F.a8,F.a22,F.a23,F.a22/F.a8);

R=subdomain_mec(M,K,nm,'nl',0);
fprintf('\n--- Etat magnetique local (non lineaire) ---\n');
fprintf('mu_r des 15 dents : min %.0f  max %.0f  moy %.0f\n', ...
    min(R.muT(:,1)),max(R.muT(:,1)),mean(R.muT(:,1)));
fprintf('B_fer dents : min %.3f  max %.3f T\n', ...
    min(abs(R.PhiT(:,1)))/R.Abody,max(abs(R.PhiT(:,1)))/R.Abody);
fprintf('iterations point fixe : %d   rcond %.1e\n',R.iters(1),R.rcond);

fprintf('\n--- BILAN non lineaire vs FEA ---\n');
pr=@(n,a,b)fprintf('  %-24s %9.4f %9.4f %+8.1f%%\n',n,a,b,100*(a-b)/b);
fprintf('  %-24s %9s %9s %9s\n','grandeur','MEC','FEA','ecart');
pr('Bg1 (T)',R.Bg1,F.Bg1); pr('B moyen (T)',R.Bmean,F.Bmean);
pr('B crete (T)',R.Bpeak,F.Bpeak); pr('Bt RMS (T)',R.Btrms,F.Btrms);
pr('a8  |Ns-p| (T)',R.a8,F.a8); pr('a22  Ns+p (T)',R.a22,F.a22);
pr('a23 |2Ns-p| (T)',R.a23,F.a23);
