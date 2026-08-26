%% RUN_SUBDOMAIN  -  Validation du modele HYBRIDE sous-domaines <-> reseau
%  Objectif : supprimer le parametre ajuste kfringe et corriger le spectre de
%  denture (a8, a22, a23) qui pilote les pertes aimant.
clear; clc;
M=machine_bldc();

% --- reference FEA : profil d'entrefer magnetostatique ---
D=readmatrix(fullfile(M.FEA.dir,'magnetostique(Magnetic_loading)', ...
    'Calculator Expressions Plot 4.tab'),'FileType','text','Delimiter','\t','NumHeaderLines',1);
thF=D(:,2)*pi/180; BrF=D(:,3); [thF,ix]=sort(thF); BrF=BrF(ix);
fu=@(y,t,k)hypot(2*mean(y(:).*cos(k*t(:))), 2*mean(y(:).*sin(k*t(:))));
F.Bg1=fu(BrF,thF,M.p); F.a8=fu(BrF,thF,8); F.a22=fu(BrF,thF,22); F.a23=fu(BrF,thF,23);
F.Bmean=mean(abs(BrF)); F.Bpeak=max(BrF);
D2=readmatrix(fullfile(M.FEA.dir,'magnetostique(Magnetic_loading)', ...
    'Calculator Expressions Plot 2.tab'),'FileType','text','Delimiter','\t','NumHeaderLines',1);
F.Btrms=sqrt(mean(D2(:,end).^2));

fprintf('=== MODELE HYBRIDE SOUS-DOMAINES (aucun parametre ajuste) ===\n');
fprintf('ouverture delta = %.4f rad = %.2f deg ; lam_1 = %.1f\n\n', ...
    M.ws0/M.Rsi,(M.ws0/M.Rsi)*180/pi,pi/(M.ws0/M.Rsi));

% --- convergence en nombre de modes K ---
fprintf('%4s %7s | %8s %8s %8s %8s | %9s\n','K','numax','Bg1','a8','a22','a23','balance');
Klist=[2 4 6 8 10];
for K=Klist
    nm=ceil(3.2*K*pi/(M.ws0/M.Rsi));
    R=subdomain_mec(M,K,nm,[],0);
    fprintf('%4d %7d | %8.4f %8.5f %8.5f %8.5f | %9.1e\n', ...
        K,nm,R.Bg1,R.a8,R.a22,R.a23,R.balance);
end
fprintf('%4s %7s | %8.4f %8.5f %8.5f %8.5f |\n','FEA','-',F.Bg1,F.a8,F.a22,F.a23);

% --- convergence en numax a K fixe ---
K=8;
fprintf('\n%4s %7s | %8s %8s %8s\n','K','numax','a8','a22','a23');
for fac=[2 3 4 5]
    nm=ceil(fac*K*pi/(M.ws0/M.Rsi));
    R=subdomain_mec(M,K,nm,[],0);
    fprintf('%4d %7d | %8.5f %8.5f %8.5f\n',K,nm,R.a8,R.a22,R.a23);
end

% --- bilan final ---
K=8; nm=ceil(4*K*pi/(M.ws0/M.Rsi));
R=subdomain_mec(M,K,nm,[],0);
fprintf('\n--- BILAN (K=%d, numax=%d) ---\n',K,nm);
pr=@(n,a,b)fprintf('  %-26s %9.4f %9.4f %+8.1f%%\n',n,a,b,100*(a-b)/b);
fprintf('  %-26s %9s %9s %9s\n','grandeur','MEC','FEA','ecart');
pr('Bg1 fondamental (T)',R.Bg1,F.Bg1);
pr('B entrefer moyen (T)',R.Bmean,F.Bmean);
pr('B entrefer crete (T)',R.Bpeak,F.Bpeak);
pr('Bt tangentiel RMS (T)',R.Btrms,F.Btrms);
pr('harmonique 8  |Ns-p| (T)',R.a8,F.a8);
pr('harmonique 22  Ns+p (T)',R.a22,F.a22);
pr('harmonique 23 |2Ns-p|(T)',R.a23,F.a23);
fprintf('  rapport a22/a8 : MEC %.3f   FEA %.3f\n',R.a22/R.a8,F.a22/F.a8);
fprintf('  conservation de flux dans les ouvertures : %.1e\n',R.balance);
