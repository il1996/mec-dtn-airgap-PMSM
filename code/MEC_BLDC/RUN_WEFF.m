%% RUN_WEFF  -  Hypothese UNIQUE : ouverture d'encoche EFFECTIVE (bec sature)
%
%  Diagnostic etabli : le champ d'AIMANT est exact (ordres 7/21/35 a +-1.4 %),
%  mais la 1ere bande de denture est -28 % et la 2e +12 %. Geometrie, profondeur,
%  mu et saturation lumped sont ecartees. Trois indicateurs INDEPENDANTS
%  (amplitude a8, rapport 2e/1ere bande, niveau Bg1 via Carter) pointent tous
%  vers une OUVERTURE EFFECTIVE plus large que les 2 mm geometriques.
%
%  Mecanisme : les CORNES du bec de dent (vives dans le modele EF) saturent
%  localement -- l'induction de dent est deja 1.46 T pour un M350 qui sature a
%  1.83 T -- et le fer sature s'y comporte comme de l'air : le champ "voit"
%  une gorge elargie. C'est un effet 2D local qu'un reseau a dent lumped ne
%  peut pas produire.
%
%  Test : UN SEUL parametre (ws0_eff) doit corriger SIMULTANEMENT quatre
%  grandeurs independantes. Si oui, l'hypothese est validee ; sinon, refutee.
clear; clc;
M=machine_bldc();
D=readmatrix(fullfile(M.FEA.dir,'magnetostique(Magnetic_loading)', ...
    'Calculator Expressions Plot 4.tab'),'FileType','text','Delimiter','\t','NumHeaderLines',1);
thF=D(:,2)*pi/180; BrF=D(:,3); [thF,ix]=sort(thF); BrF=BrF(ix);
D2=readmatrix(fullfile(M.FEA.dir,'magnetostique(Magnetic_loading)', ...
    'Calculator Expressions Plot 2.tab'),'FileType','text','Delimiter','\t','NumHeaderLines',1);
fu=@(y,t,k)hypot(2*mean(y(:).*cos(k*t(:))),2*mean(y(:).*sin(k*t(:))));
F.Bg1=fu(BrF,thF,M.p); F.a8=fu(BrF,thF,8); F.a22=fu(BrF,thF,22);
F.a23=fu(BrF,thF,23); F.a37=fu(BrF,thF,37);
F.Bmean=mean(abs(BrF)); F.Bpeak=max(BrF); F.Btrms=sqrt(mean(D2(:,end).^2));

K=8;
fprintf('=== Balayage de l''ouverture EFFECTIVE (fer non lineaire M350) ===\n');
fprintf('geometrique ws0 = 2.000 mm ; bec sature -> ws0_eff > ws0\n\n');
fprintf('%8s | %8s %8s %8s %8s | %7s %8s\n', ...
    'ws0_eff','Bg1','a8','a22','a23','a22/a8','err4');
best=inf; wbest=NaN;
wl=2.0:0.1:3.2;
res=zeros(numel(wl),5);
for q=1:numel(wl)
    M2=M; M2.ws0=wl(q)*1e-3;
    nm=ceil(4*K*pi/(M2.ws0/M2.Rsi));
    R=subdomain_mec(M2,K,nm,'nl',0);
    % erreur composite sur 4 grandeurs INDEPENDANTES du spectre
    e=[ (R.Bg1-F.Bg1)/F.Bg1, (R.a8-F.a8)/F.a8, (R.a22-F.a22)/F.a22, (R.a23-F.a23)/F.a23 ];
    err=sqrt(mean(e.^2))*100;
    res(q,:)=[R.Bg1 R.a8 R.a22 R.a23 err];
    fprintf('%8.2f | %8.4f %8.5f %8.5f %8.5f | %7.3f %7.1f%%\n', ...
        wl(q),R.Bg1,R.a8,R.a22,R.a23,R.a22/R.a8,err);
    if err<best, best=err; wbest=wl(q); end
end
fprintf('%8s | %8.4f %8.5f %8.5f %8.5f | %7.3f\n','FEA',F.Bg1,F.a8,F.a22,F.a23,F.a22/F.a8);
fprintf('\n>> OPTIMUM : ws0_eff = %.2f mm  (erreur composite %.1f %%)\n',wbest,best);
fprintf('   elargissement = %.2f mm, soit %.2f mm par corne de bec\n', ...
    wbest-2.0,(wbest-2.0)/2);

% --- bilan complet a l'optimum ---
M2=M; M2.ws0=wbest*1e-3; nm=ceil(4*K*pi/(M2.ws0/M2.Rsi));
R=subdomain_mec(M2,K,nm,'nl',0);
fprintf('\n--- BILAN au ws0_eff optimal ---\n');
pr=@(n,a,b)fprintf('  %-24s %9.4f %9.4f %+8.1f%%\n',n,a,b,100*(a-b)/b);
fprintf('  %-24s %9s %9s %9s\n','grandeur','MEC','FEA','ecart');
pr('Bg1 (T)',R.Bg1,F.Bg1); pr('B moyen (T)',R.Bmean,F.Bmean);
pr('B crete (T)',R.Bpeak,F.Bpeak); pr('Bt RMS (T)',R.Btrms,F.Btrms);
pr('a8   |Ns-p| (T)',R.a8,F.a8); pr('a22   Ns+p (T)',R.a22,F.a22);
pr('a23 |2Ns-p| (T)',R.a23,F.a23); pr('a37  2Ns+p (T)',hypot(R.Brc(37,1),R.Brs(37,1)),F.a37);
fprintf('\n  (a8 et a22 sont les harmoniques qui pilotent les PERTES AIMANT)\n');
save('weff_opt.mat','wbest');
