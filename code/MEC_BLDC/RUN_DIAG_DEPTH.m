%% RUN_DIAG_DEPTH  -  Le "trou magnetique" de l'encoche est-il trop peu profond ?
%  Le sous-domaine d'ouverture est ferme par Neumann a r = Rs+hs0 (1 mm).
%  Or l'encoche reelle s'evase (ws0=2 -> ws1=6.7 mm) et se prolonge sur 21 mm :
%  le "trou" vu par le champ est bien plus profond. On balaie ici la
%  profondeur (et la largeur) de la cavite, a fer INFINIMENT permeable, pour
%  isoler l'effet geometrique de tout effet de permeabilite.
clear; clc;
M=machine_bldc();
D=readmatrix(fullfile(M.FEA.dir,'magnetostique(Magnetic_loading)', ...
    'Calculator Expressions Plot 4.tab'),'FileType','text','Delimiter','\t','NumHeaderLines',1);
thF=D(:,2)*pi/180; BrF=D(:,3); [thF,ix]=sort(thF); BrF=BrF(ix);
fu=@(y,t,k)hypot(2*mean(y(:).*cos(k*t(:))),2*mean(y(:).*sin(k*t(:))));
Fa8=fu(BrF,thF,8); Fa22=fu(BrF,thF,22); Fa23=fu(BrF,thF,23); FBg1=fu(BrF,thF,M.p);
mu=1e7; K=8;

fprintf('=== Effet de la PROFONDEUR de la cavite d''encoche (fer infini) ===\n');
fprintf('largeur ws0 = %.2f mm constante\n',M.ws0*1e3);
fprintf('%9s | %8s %8s %8s %8s | %7s\n','hs0 [mm]','Bg1','a8','a22','a23','a22/a8');
for hd=[0.5 1 2 4 8 16]
    M2=M; M2.hs0=hd*1e-3;
    nm=ceil(4*K*pi/(M2.ws0/M2.Rsi));
    R=subdomain_mec(M2,K,nm,mu,0);
    fprintf('%9.1f | %8.4f %8.5f %8.5f %8.5f | %7.3f\n',hd,R.Bg1,R.a8,R.a22,R.a23,R.a22/R.a8);
end
fprintf('%9s | %8.4f %8.5f %8.5f %8.5f | %7.3f\n','FEA',FBg1,Fa8,Fa22,Fa23,Fa22/Fa8);

fprintf('\n=== Effet de la LARGEUR de la cavite (profondeur 4 mm) ===\n');
fprintf('%9s | %8s %8s %8s %8s | %7s\n','ws0 [mm]','Bg1','a8','a22','a23','a22/a8');
for wd=[2 3 4 5 6.7]
    M2=M; M2.ws0=wd*1e-3; M2.hs0=4e-3;
    nm=ceil(4*K*pi/(M2.ws0/M2.Rsi));
    R=subdomain_mec(M2,K,nm,mu,0);
    fprintf('%9.2f | %8.4f %8.5f %8.5f %8.5f | %7.3f\n',wd,R.Bg1,R.a8,R.a22,R.a23,R.a22/R.a8);
end
fprintf('%9s | %8.4f %8.5f %8.5f %8.5f | %7.3f\n','FEA',FBg1,Fa8,Fa22,Fa23,Fa22/Fa8);
fprintf('\n(la vraie encoche : gorge 2 mm sur 1 mm, PUIS evasement a 6.7 mm\n');
fprintf(' et corps jusqu''a 21.1 mm -> le "trou" effectif est entre les deux)\n');
