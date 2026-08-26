%% RUN_IDENT_H8 - re-identification de kfringe sur les BANDES DE DENTURE
%
%  CE QUI CHANGE ET POURQUOI. kfringe est la conductance de frange de la
%  bouche d'encoche : c'est le seul parametre du reseau qui gouverne la
%  MODULATION D'ENCOCHE, donc les bandes laterales nu = |Ns -+ p| = 8 et 22.
%  Il avait ete identifie a 0.75 en minimisant l'ecart RMS sur la forme
%  d'onde complete du champ d'entrefer. Ce critere est mal pose pour un
%  reseau a branches localisees : l'ecart RMS est domine par la forme
%  POINTUE des encoches, qu'aucun modele a branches ne reproduit point par
%  point, et il se laisse minimiser en APLATISSANT la modulation. Resultat :
%  RMS optimal, mais rang 8 a -40 % -- et comme ce rang porte a lui seul la
%  denture vue par le rotor, il entraine avec lui la dispersion de FMM
%  (-29.5 %) et les DEUX pertes aimant (-50 %, en B^2).
%
%  On identifie donc kfringe sur les grandeurs qu'il represente : les deux
%  bandes laterales de denture. C'est un critere PHYSIQUE et non global, et
%  il ne comporte toujours qu'UN seul parametre.
clear; clc;
M=machine_bldc(); Ns=M.Ns; p=M.p; fea=M.FEA.dir;
rd=@(f)readmatrix(f,'FileType','text','Delimiter','\t','NumHeaderLines',1);
d4=rd(fullfile(fea,'magnetostique(Magnetic_loading)','Calculator Expressions Plot 4.tab'));
thu=linspace(0,2*pi,3601); thu(end)=[];
BrFu=interp1([d4(:,2)*pi/180;2*pi],[d4(:,3);d4(1,3)],thu,'linear','extrap');
amp=@(B,th,k) 2*abs(mean(B(:).'.*exp(-1i*k*th(:).')));
t8=amp(BrFu,thu,8); t22=amp(BrFu,thu,22); t1=amp(BrFu,thu,p);

kv=0.20:0.025:0.60; J=zeros(size(kv)); A1=J; A8=J; A22=J;
for q=1:numel(kv)
    R=cogging_mec(M,1260,0,3,M.muI,kv(q),1e-6);
    th=R.thq(1:end-1); Br=R.Br(1:end-1);
    A1(q)=amp(Br,th,p); A8(q)=amp(Br,th,8); A22(q)=amp(Br,th,22);
    J(q)=sqrt(((A8(q)-t8)/t8)^2+((A22(q)-t22)/t22)^2)/sqrt(2);
end
[Jb,qb]=min(J); kbest=kv(qb);
fprintf('=== Identification de kfringe sur les bandes de denture ===\n');
fprintf('  cibles FEA : rang 8 = %.5f T, rang 22 = %.5f T\n',t8,t22);
fprintf('  %8s %10s %9s %10s %9s %10s\n','kfringe','rang 8','ecart','rang 22','ecart','critere');
for q=1:2:numel(kv)
    mk=''; if q==qb, mk='  <--'; end
    fprintf('  %8.3f %10.5f %8.1f %% %10.5f %8.1f %% %9.1f %%%s\n', ...
        kv(q),A8(q),100*(A8(q)-t8)/t8,A22(q),100*(A22(q)-t22)/t22,100*J(q),mk);
end
fprintf('\n  kfringe identifie = %.3f (critere %.1f %%)\n',kbest,100*Jb);
fprintf('  fondamental Bg1 : %.5f T contre %.5f T FEA (%+.1f %%)\n', ...
    A1(qb),t1,100*(A1(qb)-t1)/t1);
fprintf('  ancienne valeur 0.750 : rang 8 -40.0 %%, rang 22 -25.2 %%, Bg1 +0.8 %%\n');
fprintf('\n  gain attendu sur les pertes aimant (elles varient en B^2) :\n');
r=(A8(qb)/0.01180)^2;
fprintf('    facteur %.2f -> a vide %.4f -> %.4f W (cible %.4f)\n',r,0.1679,0.1679*r,0.3335);
fprintf('                  -> en charge %.3f -> %.3f W (cible %.3f)\n',1.529,1.529*r,3.183);

if isfile('kfringe_ident.mat')
    S=load('kfringe_ident.mat');
    save('kfringe_ident_RMS.mat','-struct','S');   % on conserve l'ancienne
    fprintf('\n  ancienne identification (critere RMS, kbest = %.3f) sauvegardee\n',S.kbest);
    fprintf('  dans kfringe_ident_RMS.mat\n');
end
kbest_old=0.75; critere='bandes de denture nu = |Ns -+ p| = 8 et 22'; %#ok<NASGU>
save('kfringe_ident.mat','kbest','kv','J','A1','A8','A22','t8','t22','critere','kbest_old');
fprintf('  -> kfringe_ident.mat mis a jour (kbest = %.3f)\n',kbest);
