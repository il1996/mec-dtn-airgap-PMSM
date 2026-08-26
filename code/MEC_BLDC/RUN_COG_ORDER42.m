%% RUN_COG_ORDER42  -  Les 220 mN.m a l'ordre 42 sont-ils REELS ?
%
%  TEST DECISIF. Les deux essais FEA portent sur LA MEME machine :
%    - magnetostatique : 360 positions, REMAILLEE a chaque position ;
%    - transitoire     : bande glissante, maillage PRESERVE.
%  Si le couple de 220 mN.m a l'ordre 42 etait un couple PHYSIQUE, il
%  apparaitrait dans LES DEUX. L'ordre 42 est parfaitement resolu par
%  l'essai transitoire (fenetre 1/7 de tour -> ordres multiples de 7,
%  42 = 6x7, Nyquist ~350). Son absence dans le transitoire prouverait
%  que c'est un artefact du REMAILLAGE.
clear; clc;
M=machine_bldc(); Ns=M.Ns; Nm=M.Nm; LCMv=lcm(Ns,Nm);
rd=@(f)readmatrix(f,'FileType','text','Delimiter','\t','NumHeaderLines',1);

% ---- (1) magnetostatique ----
dm=rd(fullfile(M.FEA.dir,'magnetostique(Magnetic_loading)','Torque Plot 1.tab'));
am=dm(:,1); Tm=dm(:,2); if abs(am(end)-360)<1e-6, am(end)=[]; Tm(end)=[]; end
Tm=Tm-mean(Tm); thm=am*pi/180;
% ---- (2) transitoire ----
dt=rd(fullfile(M.FEA.dir,'transitoire (Back_emf)','Torque Plot.tab'));
tht=dt(:,1)*1e-3*(M.speed/60)*2*pi; Tt=dt(:,2); Tt=Tt-mean(Tt);

fprintf('=== Les 220 mN.m a l''ordre 42 sont-ils physiques ? ===\n');
fprintf('  magnetostatique : %d positions sur %.0f deg (remaillee)\n',numel(thm),max(am));
fprintf('  transitoire     : %d positions sur %.2f deg (maillage preserve)\n', ...
    numel(tht),max(tht)*180/pi);
fprintf('  fenetre du transitoire = 1/%.0f de tour -> ordres multiples de %d\n', ...
    2*pi/max(tht),round(2*pi/max(tht)));

% ---- amplitude d'un ordre donne par moindres carres, avec incertitude ----
    function [A,s]=fitord(th,T,n)
        X=[cos(n*th) sin(n*th) ones(size(th)) th th.^2];
        c=X\T; r=T-X*c; A=hypot(c(1),c(2)); s=std(r)*sqrt(2/numel(th));
    end

ords=[42 84 126 210 150];
fprintf('\n  %6s | %22s | %22s\n','ordre','magnetostatique','transitoire');
fprintf('  %6s | %10s %11s | %10s %11s\n','','ampl.','incert.','ampl.','incert.');
for n=ords
    [A1,s1]=fitord(thm,Tm,n);
    [A2,s2]=fitord(tht,Tt,n);
    fprintf('  %6d | %10.3f +/-%7.3f | %10.3f +/-%7.3f\n',n,A1,s1,A2,s2);
end

[A42m,s42m]=fitord(thm,Tm,42);
[A42t,s42t]=fitord(tht,Tt,42);
fprintf('\n--- VERDICT ---\n');
fprintf('  ordre 42 dans la magnetostatique : %.2f mN.m  (%.0f sigma)\n',A42m,A42m/s42m);
fprintf('  ordre 42 dans le transitoire     : %.2f mN.m  (%.1f sigma)\n',A42t,A42t/s42t);
fprintf('  rapport magnetostatique/transitoire : %.0f x\n',A42m/max(A42t,eps));
if A42t < 3*s42t
    fprintf('  => L''ordre 42 est ABSENT du transitoire (sous le bruit).\n');
    fprintf('     Meme machine, meme physique, seule difference : le REMAILLAGE.\n');
    fprintf('     Les 220 mN.m sont donc un ARTEFACT NUMERIQUE, non un couple.\n');
    fprintf('     Il n''existe aucun modele physique qui puisse les reproduire.\n');
else
    fprintf('  => L''ordre 42 est present dans les DEUX : a re-examiner.\n');
end

% ---- controle de coherence : les deux essais voient-ils la meme detente ? ----
[A210m,~]=fitord(thm,Tm,LCMv); [A210t,s210t]=fitord(tht,Tt,LCMv);
fprintf('\n  detente (ordre %d) : magnetostatique %.3f (hors bande, Nyquist 180)\n', ...
    LCMv,A210m);
fprintf('                        transitoire     %.3f +/- %.3f mN.m  <-- reference\n', ...
    A210t,s210t);
