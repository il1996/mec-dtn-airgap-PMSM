%% RUN_DIAG_SEQ - ordre des phases et calage de commutation
%  Le balayage de RUN_DIAG_LOAD montre un couple maximal a 90 deg du calage
%  etabli sur le flux de la SEULE phase A. Un calage correct doit tomber du
%  premier coup : on verifie donc l'ordre des phases (sens de rotation) en
%  confrontant les TROIS flux, puis on relit la position de commutation.
clear; clc;
M=machine_bldc(); p=M.p; Ntc=M.Ntc; fea=M.FEA.dir;
PA=[1 -2 -15 3 14]; PB=[6 -7 -5 8 4]; PC=[11 -12 -10 13 9];
if isfile('kfringe_ident.mat'), S0=load('kfringe_ident.mat'); kfr=S0.kbest; else, kfr=0.75; end
rd=@(f)readmatrix(f,'FileType','text','Delimiter','\t','NumHeaderLines',1);

R=cogging_mec(M,1260,0,721,M.muI,kfr,2*pi/p);
the=R.phis*p; lamf=@(P) Ntc*sum(sign(P(:)).*R.PhiT(abs(P(:)),:),1);
lam=[lamf(PA); lamf(PB); lamf(PC)];

d=rd(fullfile(fea,'transitoire (Back_emf)','Output Variables Plot 4.tab'));
pos=d(:,2)*p*pi/180;                     % angle ELECTRIQUE ANSYS (rad)
FL=[d(:,3) d(:,7) d(:,5)].';             % ordre du fichier : a, c, b -> a, b, c

fprintf('=== Phase du fondamental du flux totalise ===\n');
fprintf('  %-6s %14s %14s %12s\n','phase','MEC (deg)','FEA (deg)','MEC-FEA');
nm='abc'; dM=zeros(1,3); dF=dM;
for k=1:3
    cM=2*mean(lam(k,:).*exp(-1i*the));
    cF=2*mean(FL(k,:).*exp(-1i*pos.'));
    dM(k)=angle(cM)*180/pi; dF(k)=angle(cF)*180/pi;
    fprintf('  %-6s %14.2f %14.2f %12.2f\n',nm(k),dM(k),dF(k), ...
        mod(dM(k)-dF(k)+180,360)-180);
end
fprintf('\n  ecarts entre phases (doivent valoir -120 deg pour une sequence directe)\n');
fprintf('  MEC : b-a = %+7.1f   c-b = %+7.1f\n', ...
    mod(dM(2)-dM(1)+180,360)-180,mod(dM(3)-dM(2)+180,360)-180);
fprintf('  FEA : b-a = %+7.1f   c-b = %+7.1f\n', ...
    mod(dF(2)-dF(1)+180,360)-180,mod(dF(3)-dF(2)+180,360)-180);
if sign(mod(dM(2)-dM(1)+180,360)-180) ~= sign(mod(dF(2)-dF(1)+180,360)-180)
    fprintf('  => SEQUENCE INVERSEE entre le MEC et ANSYS.\n');
else
    fprintf('  => meme sequence de phases.\n');
end

%% ---- calage de commutation LU sur l'essai a vide ----------------------
%  La commande six-step d'ANSYS met A+ a ON pour theta_e dans [0,120).
%  Le couple est maximal si, sur cet intervalle, (psi_a - psi_b) est
%  maximal. On mesure directement OU se trouve ce maximum, en angle ANSYS.
psiF=zeros(3,numel(pos));
for k=1:3
    c=2*mean(FL(k,:).*exp(-1i*pos.'));   % fondamental
    psiF(k,:)=-p*abs(c)*sin(pos.'+angle(c));
end
dab=psiF(1,:)-psiF(2,:);
[~,j]=max(dab);
fprintf('\n  (psi_a - psi_b) maximal a theta_e = %.1f deg (convention ANSYS)\n', ...
    mod(pos(j)*180/pi,360));
fprintf('  fenetre de conduction A+/B- : theta_e dans [0,60) -> centre 30 deg\n');
fprintf('  => avance de commutation ANSYS : %+.1f deg elec\n', ...
    mod(30-mod(pos(j)*180/pi,360)+180,360)-180);
