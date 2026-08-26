%% RUN_DIAG_SIGNE  -  Pourquoi FEM et FLUX demandent-ils des decalages differents ?
%  e_a = -omega*dlambda_a/dtheta : les deux grandeurs viennent de LA MEME
%  solution. Si elles exigent des recalages differant de ~180 deg, c'est une
%  CONVENTION DE SIGNE, pas un decalage a ajuster.
clear; clc;
M=machine_bldc(); fea=M.FEA.dir; p=M.p; Ntc=M.Ntc;
PA=[1 -2 -15 3 14];
rd=@(f)readmatrix(f,'FileType','text','Delimiter','\t','NumHeaderLines',1);
if isfile('kfringe_ident.mat'), S=load('kfringe_ident.mat'); kfr=S.kbest; else, kfr=0.75; end
R=cogging_mec(M,1260,0,721,M.muI,kfr,2*pi/p); phis=R.phis; om=M.speed*2*pi/60;
lamA=Ntc*sum(sign(PA(:)).*R.PhiT(abs(PA(:)),:),1);
eA=-gradient(lamA,phis)*om;

dV=rd(fullfile(fea,'transitoire (Back_emf)','NodeVoltage Plot.tab'));
dFL=rd(fullfile(fea,'transitoire (Back_emf)','Output Variables Plot 4.tab'));
posV=dV(:,2)*p*pi/180; EaF=dV(:,3);
posFL=dFL(:,2)*p*pi/180; FLaF=dFL(:,3);

f1=@(y,t)2*mean(y(:).*exp(-1i*t(:)));
cLM=f1(lamA,phis*p); cLF=f1(FLaF,posFL);
cEM=f1(eA,phis*p);   cEF=f1(EaF,posV);

fprintf('=== Phases des fondamentaux (deg electriques) ===\n');
fprintf('  flux totalise : MEC %+8.2f   FEA %+8.2f  -> decalage %+8.2f\n', ...
    angle(cLM)*180/pi,angle(cLF)*180/pi,angle(cLM/cLF)*180/pi);
fprintf('  FEM de phase  : MEC %+8.2f   FEA %+8.2f  -> decalage %+8.2f\n', ...
    angle(cEM)*180/pi,angle(cEF)*180/pi,angle(cEM/cEF)*180/pi);
d=mod(angle(cEM/cEF)*180/pi - angle(cLM/cLF)*180/pi + 180,360)-180;
fprintf('\n  DIFFERENCE des deux decalages : %+.2f deg\n',d);
if abs(abs(d)-180)<25
    fprintf('  => ~180 deg : CONVENTION DE SIGNE de la FEM, pas un decalage.\n');
end

fprintf('\n=== Verification directe sur les donnees FEA ===\n');
dl=gradient(FLaF,posFL);
fprintf('  correlation( e_FEA , +dlambda_FEA/dtheta ) = %+.4f\n', ...
    corr(EaF,dl));
fprintf('  correlation( e_FEA , -dlambda_FEA/dtheta ) = %+.4f\n',corr(EaF,-dl));
dlM=gradient(lamA,phis*p);
fprintf('  correlation( e_MEC , -dlambda_MEC/dtheta ) = %+.4f  (ma convention)\n', ...
    corr(eA(:),-dlM(:)));
fprintf(['\n  => ANSYS (NodeVoltage Ie_a) fournit +dlambda/dt ; le MEC calcule\n' ...
    '     e = -omega*dlambda/dtheta. Les deux flux, eux, sont de MEME signe.\n' ...
    '     Il faut donc UN SEUL decalage (celui du flux) + le signe de la FEM.\n']);
