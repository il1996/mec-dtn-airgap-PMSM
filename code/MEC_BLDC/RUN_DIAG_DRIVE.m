%% RUN_DIAG_DRIVE - identification EXACTE du convertisseur de l'essai en charge
%  Objectif : etablir, SANS aucun ajustement, la loi de commutation utilisee
%  par ANSYS, pour que le modele MEC transitoire soit pilote a l'identique.
%  Source : netlist 'vf-15-14.ckt' extraite de BLDC.aedt
%     Vdc=500 (bus fractionne +/-250), Rph=10.22, etoile a neutre isole,
%     6 interrupteurs commandes par PULSE de periode 720/Nm et de largeur
%     720/(3*Nm) -> conduction 120 deg. Reste a savoir si l'argument des
%     PULSE est un TEMPS ou une POSITION.
clear; clc;
M=machine_bldc(); Nm=M.Nm; p=Nm/2;
fea=M.FEA.dir; ol=fullfile(fea,'transitoire (en charge)');
rd=@(f)readmatrix(f,'FileType','text','Delimiter','\t','NumHeaderLines',1);

%% ---- vitesse mesuree -> position mecanique ----
ds=rd(fullfile(ol,'Speed Plot 1.tab'));  ts=ds(:,1)*1e-3; ns=ds(:,2);
oms=ns*2*pi/60;
thm=cumtrapz(ts,oms)*180/pi;                       % position mecanique (deg)

%% ---- signaux de commande ----
dg=rd(fullfile(ol,'NodeVoltage Plot 1.tab'));      % temps en ns
tg=dg(:,1)*1e-9;  G=dg(:,2:7)>0;                   % Ivca1 Ivca2 Ivcb1 Ivcb2 Ivcc1 Ivcc2
nmg={'Ivca1(A+)','Ivca2(A-)','Ivcb1(B+)','Ivcb2(B-)','Ivcc1(C+)','Ivcc2(C-)'};
thg=interp1(ts,thm,tg,'linear','extrap');

fprintf('=== Loi de commutation de l''essai en charge ===\n');
fprintf('  duree simulee %.2f ms, %d pas ; vitesse finale %.0f tr/min\n', ...
    ts(end)*1e3,numel(ts),ns(end));
fprintf('  position finale %.0f deg mec (%.1f tours)\n',thm(end),thm(end)/360);
fprintf('  periode annoncee par la netlist : 720/Nm = %.4f\n',720/Nm);

fprintf('\n  %-11s | %-28s | %-28s\n','signal','fronts montants en TEMPS (ms)','... en POSITION (deg mec)');
for k=1:6
    j=find(diff(G(:,k))>0);
    dt_=diff(tg(j+1))*1e3; dth=diff(thg(j+1));
    fprintf('  %-11s | ecart %6.3f +/- %6.3f     | ecart %7.3f +/- %6.3f\n', ...
        nmg{k},mean(dt_),std(dt_),mean(dth),std(dth));
end

%% ---- verdict ----
j=find(diff(G(:,1))>0); dth=diff(thg(j+1)); dt_=diff(tg(j+1))*1e3;
fprintf('\n--- VERDICT ---\n');
fprintf('  dispersion relative en TEMPS    : %.1f %%\n',100*std(dt_)/mean(dt_));
fprintf('  dispersion relative en POSITION : %.1f %%\n',100*std(dth)/mean(dth));
if std(dth)/mean(dth) < std(dt_)/mean(dt_)
    fprintf('  => commutation pilotee par la POSITION (equivalent capteurs Hall).\n');
    fprintf('     periode mesuree %.3f deg mec  vs  720/Nm = %.3f deg  (ecart %.2f %%)\n', ...
        mean(dth),720/Nm,100*(mean(dth)-720/Nm)/(720/Nm));
else
    fprintf('  => commutation pilotee par le TEMPS (boucle ouverte V/f).\n');
end

%% ---- phase de commutation : position du 1er front de chaque signal ----
fprintf('\n  --- calage angulaire (deg ELECTRIQUES, modulo 360) ---\n');
for k=1:6
    j=find(diff(G(:,k))>0);
    if isempty(j), th1=thg(1); else, th1=thg(j(1)+1); end
    fprintf('  %-11s premier front a %8.3f deg mec  ->  %7.2f deg elec\n', ...
        nmg{k},th1,mod(th1*p,360));
end
fprintf('  netlist : A+ 0, C- 60, B+ 120, A- 180, C+ 240, B- 300 deg elec\n');

%% ---- sequence observee sur un pas : quelles phases conduisent ? ----
fprintf('\n  --- sequence des 6 etats sur une periode electrique ---\n');
i0=find(tg>tg(end)*0.8,1);
lab='ABC';
for k=i0:i0+11
    up=find(G(k,[1 3 5])); dn=find(G(k,[2 4 6]));
    su=''; if ~isempty(up), su=lab(up); end
    sd=''; if ~isempty(dn), sd=lab(dn); end
    fprintf('   t=%6.3f ms  th_e=%6.1f deg  +[%-3s] -[%-3s]\n', ...
        tg(k)*1e3,mod(thg(k)*p,360),su,sd);
end

%% ---- parametres mecaniques (BLDC.aedt, MotionSetup1) ----
fprintf('\n  --- parametres mecaniques lus dans BLDC.aedt ---\n');
fprintf('   J = 0.001 kg.m2 | frottement = 6.07927e-4 N.m.s | couple resistant = 4.8 N.m\n');
fprintf('   vitesse initiale 0 rpm ; regime : T = 4.8 + 6.079e-4*om\n');
omf=ns(end)*2*pi/60;
fprintf('   -> couple attendu en regime : %.4f N.m a %.0f tr/min\n',4.8+6.07927101854027e-4*omf,ns(end));
