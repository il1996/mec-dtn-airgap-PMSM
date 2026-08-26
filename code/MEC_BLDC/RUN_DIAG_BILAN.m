%% RUN_DIAG_BILAN - l'essai FEA en charge conserve-t-il l'energie ?
%  Avant de comparer un transitoire MEC a un transitoire FEA, il faut savoir
%  si la reference est exploitable. Test independant de tout modele : le bus
%  fournit une energie, elle doit se retrouver en pertes + energie cinetique
%  + travail sur la charge, a l'energie magnetique pres (nulle en moyenne).
clear; clc;
M=machine_bldc(); ol=fullfile(M.FEA.dir,'transitoire (en charge)');
rd=@(f)readmatrix(f,'FileType','text','Delimiter','\t','NumHeaderLines',1);
J=1e-3; Bf=6.07927101854027e-4; Tl=4.8; Vdc=500;

dL =rd(fullfile(ol,'Plot 1_loss.tab'));      % t[ns] Solid Core p_cu P_out Pro Pem
dIt=rd(fullfile(ol,'BranchCurrent Plot 2.tab'));
dS =rd(fullfile(ol,'Speed Plot 1.tab'));
dI =rd(fullfile(ol,'BranchCurrent Plot 1.tab'));
t=dL(:,1)*1e-9;                               % s
om=interp1(dS(:,1)*1e-3,dS(:,2),t,'linear','extrap')*2*pi/60;
Ib=interp1(dIt(:,1)*1e-3,dIt(:,2),t,'linear','extrap');
Pin=Vdc*Ib; Pfe=dL(:,3)/1e3; Ppm=dL(:,2)/1e3; Pem=dL(:,7);
%  la colonne p_cu du fichier a sa premiere valeur vide -> on la recalcule
%  avec la definition exacte d'ANSYS : p_cu = Rph*(ia^2+ib^2+ic^2)
iph=interp1(dI(:,1)*1e-3,dI(:,2:4),t,'linear','extrap');
Pcu=M.Rph*sum(iph.^2,2);
fprintf('  (p_cu recalcule : ecart max %.2f W avec la colonne du fichier)\n', ...
    max(abs(Pcu(2:end)-dL(2:end,4))));

fprintf('=== Bilan energetique de l''essai FEA en charge (0 -> %.2f ms) ===\n',t(end)*1e3);
Ein=trapz(t,Pin); Ecu=trapz(t,Pcu); Efe=trapz(t,Pfe); Epm=trapz(t,Ppm);
Ekin=0.5*J*om(end)^2; Efr=trapz(t,Bf*om.^2); Elo=trapz(t,Tl*om);
Eem=trapz(t,Pem);
fprintf('  energie fournie par le bus        %9.3f J\n',Ein);
fprintf('  pertes Joule                      %9.3f J\n',Ecu);
fprintf('  pertes fer + aimant               %9.3f J\n',Efe+Epm);
fprintf('  -> disponible pour la conversion  %9.3f J\n',Ein-Ecu-Efe-Epm);
fprintf('  energie cinetique acquise         %9.3f J\n',Ekin);
fprintf('  travail sur la charge (4.8 N.m)   %9.3f J\n',Elo);
fprintf('  travail contre le frottement      %9.3f J\n',Efr);
fprintf('  -> total mecanique                %9.3f J\n',Ekin+Elo+Efr);
fprintf('  integrale de Pem                  %9.3f J\n',Eem);
fprintf('\n  ECART conversion - mecanique      %+9.3f J  (%+.1f %%)\n', ...
    (Ein-Ecu-Efe-Epm)-(Ekin+Elo+Efr),100*((Ein-Ecu-Efe-Epm)-(Ekin+Elo+Efr))/(Ekin+Elo+Efr));
fprintf('  ECART Pem - mecanique             %+9.3f J  (%+.1f %%)\n', ...
    Eem-(Ekin+Elo+Efr),100*(Eem-(Ekin+Elo+Efr))/(Ekin+Elo+Efr));

%% ---- ou l'ecart se cree-t-il ? bilan instantane ------------------------
fprintf('\n=== Ou l''incoherence apparait-elle ? ===\n');
Pdisp=Pin-Pcu-Pfe-Ppm;                        % puissance convertible
fprintf('  %8s %10s %10s %10s %10s %10s\n', ...
    't (ms)','n (tr/min)','P_bus','P_cu','dispo','Pem');
for k=round(linspace(2,numel(t),16))
    fprintf('  %8.2f %10.0f %10.0f %10.0f %10.0f %10.0f\n', ...
        t(k)*1e3,om(k)*60/(2*pi),Pin(k),Pcu(k),Pdisp(k),Pem(k));
end
m1=t<12e-3; m2=t>25e-3;
fprintf('\n  demarrage (t < 12 ms) : Pem / dispo = %.2f   (doit valoir ~1)\n', ...
    trapz(t(m1),Pem(m1))/trapz(t(m1),Pdisp(m1)));
fprintf('  regime    (t > 25 ms) : Pem / dispo = %.2f\n', ...
    trapz(t(m2),Pem(m2))/trapz(t(m2),Pdisp(m2)));

%% ---- controle : la loi des mailles est-elle verifiee ? ----------------
dV=rd(fullfile(ol,'NodeVoltage Plot 2.tab'));
tv=dV(:,1)*1e-9; iv=interp1(dI(:,1)*1e-3,dI(:,2:4),tv,'linear','extrap');
Ie=dV(:,2:4)/1e3; Iv=dV(:,5:7)/1e3;
res=Iv-Ie-M.Rph*iv;                           % doit etre nul (loi d'Ohm sur Rph)
fprintf('\n  loi d''Ohm sur Rph : residu max %.3f V sur %.0f V  -> circuit coherent\n', ...
    max(abs(res(:))),max(abs(Iv(:))));
%  puissance entrant reellement dans les enroulements
Pw=sum(Ie.*iv,2);
fprintf('  puissance entrant dans les enroulements (somme e_k i_k) :\n');
fprintf('     au demarrage (t=5 ms)  %8.0f W   alors que Pem = %.0f W\n', ...
    interp1(tv,Pw,5e-3),interp1(t,Pem,5e-3));
fprintf('     en regime    (t=32 ms) %8.0f W   alors que Pem = %.0f W\n', ...
    interp1(tv,Pw,32e-3),interp1(t,Pem,32e-3));
