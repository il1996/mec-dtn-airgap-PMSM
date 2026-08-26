%% RUN_DIAG_ETA - definition EXACTE du rendement et du bilan de puissance FEA
clear; clc;
M=machine_bldc(); ol=fullfile(M.FEA.dir,'transitoire (en charge)');
rd=@(f)readmatrix(f,'FileType','text','Delimiter','\t','NumHeaderLines',1);
w=@(x)x(max(1,round(0.75*numel(x))):end);          % dernier quart

dL =rd(fullfile(ol,'Plot 1_loss.tab'));            % t[ns] Solid Core p_cu P_out Pro Pem
dE =rd(fullfile(ol,'Output Variables Plot 3.tab'));% t[ms] efficiency
dS =rd(fullfile(ol,'Speed Plot 1.tab'));
dIt=rd(fullfile(ol,'BranchCurrent Plot 2.tab'));   % t[ms] I(Vi_t)
dI =rd(fullfile(ol,'BranchCurrent Plot 1.tab'));
dV =rd(fullfile(ol,'NodeVoltage Plot 2.tab'));     % Ie_a..c [mV], Iv_a..c [mV]
dFL=rd(fullfile(ol,'Output Variables Plot 1.tab'));
dB =rd(fullfile(ol,'Output Variables Plot 2.tab'));% BEMF

Ps=mean(w(dL(:,2)))/1e3; Pc=mean(w(dL(:,3)))/1e3; Pj=mean(w(dL(:,4)));
Po=mean(w(dL(:,5)));     Pr=mean(w(dL(:,6)));     Pe=mean(w(dL(:,7)));
n =mean(w(dS(:,2))); om=n*2*pi/60;
It=mean(w(dIt(:,2))); eta=mean(w(dE(:,2)));
Irms=sqrt(mean(w(dI(:,2)).^2));

fprintf('=== Bilan de puissance FEA en regime (dernier quart) ===\n');
fprintf('  vitesse                  %10.1f tr/min  (%.2f rad/s)\n',n,om);
fprintf('  Pem  (couple elm x om)   %10.2f W\n',Pe);
fprintf('  P_out                    %10.2f W\n',Po);
fprintf('  Pro  (frottement)        %10.2f W\n',Pr);
fprintf('  p_cu (Joule)             %10.2f W\n',Pj);
fprintf('  CoreLoss (fer)           %10.3f W\n',Pc);
fprintf('  SolidLoss (aimant)       %10.3f W\n',Ps);
fprintf('  I(Vi_t) moyen (bus)      %10.4f A   -> Vdc*I = %.2f W\n',It,500*It);
fprintf('  I phase rms              %10.4f A   -> 3*R*I^2 = %.2f W\n',Irms,3*M.Rph*Irms^2);
fprintf('  rendement FEA            %10.2f %%\n',100*eta);

fprintf('\n=== Quelle formule reproduit ce rendement ? ===\n');
cand={ ...
 'P_out/(P_out+Pj+Pc+Ps)',            Po/(Po+Pj+Pc+Ps); ...
 'P_out/(P_out+Pj+Pc+Ps+Pro)',        Po/(Po+Pj+Pc+Ps+Pr); ...
 'Pem/(Pem+Pj+Pc+Ps)',                Pe/(Pe+Pj+Pc+Ps); ...
 '(Pem-Pro)/(Pem+Pj+Pc+Ps)',          (Pe-Pr)/(Pe+Pj+Pc+Ps); ...
 'P_out/(Vdc*I_bus)',                 Po/(500*It); ...
 'Pem/(Vdc*I_bus)',                   Pe/(500*It); ...
 '(Pem-Pro)/(Vdc*I_bus)',             (Pe-Pr)/(500*It)};
for k=1:size(cand,1)
    fprintf('  %-32s = %7.3f %%   (ecart %+6.3f pt)\n', ...
        cand{k,1},100*cand{k,2},100*(cand{k,2}-eta));
end

fprintf('\n=== Coherence interne ===\n');
fprintf('  Pem - Pro - P_out = %+.3f W  (doit etre ~0 en regime)\n',Pe-Pr-Po);
fprintf('  Vdc*Ibus - (Pem+Pj+Pc+Ps) = %+.3f W\n',500*It-(Pe+Pj+Pc+Ps));
fprintf('  couple  Pem/om = %.4f N.m ; 4.8+6.07927e-4*om = %.4f N.m\n', ...
    Pe/om,4.8+6.07927101854027e-4*om);

fprintf('\n=== Formule EXACTE lue dans BLDC.aedt ===\n');
fprintf('  p_cu       = $Rph*(i_a^2+i_b^2+i_c^2)\n');
fprintf('  P_out      = Moving1.Torque * Moving1.Speed\n');
fprintf('  efficiency = P_out/(P_out*1.02 + p_cu + CoreLoss + SolidLoss)\n');
tl=dL(:,1)*1e-9*1e3; te=dE(:,1);
for k=[5 7]
    ek=dL(:,k)./(dL(:,k)*1.02+dL(:,4)+dL(:,3)/1e3+dL(:,2)/1e3);
    ei=interp1(tl,ek,te,'linear','extrap');
    m=te>te(end)*0.5;
    fprintf('  colonne %d -> ecart RMS avec la trace efficiency : %.4f pt\n', ...
        k,100*sqrt(mean((ei(m)-dE(m,2)).^2)));
end
% puissance convertie Somme(e_k*i_k) a partir des sondes FEA
ti=dI(:,1); ev=interp1(dV(:,1)*1e-9*1e3,dV(:,2:4)/1e3,ti,'linear','extrap');
Pconv=sum(ev.*dI(:,2:4),2);
fprintf('  somme(e_k*i_k) mesuree = %.2f W  -> couple %.4f N.m\n', ...
    mean(w(Pconv)),mean(w(Pconv))/om);

fprintf('\n=== Grandeurs disponibles pour la comparaison temporelle ===\n');
fprintf('  Speed Plot 1          : vitesse(t)          %d pts\n',size(dS,1));
fprintf('  BranchCurrent Plot 1  : i_a,i_b,i_c(t)      %d pts\n',size(dI,1));
fprintf('  BranchCurrent Plot 2  : i_bus(t)            %d pts\n',size(dIt,1));
fprintf('  NodeVoltage Plot 2    : e_abc, v_abc (mV)   %d pts\n',size(dV,1));
fprintf('  Output Variables 1    : FL_a,FL_b,FL_c(t)   %d pts\n',size(dFL,1));
fprintf('  Output Variables 2    : BEMF(t)             %d pts\n',size(dB,1));
fprintf('  Output Variables 3    : rendement(t)        %d pts\n',size(dE,1));
fprintf('  Plot 1_loss           : Solid/Core/pcu/Pout/Pro/Pem  %d pts\n',size(dL,1));
fprintf('\n  crete de courant au demarrage : %.2f A a t = %.2f ms\n', ...
    max(abs(dI(:,2))),dI(find(abs(dI(:,2))==max(abs(dI(:,2))),1),1));
fprintf('  temps de montee 0->95%% de %.0f tr/min : %.2f ms\n', ...
    n,dS(find(dS(:,2)>=0.95*n,1),1));
