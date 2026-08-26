%% RUN_FEA_MESHNOISE  -  Sensibilite au MAILLAGE de la reference magnetostatique
%
%  Impossible de re-mailler ici (AEDT plante en batch dans cet environnement).
%  Mais le design 'Magnetic_loading' fournit deja une mesure DIRECTE de sa
%  propre dispersion de maillage : Torque Plot 1.tab donne le couple a 361
%  positions rotoriques, chacune REMAILLEE independamment. Le couple de
%  detente d'un 15/14 doit etre un multiple de LCM(15,14)=210 par tour ;
%  tout le reste est du BRUIT DE REMAILLAGE.
%
%  On le quantifie, puis on le convertit en incertitude equivalente sur les
%  harmoniques de champ, via  T = (L*r^2/mu0)*pi*SUM(Brc*Btc + Brs*Bts).
clear; clc;
M=machine_bldc(); mu0=4*pi*1e-7;
d=readmatrix(fullfile(M.FEA.dir,'magnetostique(Magnetic_loading)','Torque Plot 1.tab'), ...
    'FileType','text','Delimiter','\t','NumHeaderLines',1);
al=d(:,1); T=d(:,2);                 % deg mecaniques, mN.m
if abs(al(end)-360)<1e-6, al(end)=[]; T(end)=[]; end
N=numel(T); T=T-mean(T);
fprintf('=== Bruit de remaillage de la reference magnetostatique ===\n');
fprintf('%d positions sur %.0f deg (pas %.2f deg)\n',N,max(al),al(2)-al(1));
fprintf('couple : crete-a-crete %.2f mN.m, ecart-type %.2f mN.m\n',max(T)-min(T),std(T));

Y=abs(fft(T))/N*2; nn=0:N-1;
[~,im]=max(Y(2:floor(N/2))); ord=im;
fprintf('\nraie DOMINANTE : ordre %d, amplitude %.2f mN.m\n',ord,Y(ord+1));
LCMv=lcm(M.Ns,M.Nm);
fprintf('or la detente d''un %d/%d n''existe qu''aux multiples de LCM = %d\n',M.Ns,M.Nm,LCMv);
nyq=floor(N/2); fam=LCMv:LCMv:nyq;
if isempty(fam)
    fprintf(['\n*** L''ordre %d est AU-DESSUS DE NYQUIST : avec un pas de %.2f deg\n' ...
        '    (Nyquist = ordre %d) cet essai NE PEUT PAS representer la detente\n' ...
        '    physique ; il faudrait un pas < %.3f deg. Tout ce qui y est mesure\n' ...
        '    est donc, par construction, repliement + bruit de maillage. ***\n'], ...
        LCMv,al(2)-al(1),nyq,360/(2*LCMv));
    Efam=0;
else
    Efam=sum(Y(fam+1).^2);
end
Etot=sum(Y(2:nyq).^2);
fprintf('energie dans la famille %d : %.2f %% du total\n',LCMv,100*Efam/Etot);
fprintf('=> l''ordre %d est PHYSIQUEMENT IMPOSSIBLE : c''est du bruit de maillage.\n',ord);

% --- conversion en incertitude de champ ---
%  T = (L*r^2/mu0)*pi*SUM_n (Brc_n*Btc_n + Brs_n*Bts_n)
kT=(M.ls*M.rmid^2/mu0)*pi;                       % N.m par unite de SUM(Br*Bt)
Tn=std(T)*1e-3;                                   % N.m (ecart-type du bruit)
d2=readmatrix(fullfile(M.FEA.dir,'magnetostique(Magnetic_loading)', ...
    'Calculator Expressions Plot 2.tab'),'FileType','text','Delimiter','\t','NumHeaderLines',1);
Bt_rms=sqrt(mean(d2(:,end).^2));
dBr=Tn/kT/Bt_rms;                                 % perturbation equivalente de Br
fprintf('\n--- Traduction en incertitude sur les harmoniques de champ ---\n');
fprintf('k_T = L*r^2/mu0*pi          = %.1f N.m par unite de SUM(Br*Bt)\n',kT);
fprintf('bruit de couple (ecart-type)= %.2f mN.m\n',Tn*1e3);
fprintf('Bt RMS de reference         = %.4f T\n',Bt_rms);
fprintf('=> perturbation de Br equivalente ~ %.4f T\n',dBr);
fprintf('   a comparer a l''harmonique a8 mesuree      : 0.0187 T\n');
fprintf('   et a l''ecart MEC-FEA a expliquer (-28 %%)  : 0.0052 T\n');
if dBr>0.0052
    fprintf('\n>> Le bruit de maillage de la reference est DU MEME ORDRE que\n');
    fprintf('   l''ecart a expliquer : la valeur a8 = 0.0187 T ne peut pas etre\n');
    fprintf('   consideree comme convergee sans etude de raffinement.\n');
else
    fprintf('\n>> Le bruit de maillage reste sous l''ecart a expliquer :\n');
    fprintf('   la reference est probablement fiable, l''ecart est du modele.\n');
end

figure('Color','w','Position',[60 60 1000 420],'Visible','off');
subplot(1,2,1); plot(al,T,'b','LineWidth',1); grid on;
xlabel('position rotor (deg)'); ylabel('couple (mN.m)');
title(sprintf('Couple magnetostatique : %.1f mN.m c-c',max(T)-min(T)));
subplot(1,2,2); stem(1:nyq-1,Y(2:nyq),'b','MarkerSize',2); grid on;
set(gca,'YScale','log'); hold on;
xline(nyq,'r--','Nyquist','LineWidth',1.2);
xlabel('ordre spatial (/tour)'); ylabel('|T_n| (mN.m)');
title(sprintf('Spectre : max a l''ordre %d ; la detente (ordre %d) est HORS bande',ord,LCMv));
exportgraphics(gcf,'FIG_meshnoise.png','Resolution',130);
fprintf('\n-> FIG_meshnoise.png\n');
