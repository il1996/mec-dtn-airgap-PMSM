%% RUN_COG_FEA  Reference FEA du couple de detente : ANALYSE CRITIQUE
%  Deux essais ANSYS donnent le couple de detente :
%    (a) magnetostatique parametrique : Torque Plot 1.tab, 361 pts, pas 1 deg
%    (b) transitoire back-EMF         : Torque Plot.tab,   pas 0.0571 ms
%
%  L'ordre de detente vaut LCM(15,14) = 210 -> periode 1.714 deg MECANIQUE.
%  Un pas de 1 deg mecanique serait SOUS-NYQUIST (il faudrait < 0.857 deg).
%  Deux hypotheses sur la variable balayee $alfa_R :
%    H1 : alfa_R MECANIQUE -> 210 alias sur l'ordre |210-360| = 150
%         (tous les alias de la famille 210 tombent sur des multiples de
%          gcd(210,360) = 30)
%    H2 : alfa_R ELECTRIQUE -> 360 elec = 51.43 mec = 2 pas polaires
%         -> 30 periodes de detente sur le record -> raie NETTE a l'ordre 30
%          (12 points/periode : PARFAITEMENT resolu)
%  Le spectre tranche. Ce script produit le verdict et la reference retenue.
clear; clc;
M=machine_bldc(); fea=M.FEA.dir;
Ns=M.Ns; Nm=M.Nm; p=M.p; NL=lcm(Ns,Nm);

fprintf('=== REFERENCE FEA DU COUPLE DE DETENTE : ANALYSE CRITIQUE ===\n');
fprintf('  Ordre de detente attendu : LCM(%d,%d) = %d -> periode %.4f deg mec\n',Ns,Nm,NL,360/NL);

%% ---------- (a) magnetostatique parametrique ----------
d = readmatrix(fullfile(fea,'magnetostique(Magnetic_loading)','Torque Plot 1.tab'), ...
               'FileType','text','NumHeaderLines',1);
al=d(:,1); Tm=d(:,2);
al=al(1:360); Tm=Tm(1:360);                 % 0..359 (360 est le repli de 0)
N=360;
fprintf('\n[a] MAGNETOSTATIQUE : %d points, pas 1 deg, p-p = %.2f mN.m, RMS = %.2f\n',...
        N, max(Tm)-min(Tm), std(Tm));
fprintf('    moyenne = %.3f mN.m (doit etre ~0 : pas de courant)\n', mean(Tm));

X=fft(Tm)/N; amp=2*abs(X(2:181));           % ordres 1..180
[~,ki]=sort(amp,'descend');
fprintf('    8 raies dominantes :\n');
for t=1:8, fprintf('      ordre %3d : %7.3f mN.m\n',ki(t),amp(ki(t))); end

% --- test de coherence : energie sur les multiples de 30 ---
bins=(1:180).'; is30=mod(bins,30)==0; isodd=mod(bins,2)==1;
E30=sum(amp(is30).^2); Etot=sum(amp.^2);
fprintf('    Energie sur les multiples de 30 : %.1f %% (bruit blanc -> %.1f %%)\n',...
        E30/Etot*100, sum(is30)/180*100);
fprintf('    Energie sur les ordres IMPAIRS  : %.2e %% (symetrie 180 deg)\n',...
        sum(amp(isodd).^2)/Etot*100);

% --- verdict H1 / H2 ---
a150=amp(150); a30=amp(30); a60=amp(60); a90=amp(90); a120=amp(120);
fprintf('    Raies cles : ordre 30 = %.3f | 60 = %.3f | 90 = %.3f | 120 = %.3f | 150 = %.3f mN.m\n',...
        a30,a60,a90,a120,a150);
if a30 > 3*max([a60 a90 a120 a150])
    hyp='H2 (alfa_R ELECTRIQUE : detente RESOLUE a l''ordre 30 du record)';
    cog_pp_mag = (max(Tm)-min(Tm));
    cog_amp_mag = a30;
elseif a150 > 3*max([a30 a60 a90 a120])
    hyp='H1 (alfa_R MECANIQUE : detente ALIASEE 210->150, amplitude conservee)';
    cog_amp_mag = a150; cog_pp_mag = 2*a150;
else
    hyp='AUCUNE raie dominante -> spectre LARGE BANDE = BRUIT DE MAILLAGE';
    cog_amp_mag = NaN; cog_pp_mag = NaN;
end
fprintf('    >> VERDICT : %s\n',hyp);

%% ---------- (b) transitoire back-EMF ----------
dt = readmatrix(fullfile(fea,'transitoire (Back_emf)','Torque Plot.tab'), ...
                'FileType','text','NumHeaderLines',1);
tms=dt(:,1); Tt=dt(:,2);
angm = tms*(M.speed*360/60000);             % deg mecaniques (1500 tr/min -> 9 deg/ms)
dang = mean(diff(angm));
fprintf('\n[b] TRANSITOIRE : %d points, pas %.4f deg mec (%.2f pts/periode de detente)\n',...
        numel(Tt), dang, (360/NL)/dang);
fprintf('    span = %.2f deg mec = %.1f periodes | p-p brut = %.2f mN.m, moyenne = %.2f\n',...
        angm(end)-angm(1), (angm(end)-angm(1))/(360/NL), max(Tt)-min(Tt), mean(Tt));

% ajustement moindres carres aux ordres 210 et 420 (record non entier -> pas de FFT)
th=angm*pi/180;
A=[ones(size(th)), th, cos(NL*th), sin(NL*th), cos(2*NL*th), sin(2*NL*th)];
c=A\Tt;
amp210=hypot(c(3),c(4)); amp420=hypot(c(5),c(6));
resid=Tt-A*c;
fprintf('    Ajustement LSQ : |ordre %d| = %.2f mN.m, |ordre %d| = %.2f mN.m\n',NL,amp210,2*NL,amp420);
fprintf('    Residu RMS = %.2f mN.m (bruit) | signal RMS ajuste = %.2f mN.m\n',...
        std(resid), std(A(:,3:6)*c(3:6)));
cog_pp_tr = 2*hypot(amp210,amp420);

%% ---------- synthese ----------
fprintf('\n--- SYNTHESE DE LA REFERENCE FEA ---\n');
fprintf('  magnetostatique : amplitude ordre %d = %.2f mN.m -> c-c ~ %.2f mN.m\n',...
        NL, cog_amp_mag, cog_pp_mag);
fprintf('  transitoire     : amplitude ordre %d = %.2f mN.m -> c-c ~ %.2f mN.m\n',...
        NL, amp210, cog_pp_tr);
fprintf('  couple nominal  : %.3f N.m -> detente = %.3f %% du nominal\n',...
        M.FEA.T_load, cog_pp_mag/1000/M.FEA.T_load*100);

figure('Name','Reference FEA du couple de detente','Color','w','Position',[60 60 1150 620]);
subplot(2,2,1);
plot(al,Tm,'r','LineWidth',1.0); grid on; xlabel('\alpha_R (deg)'); ylabel('T (mN.m)');
title(sprintf('(a) Magnetostatique : p-p %.1f mN.m',max(Tm)-min(Tm)));
subplot(2,2,2);
stem(1:180,amp,'filled','MarkerSize',2.5); hold on;
stem(bins(is30),amp(is30),'r','filled','MarkerSize',4); grid on;
xlabel('ordre spatial du record'); ylabel('|T_k| (mN.m)');
title('Spectre (rouge = multiples de 30)'); xlim([0 180]);
subplot(2,2,3);
plot(al(1:60),Tm(1:60),'r-o','MarkerSize',3,'LineWidth',1.1); grid on;
xlabel('\alpha_R (deg)'); ylabel('T (mN.m)'); title('(a) zoom 0-60 deg');
subplot(2,2,4);
plot(angm,Tt,'b-o','MarkerSize',3,'LineWidth',1.0); hold on;
plot(angm,A*c,'k--','LineWidth',1.4); grid on;
xlabel('angle mec (deg)'); ylabel('T (mN.m)');
legend('FEA transitoire','ajustement ordre 210+420','Location','best');
title(sprintf('(b) Transitoire : |210| = %.1f mN.m',amp210));
sgtitle('Couple de detente : critique des deux references FEA');
saveas(gcf,fullfile(fileparts(mfilename('fullpath')),'FIG_cog_fea.png'));
fprintf('  (figure FIG_cog_fea.png sauvee)\n');
