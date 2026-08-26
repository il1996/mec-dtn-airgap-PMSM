%% RUN_FIELD  Validation du champ a vide analytique (DtN etendu) vs FEA.
clear; clc;
M = machine_bldc();

% --- champ analytique sur un tour mecanique ---
th = linspace(0, 2*pi, 3601);
F  = magnet_subdomain(M, th, 40);          % 40 harmoniques impairs
Br = F.Br;  Bt = F.Bt;

Bmean_a = mean(abs(Br));
Bpeak_a = max(Br);
% fondamental (ordre spatial p=7) par projection
c7 = 2*mean(Br.*cos(M.p*th));  s7 = 2*mean(Br.*sin(M.p*th));
Bg1_a = hypot(c7,s7);
Bt_a  = max(abs(Bt));

% --- FEA : Calculator Expressions Plot 4 (Br) et Plot 2 (Bt) ---
fea = M.FEA.dir;
d4 = readmatrix(fullfile(fea,'magnetostique(Magnetic_loading)','Calculator Expressions Plot 4.tab'), ...
                'FileType','text','NumHeaderLines',1);
angF = d4(:,2)*pi/180;  BrF = d4(:,3);
d2 = readmatrix(fullfile(fea,'magnetostique(Magnetic_loading)','Calculator Expressions Plot 2.tab'), ...
                'FileType','text','NumHeaderLines',1);
BtF = d2(:,2);

Bmean_f = mean(abs(BrF));
Bpeak_f = max(BrF);
% fondamental FEA (l'angle FEA couvre 'span' deg mec)
spanF = angF(end)-angF(1);
c7f = 2*mean(BrF.*cos(M.p*angF));  s7f = 2*mean(BrF.*sin(M.p*angF));
Bg1_f = hypot(c7f,s7f);
Bt_f  = max(abs(BtF));

fprintf('\n=== CHAMP A VIDE : sous-domaine analytique (DtN etendu) vs FEA ===\n');
fprintf('  FEA : angle couvert %.1f deg mec, %d points\n', spanF*180/pi, numel(angF));
fprintf('  %-26s %10s %10s %9s\n','Grandeur','MEC(ana)','FEA','Ecart');
row=@(n,a,f) fprintf('  %-26s %10.4f %10.4f %+8.1f%%\n',n,a,f,(a-f)/f*100);
row('B entrefer moyen |Br| (T)', Bmean_a, Bmean_f);
row('B entrefer crete (T)',      Bpeak_a, Bpeak_f);
row('Fondamental Bg1 (T)',       Bg1_a,  Bg1_f);
row('B tangentiel crete (T)',    Bt_a,   Bt_f);

% --- figure de superposition ---
figure('Name','Champ a vide MEC vs FEA','Color','w');
plot(th*180/pi, Br,'b','LineWidth',1.4); hold on;
plot(angF*180/pi, BrF,'r--','LineWidth',1.1); grid on;
xlabel('\theta (deg mec)'); ylabel('B_r entrefer (T)'); xlim([0 360]);
legend('MEC sous-domaine (lisse)','FEA (avec encoches)','Location','best');
title(sprintf('B_r mi-entrefer : moy %.3f/%.3f T, Bg1 %.3f/%.3f T (MEC/FEA)',...
      Bmean_a,Bmean_f,Bg1_a,Bg1_f));
saveas(gcf, fullfile(fileparts(mfilename('fullpath')),'FIG_field.png'));
fprintf('  (figure FIG_field.png sauvee)\n');
