%% RUN_DIAG_WAVE  -  Ou vit l'ecart ? Comparaison des FORMES D'ONDE et du spectre complet
clear; clc;
M=machine_bldc();
D=readmatrix(fullfile(M.FEA.dir,'magnetostique(Magnetic_loading)', ...
    'Calculator Expressions Plot 4.tab'),'FileType','text','Delimiter','\t','NumHeaderLines',1);
fprintf('fichier Br : %d lignes, %d colonnes ; angle %.3f .. %.3f deg\n', ...
    size(D,1),size(D,2),min(D(:,2)),max(D(:,2)));
thF=D(:,2)*pi/180; BrF=D(:,3); [thF,ix]=sort(thF); BrF=BrF(ix);
% retirer un eventuel point duplique a 360 deg
if abs(thF(end)-2*pi)<1e-9 || (thF(end)-thF(1))>2*pi-1e-6, end
fu=@(y,t,k)hypot(2*mean(y(:).*cos(k*t(:))),2*mean(y(:).*sin(k*t(:))));

K=8; nm=ceil(4*K*pi/(M.ws0/M.Rsi));
R=subdomain_mec(M,K,nm,1e7,0);

fprintf('\n=== SPECTRE COMPLET de Br au mi-entrefer (ordres 1..40) ===\n');
fprintf('%5s %10s %10s %9s   %s\n','ordre','FEA','MEC','ecart%','identification');
lbl=containers.Map('KeyType','double','ValueType','char');
lbl(7)='p (fondamental)'; lbl(21)='3p'; lbl(35)='5p';
lbl(8)='|Ns-p|'; lbl(22)='Ns+p'; lbl(23)='|2Ns-p|'; lbl(37)='2Ns+p';
lbl(6)='|Ns-3p|'; lbl(36)='Ns+3p'; lbl(9)='|2Ns-3p|'; lbl(1)='|2Ns-Nm-p|';
for k=1:40
    aF=fu(BrF,thF,k); aM=fu(R.Br,R.thq,k);
    if aF>0.002 || aM>0.002
        s=''; if isKey(lbl,k), s=lbl(k); end
        fprintf('%5d %10.5f %10.5f %+8.1f%%   %s\n',k,aF,aM,100*(aM-aF)/aF,s);
    end
end

% --- forme d'onde sur un pas d'encoche, recalee ---
[~,i0]=max(BrF); th0=thF(i0);
w=abs(angle(exp(1i*(thF-th0))))<=(2*pi/M.Ns);
[~,j0]=max(R.Br); thm=R.thq(j0);
wm=abs(angle(exp(1i*(R.thq-thm))))<=(2*pi/M.Ns);
fprintf('\n=== Creux d''encoche (autour du max de Br) ===\n');
fprintf('FEA : max %.4f  min local %.4f  -> creux %.1f %%\n', ...
    max(BrF(w)),min(BrF(w)),100*(1-min(BrF(w))/max(BrF(w))));
fprintf('MEC : max %.4f  min local %.4f  -> creux %.1f %%\n', ...
    max(R.Br(wm)),min(R.Br(wm)),100*(1-min(R.Br(wm))/max(R.Br(wm))));

% --- energie de l'ecart : localisee aux encoches ? ---
BrM=interp1(R.thq,R.Br,thF,'linear','extrap');
% recalage de phase sur le fondamental
cF=2*mean(BrF.*exp(-1i*M.p*thF)); cM=2*mean(BrM.*exp(-1i*M.p*thF));
sh=angle(cM/cF)/M.p;
BrM=interp1(R.thq,R.Br,mod(thF+sh,2*pi),'linear','extrap');
d=BrF-BrM;
fprintf('\necart RMS global = %.4f T (%.1f %% de Bg1)\n',sqrt(mean(d.^2)),100*sqrt(mean(d.^2))/1.075);

f=figure('Position',[60 60 1150 760],'Color','w','Visible','off');
subplot(2,2,1);
plot(thF*180/pi,BrF,'r','LineWidth',1.3); hold on;
plot(mod(thF+sh,2*pi)*180/pi,BrM,'b','LineWidth',1.1); grid on; xlim([0 360]);
xlabel('\theta [deg]'); ylabel('B_r [T]'); legend('FEA','MEC'); title('(a) B_r mi-entrefer');
subplot(2,2,2);
m=thF*180/pi<=72;
plot(thF(m)*180/pi,BrF(m),'r','LineWidth',1.5); hold on;
plot(thF(m)*180/pi,BrM(m),'b','LineWidth',1.3); grid on;
xlabel('\theta [deg]'); ylabel('B_r [T]'); legend('FEA','MEC'); title('(b) zoom : 3 encoches');
subplot(2,2,3);
plot(thF*180/pi,d,'k','LineWidth',1); grid on; xlim([0 360]);
xlabel('\theta [deg]'); ylabel('\Delta B_r [T]'); title('(c) ecart FEA - MEC');
subplot(2,2,4);
kk=1:40; aF=arrayfun(@(k)fu(BrF,thF,k),kk); aM=arrayfun(@(k)fu(R.Br,R.thq,k),kk);
bar(kk-0.2,aF,0.4,'r'); hold on; bar(kk+0.2,aM,0.4,'b'); grid on; set(gca,'YScale','log');
ylim([1e-4 2]); xlabel('ordre'); ylabel('|B_{r,n}| [T]'); legend('FEA','MEC'); title('(d) spectre');
exportgraphics(f,'FIG_diag_wave.png','Resolution',130);
fprintf('-> FIG_diag_wave.png\n');
