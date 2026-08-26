%% RUN_DIAG_ANGLE - ou tombe la fenetre de conduction FEA par rapport a la FEM ?
%  Le MEC donne 30 % de courant en moins que la FEA au meme point, alors que
%  R, Vdc, L et lambda_pm concordent a 2 % pres. Il ne reste qu'une variable :
%  la position de la fenetre de conduction par rapport a la FEM. On la mesure
%  directement sur l'essai en charge, sans passer par le fondamental.
clear; clc;
M=machine_bldc(); p=M.p; fea=M.FEA.dir; ol=fullfile(fea,'transitoire (en charge)');
rd=@(f)readmatrix(f,'FileType','text','Delimiter','\t','NumHeaderLines',1);
dI=rd(fullfile(ol,'BranchCurrent Plot 1.tab')); dS=rd(fullfile(ol,'Speed Plot 1.tab'));
dL=rd(fullfile(ol,'Plot 1_loss.tab'));
t=dI(:,1)*1e-3; i=dI(:,2:4);
om=interp1(dS(:,1)*1e-3,dS(:,2),t,'linear','extrap')*2*pi/60;
th=cumtrapz(t,om)*p;                                   % angle ELECTRIQUE ANSYS

%% ---- FEM a vide de reference (essai Back_emf, 1500 tr/min) -------------
d4=rd(fullfile(fea,'transitoire (Back_emf)','Output Variables Plot 4.tab'));
pe=mod(d4(:,2)*p,360); FLa=d4(:,3); FLb=d4(:,7); FLc=d4(:,5);
[pe,o]=sort(pe); FLa=FLa(o); FLb=FLb(o); FLc=FLc(o);
[pe,u]=unique(pe); FLa=FLa(u); FLb=FLb(u); FLc=FLc(u);
g=linspace(0,360,721); g(end)=[];
La=interp1([pe;360],[FLa;FLa(1)],g,'pchip');
Lb=interp1([pe;360],[FLb;FLb(1)],g,'pchip');
Lc=interp1([pe;360],[FLc;FLc(1)],g,'pchip');
dd=@(x) p*real(ifft(1i*[0:359 -360:-1].*fft(x)))/(pi/180)/  (180/pi);  % d/dth_m
sp=@(x) p*real(ifft(1i*[0:359,-360:-1].*fft(x)));       % d/dth_e * p
psa=sp(La); psb=sp(Lb); psc=sp(Lc);                     % Wb/rad mec

fprintf('=== Fenetre de conduction FEA vs FEM ===\n');
fprintf('  psi_ab = psi_a - psi_b : maximum %.3f Wb/rad a theta_e = %.1f deg\n', ...
    max(psa-psb),g(find(psa-psb==max(psa-psb),1)));

%% ---- ou conduit reellement la FEA ? -----------------------------------
m=t>0.7*t(end);
lab={'a','b','c'}; fprintf('\n  %8s %8s %7s %7s %7s   paire   psi_paire\n', ...
    't(ms)','th_e','ia','ib','ic');
sel=find(m); sel=sel(1:2:min(30,numel(sel)));
psel=zeros(numel(sel),1); k2=0;
for k=sel.'
    [~,ip]=max(i(k,:)); [~,in]=min(i(k,:));
    tk=mod(th(k)*180/pi,360); j=round(tk/0.5)+1; j=min(max(j,1),720);
    ps=[psa(j) psb(j) psc(j)]; pv=ps(ip)-ps(in);
    k2=k2+1; psel(k2)=pv;
    fprintf('  %8.2f %8.1f %7.2f %7.2f %7.2f   %s+/%s-  %8.3f\n', ...
        t(k),tk,i(k,1),i(k,2),i(k,3),lab{ip},lab{in},pv);
end
fprintf('\n  psi de la paire conductrice : moyenne %.3f Wb/rad (max possible %.3f)\n', ...
    mean(psel),max(psa-psb));
fprintf('  taux d''utilisation de la FEM : %.1f %%\n',100*mean(psel)/max(psa-psb));

%% ---- couple et courant qui en decoulent -------------------------------
nF=mean(dS(round(0.75*end):end,2)); omF=nF*2*pi/60;
TF=mean(dL(round(0.75*end):end,7))/omF;
ic=mean(abs(i(m,1))+abs(i(m,2))+abs(i(m,3)))/2;
fprintf('\n  courant conducteur moyen FEA  : %.3f A\n',ic);
fprintf('  couple FEA                    : %.3f N.m -> psi_eff = %.3f Wb/rad\n',TF,TF/ic);
fprintf('  FEM correspondante            : %.0f V a %.0f tr/min\n',TF/ic*omF,nF);
fprintf('  tension motrice Vdc - e       : %.0f V -> i = %.2f A (quasi statique)\n', ...
    500-TF/ic*omF,(500-TF/ic*omF)/(2*M.Rph));
