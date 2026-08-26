%% RUN_TORQUE2D  Le sous-modele 2D de bouche ameliore-t-il le COUPLE ?
%  Compare, sur la geometrie ANSYS corrigee et avec le COURANT REEL de la FEA :
%     (a) ouverture = scalaire kfringe identifie
%     (b) ouverture = bouche 2D maillee (sans parametre libre)
%  Reference : T = <somme(e.i)>/omega calcule a 100 % sur donnees FEA
%  (= 4.892 N.m, coherent a 0.5 % avec le P_em/omega d'ANSYS).
%
%  NB SUR L'ALIGNEMENT : en charge FL_a contient L*i, donc caler la phase sur
%  son fondamental introduit un biais. On le CONTROLE par un balayage du
%  decalage ; et la comparaison ENTRE MODELES reste robuste, le biais etant
%  quasi identique pour les deux (leurs fondamentaux sont en phase a <1 deg).
clear; clc;
M=machine_bldc(); Ns=M.Ns; p=M.p; Nm=M.Nm;
if isfile('kfringe_ident.mat'), S=load('kfringe_ident.mat'); kfr=S.kbest; else, kfr=0.5; end
PA=[1 -2 -15 3 14]; PB=[6 -7 -5 8 4]; PC=[11 -12 -10 13 9];
fea=fullfile(M.FEA.dir,'transitoire (en charge)');
om=M.speed*2*pi/60;
fprintf('=== SOUS-MODELE 2D DE BOUCHE : EFFET SUR LE COUPLE ===\n');

%% --- reference FEA -----------------------------------------------------
di=readmatrix(fullfile(fea,'BranchCurrent Plot 1.tab'),'FileType','text','NumHeaderLines',1);
dl=readmatrix(fullfile(fea,'Output Variables Plot 1.tab'),'FileType','text','NumHeaderLines',1);
ds=readmatrix(fullfile(fea,'Speed Plot 1.tab'),'FileType','text','NumHeaderLines',1);
n=size(di,1); w=round(n*0.6):n; t=di(w,1)*1e-3; ia=di(w,2); ib=di(w,3); ic=di(w,4);
nl=size(dl,1); wl=round(nl*0.6):nl; tl=dl(wl,1)*1e-3;
nsp=size(ds,1); nF=mean(ds(round(nsp*0.75):nsp,2)); fe=nF*Nm/120; omF=nF*2*pi/60;
tq=linspace(max(t(1),tl(1)),min(t(end),tl(end)),4001);
FLa=interp1(tl,dl(wl,2),tq); FLc=interp1(tl,dl(wl,3),tq); FLb=interp1(tl,dl(wl,4),tq);
iaq=interp1(t,ia,tq); ibq=interp1(t,ib,tq); icq=interp1(t,ic,tq);
T_ref=abs(mean(-gradient(FLa,tq).*iaq-gradient(FLb,tq).*ibq-gradient(FLc,tq).*icq)/omF);
% ALIGNEMENT : en charge FL_a = lambda_aimant + L_d*i_a (systeme equilibre :
% L_aa*i_a + M*(i_b+i_c) = (L_aa-M)*i_a = L_d*i_a). Avec L_d = 52.34 mH (FEA)
% et i_a crete 2.53 A, le terme L_d*i vaut 0.13 Wb pour un lambda_aimant de
% 0.26 Wb : caler la phase sur FL_a BRUT biaise l'alignement de ~27 deg.
% On retranche donc L_d*i_a (grandeur FEA) avant de prendre le fondamental.
FLa_m = FLa - M.FEA.Ld*iaq;
psiF   = angle(mean(FLa_m.*exp(-1i*2*pi*fe*tq)));
psiF_b = angle(mean(FLa  .*exp(-1i*2*pi*fe*tq)));
fprintf('  Alignement : phase corrigee %.1f deg vs brute %.1f deg -> biais %.1f deg elec\n',...
        psiF*180/pi,psiF_b*180/pi,(psiF_b-psiF)*180/pi);
% back-EMF FEA a vide (essai Back_emf) pour crete/eff
dn=readmatrix(fullfile(M.FEA.dir,'transitoire (Back_emf)','Output Variables Plot 4.tab'),'FileType','text','NumHeaderLines',1);
en=-gradient(dn(:,3),dn(:,1)*1e-9);
fprintf('  Reference FEA : T = %.3f N.m | back-EMF crete %.1f V, eff %.1f V\n',...
        T_ref,max(abs(en)),sqrt(mean(en.^2)));
fprintf('  Courant FEA : %.1f tr/min, %.3f A eff\n\n',nF,sqrt(mean(ia.^2)));

%% --- boucle sur les deux modeles ---------------------------------------
amp=@(B,th,k) 2*abs(mean(B.*exp(-1i*k*th)));
Np=721; res=struct(); off=(-30:2:30)*pi/180;      % balayage de decalage elec
for c=1:2
    if c==1, nom='(a) scalaire kfringe'; R=cogging_mec (M,840,0,Np,M.muI,kfr,2*pi/p);
    else,    nom='(b) bouche 2D (nm=8)'; R=cogging_mec2(M,840,8,Np,M.muI,2*pi/p); end
    lam=@(P) M.Ntc*sum(sign(P(:)).*R.PhiT(abs(P(:)),:),1);
    lamA=lam(PA); ang=p*R.phis;
    eA=-gradient(lamA,R.phis)*om; eB=-gradient(lam(PB),R.phis)*om; eC=-gradient(lam(PC),R.phis)*om;
    a0=ang(1:end-1); pk=@(y)[y(1:end-1) y(1)];
    ip=@(y,x) interp1([a0 2*pi],pk(y),mod(x,2*pi),'linear');
    psiM=angle(mean(lamA.*exp(-1i*ang)));
    the0=2*pi*fe*tq - psiF + psiM;
    Tv=zeros(size(off));
    for k=1:numel(off)
        th=the0+off(k);
        Tv(k)=mean(ip(eA,th).*iaq+ip(eB,th).*ibq+ip(eC,th).*icq)/om;
    end
    T0=Tv(off==0);
    thm=R.thq(1:end-1); Brm=R.Br(1:end-1);
    res(c).nom=nom; res(c).T0=T0; res(c).Tv=Tv; res(c).Tmax=max(Tv);
    res(c).epk=max(abs(eA)); res(c).erms=sqrt(mean(eA.^2));
    res(c).a8=amp(Brm,thm,8); res(c).a22=amp(Brm,thm,22); res(c).a23=amp(Brm,thm,23);
    res(c).psiM=psiM;
end

%% --- resultats ---------------------------------------------------------
aF=[0.0197 0.0328 0.0202];
fprintf('  %-22s %9s %9s %9s %9s %9s %9s\n','modele','a8','a22','a23','e_pk(V)','T(N.m)','ecart');
fprintf('  %-22s %9.4f %9.4f %9.4f %9.1f %9.3f %8s\n','FEA',aF(1),aF(2),aF(3),max(abs(en)),T_ref,'--');
for c=1:2
    fprintf('  %-22s %9.4f %9.4f %9.4f %9.1f %9.3f %+8.1f%%\n',res(c).nom,...
        res(c).a8,res(c).a22,res(c).a23,res(c).epk,res(c).T0,(res(c).T0-T_ref)/T_ref*100);
end
fprintf('\n  Harmoniques (ecart vs FEA) :\n');
for c=1:2
    fprintf('    %-22s a8 %+6.1f %% | a22 %+6.1f %% | a23 %+6.1f %%\n',res(c).nom,...
        (res(c).a8-aF(1))/aF(1)*100,(res(c).a22-aF(2))/aF(2)*100,(res(c).a23-aF(3))/aF(3)*100);
end
fprintf('\n  Sensibilite a l''alignement (balayage +/-30 deg elec) :\n');
for c=1:2
    fprintf('    %-22s T = %.3f .. %.3f N.m (max %.3f)\n',res(c).nom,min(res(c).Tv),max(res(c).Tv),res(c).Tmax);
end
fprintf('    ecart de phase entre les 2 modeles : %.2f deg elec\n',(res(2).psiM-res(1).psiM)*180/pi);
fprintf('\n  => le sous-modele 2D change le couple de %+.1f point(s) (%.3f -> %.3f N.m)\n',...
        (res(2).T0-res(1).T0)/T_ref*100,res(1).T0,res(2).T0);

figure('Name','Bouche 2D : effet sur le couple','Color','w','Position',[60 60 1000 400]);
subplot(1,2,1);
plot(off*180/pi,res(1).Tv,'b-o','LineWidth',1.3,'MarkerSize',3); hold on;
plot(off*180/pi,res(2).Tv,'g-s','LineWidth',1.3,'MarkerSize',3);
yline(T_ref,'r--','LineWidth',1.5); grid on;
xlabel('decalage d''alignement (deg elec)'); ylabel('T (N.m)');
legend('scalaire','bouche 2D','FEA','Location','best'); title('Couple vs alignement');
subplot(1,2,2);
bar([res(1).a8 res(2).a8 aF(1); res(1).a22 res(2).a22 aF(2); res(1).a23 res(2).a23 aF(3)]);
grid on; set(gca,'XTickLabel',{'a_8','a_{22}','a_{23}'}); ylabel('|B_r| (T)');
legend('scalaire','bouche 2D','FEA','Location','best'); title('Harmoniques de denture');
saveas(gcf,fullfile(fileparts(mfilename('fullpath')),'FIG_torque2d.png'));
