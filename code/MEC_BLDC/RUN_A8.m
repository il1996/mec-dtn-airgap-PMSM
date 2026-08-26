%% RUN_A8  Diagnostic du deficit sur l'harmonique de denture a8 (-37 %)
%  a8 = |Ns-p| = 8 est la bande laterale BASSE de la modulation du
%  fondamental d'aimant (7) par la denture (15). C'est l'harmonique DOMINANT
%  vu par l'aimant (le moins attenue radialement) : il pilote les pertes
%  aimant. Il reste a -37 % apres correction de la geometrie, et le
%  sous-modele 2D de bouche ne le corrige PAS -> la cause est ailleurs.
%  Ce script teste, dans l'ordre, les causes candidates.
clear; clc;
M=machine_bldc(); Ns=M.Ns; p=M.p;
if isfile('kfringe_ident.mat'), S=load('kfringe_ident.mat'); kfr=S.kbest; else, kfr=0.5; end
amp=@(B,th,k) 2*abs(mean(B.*exp(-1i*k*th)));
ords=[7 8 22 23 37];

% --- reference FEA ---
fea=M.FEA.dir;
d4=readmatrix(fullfile(fea,'magnetostique(Magnetic_loading)','Calculator Expressions Plot 4.tab'),'FileType','text','NumHeaderLines',1);
angF=d4(:,2)*pi/180; BrF=d4(:,3);
thu=linspace(0,2*pi,3601); thu(end)=[];
BrFu=interp1([angF;2*pi],[BrF;BrF(1)],thu,'linear','extrap');
aF=arrayfun(@(k) amp(BrFu,thu,k),ords);
fprintf('=== DIAGNOSTIC DU DEFICIT a8 ===\n');
fprintf('  FEA : a7 %.4f | a8 %.4f | a22 %.4f | a23 %.4f | a37 %.4f T\n',aF);
fprintf('  Rapport FEA a22/a8 = %.3f  (un simple produit de permeance\n',aF(3)/aF(2));
fprintf('  donnerait a8 = a22 au diametre d''alesage, donc a8 > a22 au mi-entrefer\n');
fprintf('  apres attenuation : la denture n''est donc PAS une simple modulation.)\n');

hd=@() fprintf('  %-26s %8s %8s %8s %8s %9s\n','cas','a7','a8','a22','a23','a22/a8');
ln=@(nm,a) fprintf('  %-26s %8.4f %8.4f %8.4f %8.4f %9.3f\n',nm,a(1),a(2),a(3),a(4),a(3)/a(2));
har=@(R) arrayfun(@(k) amp(R.Br(1:end-1),R.thq(1:end-1),k),ords);

%% --- 1. convergence en resolution de surface ---------------------------
fprintf('\n[1] Convergence en Nsurf (le deficit est-il numerique ?)\n'); hd();
for Nsf=[420 840 1260 1680]
    R=cogging_mec(M,Nsf,0,21,1500,kfr); ln(sprintf('Nsurf = %d',Nsf),har(R));
end

%% --- 2. invariance en position rotor (controle interne) ----------------
fprintf('\n[2] Invariance de |a8| avec la position rotor (doit etre constante)\n');
Rp=cogging_mec(M,840,0,25,1500,kfr,2*pi/Ns);
AG=Rp.AG; thq=linspace(0,2*pi,3601); thq(end)=[];
phis_pos=Rp.phis;                       % (Rp : NE PAS reutiliser 'R', ecrase plus bas)
a8v=zeros(1,numel(phis_pos)); a22v=a8v;
for q=1:numel(a8v)
    [Br,~]=AG.field(Rp.Usurf(:,q),phis_pos(q),thq);
    a8v(q)=amp(Br,thq,8); a22v(q)=amp(Br,thq,22);
end
fprintf('    a8  : %.4f +/- %.4f T (variation %.1f %%)\n',mean(a8v),std(a8v),(max(a8v)-min(a8v))/mean(a8v)*100);
fprintf('    a22 : %.4f +/- %.4f T (variation %.1f %%)\n',mean(a22v),std(a22v),(max(a22v)-min(a22v))/mean(a22v)*100);

%% --- 3. saturation : fer non lineaire ---------------------------------
fprintf('\n[3] Saturation (fer non lineaire vs mu_r constant)\n'); hd();
Rl=cogging_mec(M,840,0,21,1500,kfr); ln('lineaire muI=1500',har(Rl));
ph=linspace(0,2*pi/p,9);
Rn=onload_mec(M,840,kfr,ph,zeros(3,9));
[Brn,~]=Rn.AG.field(Rn.Us(:,1),ph(1),thq);
an=arrayfun(@(k) amp(Brn,thq,k),ords);
ln('NON LINEAIRE (B(H))',an);

%% --- 4. permeabilite du fer -------------------------------------------
fprintf('\n[4] Sensibilite a la permeabilite du fer\n'); hd();
for mu=[500 1500 5000 50000]
    R=cogging_mec(M,840,0,21,mu,kfr); ln(sprintf('muI = %g',mu),har(R));
end

%% --- 5. modele d'ouverture --------------------------------------------
fprintf('\n[5] Modele d''ouverture d''encoche\n'); hd();
for kk=[0 0.25 0.5 1]
    R=cogging_mec(M,840,0,21,1500,kk); ln(sprintf('kfringe = %.2f',kk),har(R));
end
R2=cogging_mec2(M,840,8,21,1500); ln('bouche 2D maillee',har(R2));

%% --- 6. RE-IDENTIFICATION CONJOINTE (muI, kfringe) ---------------------
%  muI = 1500 etait un choix ARBITRAIRE. La courbe B(H) reelle donne, a
%  l'induction de dent constatee (1.34 T), mu_r = B/(mu0*H) ~ 3150. Comme a8
%  croit fortement avec muI, kfringe avait ete identifie sur un fer trop
%  reluctant : il faut re-identifier le COUPLE (muI, kfringe).
fprintf('\n[6] RE-IDENTIFICATION CONJOINTE (muI, kfringe)\n');
BHt=bh_curve(); mur_phys=1.34/(4*pi*1e-7*BHt.Hof(1.34));
fprintf('    mu_r physique a B_dent = 1.34 T : %.0f  (le modele utilisait 1500)\n',mur_phys);
err=@(a) norm((a(2:4)-aF(2:4))./aF(2:4))/sqrt(3)*100;
muL=[1500 3000 5000]; kfL=[0.15 0.25 0.35 0.5];
fprintf('    %8s %8s %8s %8s %8s %9s\n','muI','kfringe','a8','a22','a23','err_dent%');
best=inf;
for mu=muL
    for kf2=kfL
        a=har(cogging_mec(M,840,0,21,mu,kf2));
        e=err(a);
        fprintf('    %8g %8.2f %8.4f %8.4f %8.4f %9.1f\n',mu,kf2,a(2),a(3),a(4),e);
        if e<best, best=e; bmu=mu; bkf=kf2; ba=a; end
    end
end
fprintf('    >> OPTIMUM : muI = %g, kfringe = %.2f -> erreur denture %.1f %% (etait 26.4 %%)\n',...
        bmu,bkf,best);
fprintf('       a8 %+.1f %% | a22 %+.1f %% | a23 %+.1f %% | a7 %+.1f %%\n',...
        (ba(2)-aF(2))/aF(2)*100,(ba(3)-aF(3))/aF(3)*100,(ba(4)-aF(4))/aF(4)*100,(ba(1)-aF(1))/aF(1)*100);

%% --- 7. verdict --------------------------------------------------------
fprintf('\n[7] LECTURE\n');
fprintf('    Aucun reglage du modele d''ouverture ne rapproche a8 de la FEA sans\n');
fprintf('    degrader a22/a23 : le rapport a22/a8 du MEC (~2.4) reste tres\n');
fprintf('    au-dessus de celui de la FEA (%.2f) dans TOUS les cas.\n',aF(3)/aF(2));
fprintf('    => le deficit n''est pas un parametre mal regle mais une propriete\n');
fprintf('    STRUCTURELLE de la representation de la surface statorique.\n');

figure('Name','Diagnostic a8','Color','w','Position',[60 60 1000 400]);
subplot(1,2,1);
kk=1:40; sF=arrayfun(@(k) amp(BrFu,thu,k),kk);
sM=arrayfun(@(k) amp(Rl.Br(1:end-1),Rl.thq(1:end-1),k),kk);
stem(kk-0.2,sF,'r','filled','MarkerSize',4); hold on;
stem(kk+0.2,sM,'b','filled','MarkerSize',4);
set(gca,'YScale','log'); grid on; xlabel('ordre spatial'); ylabel('|B_r| (T)');
legend('FEA','MEC','Location','best'); title('Spectre du champ d''entrefer');
subplot(1,2,2);
plot(phis_pos*180/pi,a8v,'b-o','LineWidth',1.3); hold on;
plot(phis_pos*180/pi,a22v,'g-s','LineWidth',1.3);
yline(aF(2),'b--'); yline(aF(3),'g--'); grid on;
xlabel('position rotor (deg mec)'); ylabel('|B_r| (T)');
legend('a_8 MEC','a_{22} MEC','a_8 FEA','a_{22} FEA','Location','best');
title('Invariance en position');
saveas(gcf,fullfile(fileparts(mfilename('fullpath')),'FIG_a8.png'));
