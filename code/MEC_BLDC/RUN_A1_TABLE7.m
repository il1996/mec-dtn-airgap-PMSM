%% RUN_A1_TABLE7 - Table 7 regeneree en UNE execution, config UNIQUE
%
%  A1. Aujourd'hui la colonne Mesh melange n_sh = 1 (champ, nu=8) et
%  n_sh = 2 (k_r, k_l, inconnues), et la Table 5(a) une troisieme
%  combinaison contre une reference L_d perimee. On regenere donc TOUT dans
%  une configuration unique et DECLAREE, les deux valeurs de n_sh etant
%  produites cote a cote pour que le choix soit documente et non subi.
%
%  BASE : P0. V1 a etabli que P1 ne deplace B_g1 que de 0,01 % sur le PMSM
%  (contre -5,8 % sur X_m cote MAS) : la base n'est pas un enjeu ici.
%
%  Aucune valeur n'est transcrite : toutes les references EF sont relues des
%  fichiers du projet.
clear; clc; t0=tic;
diary('A1_table7_out.txt'); diary on;
M=machine_bldc(); Ns=M.Ns; Nm=M.Nm; p=M.p; L=M.ls; Ntc=M.Ntc;
PA=[1 -2 -15 3 14]; PB=[6 -7 -5 8 4]; PC=[11 -12 -10 13 9];
if isfile('kfringe_ident.mat'), S0=load('kfringe_ident.mat'); kfr=S0.kbest; else, kfr=0.325; end
rd=@(f)readmatrix(f,'FileType','text','Delimiter','\t','NumHeaderLines',1);
fea=M.FEA.dir; mu0=4*pi*1e-7;
amp=@(B,th,k) 2*abs(mean(B(:).'.*exp(-1i*k*th(:).')));
Ms=540; Np=61; om=M.speed*2*pi/60;
thu=linspace(0,2*pi,3601); thu(end)=[];

fprintf('=== A1 : Table 7 en configuration unique ===\n');
fprintf('  machine 15/14 750 W | Ms = %d | Np = %d | base P0\n',Ms,Np);
fprintf('  kfringe = %.4f | solveur : lineaire (mu_i) ET Newton, declares\n',kfr);
fprintf('  references EF relues du projet, aucune transcription\n\n');

%% ---------- references EF ---------------------------------------------
d4=rd(fullfile(fea,'magnetostique(Magnetic_loading)','Calculator Expressions Plot 4.tab'));
d2=rd(fullfile(fea,'magnetostique(Magnetic_loading)','Calculator Expressions Plot 2.tab'));
dA=rd(fullfile(fea,'magnetostique(Magnetic_loading)','Plot 1.tab'));
BrF=interp1([d4(:,2)*pi/180;2*pi],[d4(:,3);d4(1,3)],thu,'linear','extrap');
BtF=interp1(linspace(0,2*pi,numel(d2(:,2))).',d2(:,2),thu,'linear','extrap');
AF =max(abs(dA(:,2)-mean(dA(:,2))));
dFL=rd(fullfile(fea,'transitoire (Back_emf)','Output Variables Plot 4.tab'));
dV =rd(fullfile(fea,'transitoire (Back_emf)','NodeVoltage Plot.tab'));
dE =rd(fullfile(fea,'transitoire (Back_emf)','Output Variables Plot.tab'));
FLaF=max(abs(dFL(:,3))); EaF=max(abs(dV(:,3))); envF=max(dE(:,3));
dL=rd(fullfile(fea,'magnetostique(Armature-Field)','Output Variables Table 1.tab'));
LaF=dL(1,2); MF=(dL(1,3)+dL(1,4))/2; LdF=LaF-MF;
krF=rd(fullfile(fea,'magnetostique(Magnetic_loading)','Output Variables Table 2.tab')); krF=krF(1,2);
klF=rd(fullfile(fea,'magnetostique(Magnetic_loading)','Output Variables Table 1.tab')); klF=klF(1,2);

%% ---------- balayage maillage, les DEUX n_sh --------------------------
Rm=struct();
for nsh=[1 2]
    PhiT=zeros(Ns,Np); Us=zeros(Ms,Np); U0=[];
    for q=1:Np
        phi=(q-1)*(2*pi/p)/(Np-1);
        ME=mesh_bldc(M,Ms,4,3,[],phi,kfr,M.muI,nsh);
        %  PAS DE DEMARRAGE A CHAUD ICI. Avec le solveur LINEAIRE (muI fixe),
        %  passer U0 renvoie la solution precedente inchangee : le flux
        %  totalise devient constant sur tout le balayage et la FEM sort a
        %  zero. Le demarrage a chaud n'est licite qu'avec le Newton ('nl'),
        %  ou il sert d'initialisation et non de solution acceptee.
        S=solve_bldc_mesh(ME,1e-9,30);
        Us(:,q)=S.Usurf;
        mT=ME.type==1 & ME.lay==ME.lamT & ME.iron;
        PhiT(:,q)=accumarray(ME.tooth(ME.col(mT)).',S.Phi(mT),[Ns 1]);
        if q==1, MEr=ME; Sr=S; end
    end
    phis=linspace(0,2*pi/p,Np);
    [Br,Bt]=MEr.AG.field(Us(:,1),phis(1),thu); Br=Br(:).'; Bt=Bt(:).';
    lamf=@(P) Ntc*sum(sign(P(:)).*PhiT(abs(P(:)),:),1);
    lam=[lamf(PA); lamf(PB); lamf(PC)];
    e=zeros(3,Np); for k=1:3, e(k,:)=-gradient(lam(k,:),phis)*om; end
    fprintf('  [diag n_sh=%d] lam_a : min %.5f max %.5f etendue %.3e | om %.2f\n', ...
        nsh,min(lam(1,:)),max(lam(1,:)),max(lam(1,:))-min(lam(1,:)),om);
    fprintf('  [diag n_sh=%d] e_a   : min %.3f max %.3f V\n',nsh,min(e(1,:)),max(e(1,:)));
    Ac=cumtrapz(thu,Br)*M.rmid; Ac=Ac-mean(Ac);
    [mm,mg,kr,kl,dm]=krkl(MEr.AG,Us(:,1),0,M,L,mu0);
    %  inductances : FMM de bobinage, aimants eteints
    MEa=mesh_bldc(M,Ms,4,3,[1;0;0],0,kfr,M.muI,nsh); MEa.Isrc(:)=0;
    Sa=solve_bldc_mesh(MEa);
    Fb=mesh_bldc(M,Ms,4,3,[0;1;0],0,kfr,M.muI,nsh).E;
    Fc=mesh_bldc(M,Ms,4,3,[0;0;1],0,kfr,M.muI,nsh).E;
    La=sum(MEa.E.*Sa.Phi); Mu=0.5*(sum(Fb.*Sa.Phi)+sum(Fc.*Sa.Phi));
    Rm(nsh).v=[amp(Br,thu,p) mean(abs(Br)) max(Br) sqrt(mean(Bt.^2)) ...
        amp(Br,thu,8) max(abs(Ac)) max(abs(lam(1,:))) mm mg dm ...
        max(abs(e(1,:))) max(max(e,[],1)-min(e,[],1)) ...
        La*1e3 Mu*1e3 (La-Mu)*1e3 kr kl];
    Rm(nsh).N=MEa.N;
end

%% ---------- reseau a une dent (colonne Lumped) ------------------------
R1=cogging_mec(M,1260,0,Np,M.muI,kfr,2*pi/p);
[Br1,Bt1]=R1.AG.field(R1.U(1:1260,1),0,thu); Br1=Br1(:).'; Bt1=Bt1(:).';
lam1=[Ntc*sum(sign(PA(:)).*R1.PhiT(abs(PA(:)),:),1); ...
      Ntc*sum(sign(PB(:)).*R1.PhiT(abs(PB(:)),:),1); ...
      Ntc*sum(sign(PC(:)).*R1.PhiT(abs(PC(:)),:),1)];
e1=zeros(3,Np); for k=1:3, e1(k,:)=-gradient(lam1(k,:),R1.phis)*om; end
A1c=cumtrapz(thu,Br1)*M.rmid; A1c=A1c-mean(A1c);
[m1,g1,kr1,kl1,d1]=krkl(R1.AG,R1.U(1:1260,1),0,M,L,mu0);
RI=inductance_mec(M,1260,M.muI,0);
VL=[amp(Br1,thu,p) mean(abs(Br1)) max(Br1) sqrt(mean(Bt1.^2)) ...
    amp(Br1,thu,8) max(abs(A1c)) max(abs(lam1(1,:))) m1 g1 d1 ...
    max(abs(e1(1,:))) max(max(e1,[],1)-min(e1,[],1)) ...
    RI.La*1e3 RI.M*1e3 RI.Ld*1e3 kr1 kl1];

%% ---------- reference EF ----------------------------------------------
%  FMM par aimant et par entrefer : 14 + 14 colonnes du meme fichier.
%  La 1re FMM d'entrefer est le point aberrant de la §5.4 (531 A contre
%  704-792 A ailleurs) : on l'exclut de la moyenne, comme la reference.
dmm=rd(fullfile(fea,'magnetostique(Magnetic_loading)','Calculator Expressions Table 1.tab'));
mmf_m_F=dmm(1,2:15); mmf_g_F=dmm(1,16:29);
dispF=100*std(mmf_m_F)/mean(mmf_m_F);
VF=[amp(BrF,thu,p) mean(abs(BrF)) max(BrF) sqrt(mean(BtF.^2)) ...
    amp(BrF,thu,8) AF FLaF mean(mmf_m_F) mean(mmf_g_F(2:end)) dispF EaF envF ...
    LaF*1e3 MF*1e3 LdF*1e3 krF klF];

%% ---------- tableau ----------------------------------------------------
lab={'Bg1 fondamental (T)','B entrefer moyen |Br| (T)','Br crete (T)', ...
     'Bt tangentiel rms (T)','bande de denture nu=8 (T)','ligne de flux A (Wb/m)', ...
     'flux totalise crete (Wb)','FMM moy./aimant (A)','FMM moy./entrefer (A)', ...
     'dispersion FMM (%)','FEM de phase crete (V)','enveloppe six-step (V)', ...
     'self La (mH)','mutuelle M (mH)','synchrone Ld (mH)', ...
     'facteur k_r (-)','facteur k_l (-)'};
fprintf('  inconnues : n_sh=1 -> %d | n_sh=2 -> %d | 1 dent -> 1260 surf.\n\n', ...
    Rm(1).N,Rm(2).N);
fprintf('  %-28s %12s %12s %12s %12s\n','GRANDEUR','Mesh n_sh=1','Mesh n_sh=2','Lumped','FEA');
for j=1:numel(lab)
    fprintf('  %-28s %12.5f %12.5f %12.5f %12.5f\n', ...
        lab{j},Rm(1).v(j),Rm(2).v(j),VL(j),VF(j));
end
fprintf('\n  --- ecarts a la FEA, formes sur les valeurs PLEINES ---\n');
fprintf('  %-28s %12s %12s %12s\n','GRANDEUR','Mesh n_sh=1','Mesh n_sh=2','Lumped');
for j=1:numel(lab)
    if isnan(VF(j)), fprintf('  %-28s %12s %12s %12s\n',lab{j},'--','--','--'); continue; end
    fprintf('  %-28s %11.3f %% %11.3f %% %11.3f %%\n',lab{j}, ...
        100*(Rm(1).v(j)-VF(j))/VF(j),100*(Rm(2).v(j)-VF(j))/VF(j),100*(VL(j)-VF(j))/VF(j));
end
save('A1_table7.mat','Rm','VL','VF','lab');
fprintf('\n  duree %.0f s\n=== A1 termine ===\n',toc(t0));
diary off;

% ======================================================================
function [mmf_m,mmf_g,kr,kl,disp_m]=krkl(AG,Us,phq,M,L,mu0)
    nu=AG.nu; Nm=M.Nm;
    Usc=AG.Wc*Us; Uss=AG.Ws*Us;
    pic=AG.alphau.*Usc+AG.betasrc.*cos(nu*phq);
    pis=AG.alphau.*Uss+AG.betasrc.*sin(nu*phq);
    thq=linspace(0,2*pi,2001); thq(end)=[]; dth_=thq(2)-thq(1);
    phi_i=(pic.'*cos(nu*thq)+pis.'*sin(nu*thq));
    phi_b=(Usc.'*cos(nu*thq)+Uss.'*sin(nu*thq));
    mm=zeros(1,Nm); mg=zeros(1,Nm);
    for k=1:Nm
        c=(k-1)*2*pi/Nm+phq;
        [~,j]=min(abs(angle(exp(1i*(thq-c)))));
        mm(k)=abs(phi_i(j)); mg(k)=abs(phi_i(j)-phi_b(j));
    end
    mmf_m=mean(mm); mmf_g=mean(mg); kr=mmf_m/mmf_g;
    disp_m=100*std(mm)/mmf_m;
    Brc_b=AG.g.*Usc+AG.brm.*cos(nu*phq);
    Brs_b=AG.g.*Uss+AG.brm.*sin(nu*phq);
    Br_bore=(Brc_b.'*cos(nu*thq)+Brs_b.'*sin(nu*thq));
    Phi_bore=sum(abs(Br_bore))*dth_*M.Rsi*L/2;
    Mn=AG.Mn; Pn=AG.Pn; Um=AG.Um; rmi=M.rmi; mur=M.mu_r;
    Mc=Mn.*cos(nu*phq); Msn=Mn.*sin(nu*phq);
    Pcv=Pn.*cos(nu*phq); Psv=Pn.*sin(nu*phq);
    Sm_=sinh(nu*Um); Cm_=cosh(nu*Um);
    ac =(pic-Pcv*M.Rro+Pcv*rmi.*Cm_)./Sm_;
    asv=(pis-Psv*M.Rro+Psv*rmi.*Cm_)./Sm_;
    r=0.5*(M.Rro+rmi); v=log(r/rmi);
    cmid=-mu0*mur*(Pcv+(ac .*nu.*cosh(nu*v)-Pcv*rmi.*nu.*sinh(nu*v))/r)+mu0*Mc;
    smid=-mu0*mur*(Psv+(asv.*nu.*cosh(nu*v)-Psv*rmi.*nu.*sinh(nu*v))/r)+mu0*Msn;
    Br_in=(cmid.'*cos(nu*thq)+smid.'*sin(nu*thq));
    Phi_mag=sum(abs(Br_in))*dth_*r*L/2;
    kl=Phi_bore/Phi_mag;
end
