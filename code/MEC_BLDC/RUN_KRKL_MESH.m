%% RUN_KRKL_MESH - k_r et k_l sur le maillage raffine
%
%  k_r et k_l ne dependent QUE du champ dans la couronne aimant + entrefer,
%  donc des POTENTIELS DE SURFACE. Le maillage fournit exactement la meme
%  grandeur que le reseau a une dent : on reprend telle quelle l'extraction
%  de la section 2 de BLDC_MEC_COMPLET, en lui donnant S.Usurf.
%    k_r = SUM(mmf_m)/SUM(mmf_g) le long d'une ligne radiale au centre de
%          chaque aimant  [definition exacte lue dans BLDC.aedt]
%    k_l = phi_ag / SUM(phi_Bm), flux mesure DANS l'aimant a mi-hauteur
clear; clc; t0=tic;
M=machine_bldc(); Ns=M.Ns; Nm=M.Nm; p=M.p; L=M.ls; mu0=4*pi*1e-7;
if isfile('kfringe_ident.mat'), S0=load('kfringe_ident.mat'); kfr=S0.kbest; else, kfr=0.325; end
rd=@(f)readmatrix(f,'FileType','text','Delimiter','\t','NumHeaderLines',1);
fea=M.FEA.dir;
dmm=rd(fullfile(fea,'magnetostique(Magnetic_loading)','Calculator Expressions Table 1.tab'));
mmf_m_F=dmm(1,2:15); mmf_g_F=dmm(1,16:29);
kr_F=rd(fullfile(fea,'magnetostique(Magnetic_loading)','Output Variables Table 2.tab')); kr_F=kr_F(1,2);
kl_F=rd(fullfile(fea,'magnetostique(Magnetic_loading)','Output Variables Table 1.tab')); kl_F=kl_F(1,2);
kr_F_clean=mean(mmf_m_F)/mean(mmf_g_F(2:end));

fprintf('=== k_r et k_l : maillage raffine vs reseau a une dent vs FEA ===\n');
fprintf('  FEA : k_r = %.4f (%.4f hors point aberrant) | k_l = %.4f\n', ...
    kr_F,kr_F_clean,kl_F);
fprintf('        FMM aimant %.1f A | FMM entrefer %.1f A\n\n', ...
    mean(mmf_m_F),mean(mmf_g_F(2:end)));

fprintf('  %-22s %9s %9s %9s %9s %9s %9s\n', ...
    'modele','FMM_m','FMM_g','k_r','k_l','disp. %','ec. k_l');
%  --- maillage, plusieurs resolutions ---
for Ms=[360 540]
    for muI={M.muI,'nl'}
        ME=mesh_bldc(M,Ms,4,3,[],0,kfr,muI{1},2);
        S=solve_bldc_mesh(ME,1e-9,30);
        [mm,mg,kr,kl,dm]=krkl(ME.AG,S.Usurf,0,M,L,mu0);
        lab=sprintf('maillage Ms=%d %s',Ms,ternary(ischar(muI{1}),'NL','lin'));
        fprintf('  %-22s %9.1f %9.1f %9.4f %9.4f %9.2f %8.1f %%\n', ...
            lab,mm,mg,kr,kl,dm,100*(kl-kl_F)/kl_F);
    end
end
%  --- reseau a une dent (reference interne) ---
R1=cogging_mec(M,1260,0,181,M.muI,kfr,2*pi/p);
[mm1,mg1,kr1,kl1,dm1]=krkl(R1.AG,R1.U(1:1260,1),0,M,L,mu0);
fprintf('  %-22s %9.1f %9.1f %9.4f %9.4f %9.2f %8.1f %%\n', ...
    'reseau a une dent',mm1,mg1,kr1,kl1,dm1,100*(kl1-kl_F)/kl_F);
%  --- dispersion de la reference, meme definition ---
fprintf('\n  dispersion FEA (meme definition) : %.2f %%\n', ...
    100*std(mmf_m_F)/mean(mmf_m_F));
fprintf('\n  ecarts sur k_r (reference FEA %.4f, hors point aberrant %.4f)\n',kr_F,kr_F_clean);
fprintf('\n  (%.0f s)\n',toc(t0));

%% ------------------------------------------------------------------------
function [mmf_m,mmf_g,kr,kl,disp_m]=krkl(AG,Us,phq,M,L,mu0)
%  disp_m : DISPERSION pole a pole de la FMM d'aimant, ecart-type normalise
%  par la moyenne [%]. C'est la ligne "magnet-MMF dispersion" du tableau 7 ;
%  sur la reference elle vaut 6,58 % et mesure l'encochage vu du rotor, pas
%  une dispersion de fabrication.
%  Extraction identique a la section 2 de BLDC_MEC_COMPLET, alimentee par
%  les potentiels de surface fournis en entree.
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
    %  --- k_l : flux au bore / flux DANS l'aimant a mi-hauteur ---
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
function o=ternary(c,a,b), if c, o=a; else, o=b; end, end
