%% RUN_DIAG_COGSUM - d'ou vient le facteur 11 sur la detente du maillage ?
%
%  Le couple est calcule par le tenseur de Maxwell sous forme de Parseval :
%      T = (L r^2/mu0) pi SUM_n (Brc_n Btc_n + Brs_n Bts_n)
%  C'est une SOMME DE PRODUITS D'HARMONIQUES qui se compensent : la detente
%  vaut 0.6 mN.m alors que chaque terme individuel peut peser plusieurs
%  N.m. Le maillage fournit Ms/2 = 450 harmoniques de surface contre 630
%  pour le reseau a une dent, mais surtout une distribution de potentiel
%  BEAUCOUP plus structuree (le bec est resolu). Si les rangs eleves ne se
%  compensent plus, la somme diverge.
%
%  TEST. Sommes PARTIELLES T_N = SUM_{n<=N} en fonction de N, a plusieurs
%  positions. Une somme saine converge et reste plate ; une somme
%  contaminee croit ou oscille avec N.
clear; clc;
M=machine_bldc(); Nm=M.Nm; L=M.ls; mu0=4*pi*1e-7;
if isfile('kfringe_ident.mat'), S0=load('kfringe_ident.mat'); kfr=S0.kbest; else, kfr=0.325; end
phis=linspace(0,2*pi/Nm,841); qs=[100 200 300 400];
Nlist=[20 40 60 90 120 180 260 360 450];

fprintf('=== Sommes partielles du tenseur de Maxwell ===\n');
fprintf('  MAILLAGE (Ms=900, Newton)\n');
fprintf('  %6s',"N<="); fprintf(' %9d',Nlist); fprintf('\n');
for q=qs
    ME=mesh_bldc(M,900,4,3,[],phis(q),kfr,'nl',2);
    S=solve_bldc_mesh(ME,1e-9,40);
    T=partsum(ME.AG,S.Usurf,phis(q),M,L,mu0,Nlist);
    fprintf('  q=%4d',q); fprintf(' %9.3f',T*1e3); fprintf('  mN.m\n');
end

%% ---- meme test sur le reseau a une dent -------------------------------
fprintf('\n  RESEAU A UNE DENT (Nsurf=1260)\n');
R1=cogging_mec(M,1260,0,841,M.muI,kfr);
fprintf('  %6s',"N<="); fprintf(' %9d',Nlist); fprintf('\n');
for q=qs
    T=partsum(R1.AG,R1.U(1:1260,q),R1.phis(q),M,L,mu0,Nlist);
    fprintf('  q=%4d',q); fprintf(' %9.3f',T*1e3); fprintf('  mN.m\n');
end

%% ---- poids des harmoniques : ou se joue la compensation ? -------------
ME=mesh_bldc(M,900,4,3,[],phis(200),kfr,'nl',2);
S=solve_bldc_mesh(ME,1e-9,40);
[tm,cum]=terms(ME.AG,S.Usurf,phis(200),M,L,mu0);
[t1,cu1]=terms(R1.AG,R1.U(1:1260,200),R1.phis(200),M,L,mu0);
fprintf('\n  --- amplitude des TERMES individuels (position q=200) ---\n');
fprintf('  %-18s %12s %12s\n','','maillage','1 dent');
fprintf('  %-18s %12.3f %12.3f  N.m\n','terme max',max(abs(tm)),max(abs(t1)));
fprintf('  %-18s %12.3f %12.3f  N.m\n','somme des |termes|',sum(abs(tm)),sum(abs(t1)));
fprintf('  %-18s %12.4f %12.4f  N.m\n','somme algebrique',sum(tm),sum(t1));
fprintf('  %-18s %12.1e %12.1e\n','conditionnement', ...
    sum(abs(tm))/max(abs(sum(tm)),eps),sum(abs(t1))/max(abs(sum(t1)),eps));
fprintf('\n  Le conditionnement est le rapport somme des valeurs absolues sur\n');
fprintf('  somme algebrique : c''est le facteur d''amplification de toute erreur\n');
fprintf('  relative sur un terme. Au-dela de 1e4, la somme n''a plus de sens.\n');

%% ------------------------------------------------------------------------
function T=partsum(AG,Us,phq,M,L,mu0,Nlist)
    t=terms(AG,Us,phq,M,L,mu0);
    T=zeros(1,numel(Nlist));
    for k=1:numel(Nlist), T(k)=sum(t(1:min(Nlist(k),numel(t)))); end
end
function [t,cum]=terms(AG,Us,phq,M,L,mu0)
    nu=AG.nu; Usc=AG.Wc*Us; Uss=AG.Ws*Us;
    cq=cos(nu*phq); sq=sin(nu*phq);
    Brc=AG.bru.*Usc+AG.brmq.*cq;  Brs=AG.bru.*Uss+AG.brmq.*sq;
    Btc=-AG.btu.*Uss-AG.btmq.*sq; Bts=AG.btu.*Usc+AG.btmq.*cq;
    t=(L*M.rmid^2/mu0)*pi*(Brc.*Btc+Brs.*Bts);
    cum=cumsum(t);
end
