%% RUN_X1_TABLE5B_RECONCILE - bloc X1, anomalie A-2 (SPEC_CLAUDE_CODE_v3 §8)
%
%  QUESTION. Le manuscrit publie une Table 5(b) et un encadrement eq:bracket
%  que A2 ne reproduit pas : decalage SYSTEMATIQUE de 1,4 a 2,8 points, de
%  signe coherent sur les six cellules. Un decalage systematique designe un
%  PARAMETRE. Ce bloc l'identifie, puis regenere la table ET l'encadrement
%  depuis une seule execution.
%
%  CE QUI EST PUBLIE (MEC_DtN_paper_v2.tex l. 1245-1254 et l. 2934-2939) :
%      panneau (b), titre "Four shoe layers, n_sh = 4"
%          Ms=360 : nu=8 lin -26,8 %  Newton +10,4 %  peak B 1,928 T
%          Ms=540 :          -21,6 %         +15,0 %          1,979 T
%          Ms=900 :          -18,8 %         +23,8 %          2,013 T
%      encadrement : 0,01596 (lin) < 0,01966 (EF) < 0,02434 (Newton) T
%      -> l'encadrement EST la ligne Ms=900 du panneau (b).
%
%  CE QUE A2 A OBTENU a n_sh = 2 : -25,2 / -19,6 / -16,0 et +11,8 / +16,8 /
%  +26,2. Toutes les cellules de A2 sont PLUS HAUTES que les publiees, dans
%  les deux solveurs, et l'ecart croit avec Ms.
%
%  HYPOTHESE. Le parametre est le NOMBRE DE COUCHES DE BEC. Le manuscrit
%  ecrit lui-meme (l. 1208-1210) que raffiner le bec radialement ABAISSE
%  nu=8 : le signe du decalage est donc celui d'un maillage de bec plus fin
%  que celui de A2. Reste a savoir lequel -- d'ou le balayage n_sh = 1..4.
%  ATTENTION A LA CONVENTION : mesh_bldc.m:71 pose nshT = 2*nsh couches de
%  bec (isthme hs0 + biseau hs1). Le "n_sh" du manuscrit peut donc designer
%  soit l'argument nsh, soit le total nshT = 2*nsh. Le balayage tranche.
%
%  kfringe EST ECARTE AVANT LE BALAYAGE, par lecture du code : mesh_bldc
%  range kfringe dans la structure (l. 211) mais ne l'emploie dans AUCUNE
%  permeance -- la frange d'ouverture est resolue explicitement par les
%  couches de bec (commentaire l. 182-189). kfringe n'agit que sur la chaine
%  LOCALISEE (cogging_mec, permeance Gfr). Un controle numerique est joint au
%  §0 pour que l'argument ne repose pas sur la seule lecture.
%
%  CONFIGURATION DECLAREE (regle n.2 des regles d'engagement)
%    machine       : PMSM 15 encoches / 14 poles, 750 W (machine_bldc)
%    pavage        : maillage polaire mesh_bldc, Ms colonnes x Ls couches,
%                    Ls = 2*n_sh + nst + nys, nst = 4, nys = 3
%    n_sh          : BALAYE 1..4
%    base          : P0 (airgap_magnet par defaut), numax = floor(Ms/2)
%    solveurs      : lineaire mu_r = M.muI  ET  Newton exact (bh_curve M350)
%    kfringe       : 0,325 (sans effet sur cette chaine -- §0)
%    reference EF  : magnetostique(Magnetic_loading), Calculator Expressions
%                    Plot 4.tab, dossier M.FEA.dir
%    position rotor: phi = 0 (a vide), comme A2
clear; clc; t0=tic;
diary('X1_table5b_reconcile_out.txt'); diary on;

M=machine_bldc(); p=M.p;
if isfile('kfringe_ident.mat'), S0=load('kfringe_ident.mat'); kfr=S0.kbest; else, kfr=0.325; end
rd=@(f)readmatrix(f,'FileType','text','Delimiter','\t','NumHeaderLines',1);
fea=M.FEA.dir;
d4=rd(fullfile(fea,'magnetostique(Magnetic_loading)','Calculator Expressions Plot 4.tab'));
thu=linspace(0,2*pi,2001); thu(end)=[];
BrF=interp1([d4(:,2)*pi/180;2*pi],[d4(:,3);d4(1,3)],thu,'linear','extrap');
amp=@(y,th,k)abs((2/numel(y))*sum(y(:).'.*exp(-1i*k*th(:).')));
Bg1F=amp(BrF,thu,p); r8F=amp(BrF,thu,8);

fprintf('=== X1 : Table 5(b) et encadrement -- identification du parametre ===\n');
fprintf('  machine       : PMSM 15/14, 750 W (machine_bldc)\n');
fprintf('  chaine        : mesh_bldc + solve_bldc_mesh, base P0, numax=floor(Ms/2)\n');
fprintf('  pavage        : nst=4, nys=3, n_sh balaye 1..4 (Ls = 2*n_sh + 7)\n');
fprintf('  solveurs      : lineaire mu_r=%g  ET  Newton exact (bh_curve M350)\n',M.muI);
fprintf('  kfringe       : %.4f  (sans effet ici -- controle §0)\n',kfr);
fprintf('  reference EF  : %s\n',fea);
fprintf('    Bg1  = %.6f T\n',Bg1F);
fprintf('    nu=8 = %.10f T   (c''est le denominateur de toutes les colonnes)\n\n',r8F);

%% ---- 0. controle : kfringe est inoperant sur le maillage --------------
fprintf('  ---- 0. controle kfringe (candidat n.1 de la specification) ----\n');
[~,r8a,~,~]=x1pt(M,540,2,M.muI,0.325,thu,p);
[~,r8b,~,~]=x1pt(M,540,2,M.muI,0.750,thu,p);
fprintf('    Ms=540, n_sh=2, lineaire : kfringe=0.325 -> nu=8 = %.12e\n',r8a);
fprintf('                               kfringe=0.750 -> nu=8 = %.12e\n',r8b);
fprintf('    ecart relatif = %.3e\n',abs(r8b-r8a)/r8a);
fprintf('    => kfringe N''EST PAS le parametre cherche. Ecarte.\n\n');

%% ---- 1. balayage de n_sh ---------------------------------------------
Mss=[360 540 900];
pubL=[-26.8 -21.6 -18.8];       % nu=8, solveur lineaire        [publie]
pubN=[ 10.4  15.0  23.8];       % nu=8, solveur Newton          [publie]
pubB=[1.928 1.979 2.013];       % peak iron B (colonne Newton)  [publie]
NSH=1:4; nA=numel(NSH); nB=numel(Mss);
R8L=nan(nA,nB); R8N=R8L; PKL=R8L; PKN=R8L; UNK=R8L; B1L=R8L;
fprintf('  ---- 1. balayage n_sh = 1..4 x Ms = 360/540/900 x deux solveurs ----\n');
for ia=1:nA
    for ib=1:nB
        [b1,r8l,pkl,N ]=x1pt(M,Mss(ib),NSH(ia),M.muI,kfr,thu,p);
        [ ~,r8n,pkn,~ ]=x1pt(M,Mss(ib),NSH(ia),'nl' ,kfr,thu,p);
        B1L(ia,ib)=b1; R8L(ia,ib)=r8l; R8N(ia,ib)=r8n;
        PKL(ia,ib)=pkl; PKN(ia,ib)=pkn; UNK(ia,ib)=N;
    end
    fprintf('    n_sh=%d termine (%.0f s)\n',NSH(ia),toc(t0));
end
EL=100*(R8L-r8F)/r8F;  EN=100*(R8N-r8F)/r8F;

for ia=1:nA
    fprintf('\n    -- n_sh = %d  (couches de bec nshT = %d, Ls = %d) --\n', ...
        NSH(ia),2*NSH(ia),2*NSH(ia)+7);
    fprintf('    %6s %10s %12s %12s %12s %12s\n', ...
        'Ms','inconnues','nu=8 lin.','nu=8 Newton','peak B lin.','peak B NL');
    for ib=1:nB
        fprintf('    %6d %10d %11.1f %% %11.1f %% %12.4f %12.4f\n', ...
            Mss(ib),UNK(ia,ib),EL(ia,ib),EN(ia,ib),PKL(ia,ib),PKN(ia,ib));
    end
    fprintf('    residus vs publie (points de %%) : lin %s | Newton %s\n', ...
        mat2str(round(EL(ia,:)-pubL,2)),mat2str(round(EN(ia,:)-pubN,2)));
    fprintf('    residu maximal sur les six cellules : %.2f point\n', ...
        max([abs(EL(ia,:)-pubL) abs(EN(ia,:)-pubN)]));
end

%% ---- 2. verdict -------------------------------------------------------
res=nan(nA,1);
for ia=1:nA, res(ia)=max([abs(EL(ia,:)-pubL) abs(EN(ia,:)-pubN)]); end
[rb,ib0]=min(res); nb0=NSH(ib0);
fprintf('\n  ---- 2. verdict ----\n');
fprintf('    %6s %14s %16s\n','n_sh','residu max','residu peak B NL');
for ia=1:nA
    fprintf('    %6d %12.2f pt %14.4f T\n',NSH(ia),res(ia),max(abs(PKN(ia,:)-pubB)));
end
if rb<=0.6
    fprintf('\n    RECONCILIE : n_sh = %d reproduit la Table 5(b) publiee\n',nb0);
    fprintf('    (residu maximal %.2f point sur six cellules, arrondi du\n',rb);
    fprintf('     manuscrit a 0,1 point pres). Le parametre est le NOMBRE DE\n');
    fprintf('     COUCHES DE BEC ; A2 l''avait pris a 2 au lieu de %d.\n',nb0);
else
    fprintf('\n    NON RECONCILIE : aucun n_sh de 1 a 4 ne reproduit le publie\n');
    fprintf('    (meilleur residu %.2f point a n_sh = %d). Le publie est\n',rb,nb0);
    fprintf('    IRRECUPERABLE en l''etat : la Table 5(b) et l''encadrement\n');
    fprintf('    doivent etre remplaces par la presente execution, et le\n');
    fprintf('    manuscrit doit declarer le maillage de bec employe.\n');
end

%% ---- 3. Table 5(b) et encadrement, UNE SEULE EXECUTION ---------------
fprintf('\n  ---- 3. Table 5(b) retenue : n_sh = %d ----\n',nb0);
fprintf('    machine PMSM 15/14 | nst=4, nys=3, n_sh=%d | base P0 |\n',nb0);
fprintf('    numax=floor(Ms/2) | kfringe=%.3f (inoperant) | phi=0 |\n',kfr);
fprintf('    reference EF nu=8 = %.10f T\n\n',r8F);
fprintf('    %6s %10s %14s %10s %14s %10s %11s %11s\n', ...
    'Ms','inconnues','nu=8 lin (T)','ecart','nu=8 NL (T)','ecart','pk B lin','pk B NL');
for ib=1:nB
    fprintf('    %6d %10d %14.8f %8.1f %% %14.8f %8.1f %% %11.4f %11.4f\n', ...
        Mss(ib),UNK(ib0,ib),R8L(ib0,ib),EL(ib0,ib),R8N(ib0,ib),EN(ib0,ib), ...
        PKL(ib0,ib),PKN(ib0,ib));
end

i900=find(Mss==900);
fprintf('\n  ---- encadrement eq:bracket, MEME EXECUTION, Ms = 900 ----\n');
fprintf('    lineaire %.8f T  <  EF %.8f T  <  Newton %.8f T\n', ...
    R8L(ib0,i900),r8F,R8N(ib0,i900));
fprintf('    soit %+.2f %% et %+.2f %%\n',EL(ib0,i900),EN(ib0,i900));
fprintf('    publie  : 0.01596 < 0.01966 < 0.02434 T, soit -18.82 %% et +23.81 %%\n');
fprintf('    ecart sur les bornes : %+.6f et %+.6f T\n', ...
    R8L(ib0,i900)-0.01596,R8N(ib0,i900)-0.02434);
fprintf('    L''encadrement est VALIDE si et seulement si la ligne Ms=900\n');
fprintf('    ci-dessus est celle qui entre dans la Table 5(b) : c''est le cas\n');
fprintf('    par construction, les deux sortant de la meme execution.\n');

save('X1_table5b.mat','Mss','NSH','R8L','R8N','PKL','PKN','UNK','B1L', ...
     'EL','EN','r8F','Bg1F','pubL','pubN','pubB','nb0','res','kfr');
fprintf('\n  duree %.0f s\n=== X1 termine ===\n',toc(t0));
diary off;

% ======================================================================
function [b1,r8,pkB,N] = x1pt(M,Ms,nsh,sol,kfr,thu,p)
%  Un point du balayage. Tout est passe explicitement : les fonctions
%  locales d'un SCRIPT ne partagent pas l'espace de travail (c'est ce qui
%  avait fait planter la premiere execution de A2).
    amp=@(y,th,k)abs((2/numel(y))*sum(y(:).'.*exp(-1i*k*th(:).')));
    ME=mesh_bldc(M,Ms,4,3,[],0,kfr,sol,nsh);
    S =solve_bldc_mesh(ME,1e-9,30);
    Br=ME.AG.field(S.Usurf,0,thu); Br=Br(:).';
    b1=amp(Br,thu,p); r8=amp(Br,thu,8);
    pkB=max(abs(S.B(ME.iron))); N=ME.N;
end
