%% RUN_G1_GAMMA - convention sur gamma : regenerer la croissance du crochet
%
%  QUESTION (message v7, section 2.2). Le manuscrit porte DEUX quantites sous
%  le meme intitule "bracket" :
%      Table 6, panneau bas : ln N + ln|2 sin(d/2)|          -> ln(pi)
%      Table 7, panneau (a) : ln N + gamma + ln|2 sin(d/2)|  -> ln(pi)+gamma
%  L'equation (28) est la forme exacte : le developpement de
%  sum_{n<=N} (1/n) sin^2(nd/2) = 1/2 [ ln N + gamma + ln|2 sin(d/2)| ] + o(1)
%  contient gamma. Retenir gamma partout fait passer la croissance annoncee
%  de +182 % a une autre valeur, qui DOIT etre regeneree ici et non formee
%  par arithmetique sur les valeurs imprimees (regle 2 du dossier).
%
%  CE QUE CE BLOC PRODUIT
%    1. le balayage a troncature relachee, M_s = 1080 FIXE, N_h libre,
%       chaine identique a RUN_R4_UNLOCK (cogging_mec + airgap_magnet, P0) ;
%    2. le crochet SANS gamma et AVEC gamma, a pleine precision ;
%    3. les deux croissances, formees sur les valeurs pleine precision et
%       arrondies une seule fois a l'affichage ;
%    4. les deux pentes par decade, contre ln 10 ;
%    5. le panneau LaTeX de la Table 6 sous la convention AVEC gamma.
%
%  GARDE. gamma est une constante additive : la PENTE par decade doit etre
%  IDENTIQUE dans les deux conventions, au bit pres. Si les deux pentes
%  different, le calcul est faux et le resultat doit etre rejete.
%  GARDE 2. A N_h = M_s/2 exactement, la chaine deverrouillee doit redonner
%  la chaine verrouillee au bit pres (garde de R4, rejouee).
clear; clc; t0=tic;
diary('G1_gamma_out.txt'); diary on;

gam = 0.57721566490153286;      % constante d'Euler-Mascheroni

M=machine_bldc(); p=M.p;
if isfile('kfringe_ident.mat'), S0=load('kfringe_ident.mat'); kfr=S0.kbest; else, kfr=0.325; end
thu=linspace(0,2*pi,2001); thu(end)=[];
amp=@(y,th,k)abs((2/numel(y))*sum(y(:).'.*exp(-1i*k*th(:).')));

fprintf('=== G1 : convention sur gamma, croissance du crochet regeneree ===\n\n');
fprintf('  CONFIGURATION\n');
fprintf('    machine     : PMSM 15/14, 750 W (machine_bldc)\n');
fprintf('    chaine      : cogging_mec + airgap_magnet (identique a R4)\n');
fprintf('    base        : P0 (constante par morceaux)\n');
fprintf('    solveur     : lineaire, mu_i = %g\n',M.muI);
fprintf('    k_fringe    : %.4f\n',kfr);
fprintf('    pavage      : M_s = 1080 FIXE\n');
fprintf('    gamma       : %.17f\n',gam);

Msf=1080; nu_lock=floor(Msf/2); d=2*pi/Msf;
nul=nu_lock*[1 2 4 8];          % le balayage PUBLIE (Table 6, panneau bas)

fprintf('\n  ---- 1. BALAYAGE A TRONCATURE RELACHEE (le balayage publie) ----\n');
fprintf('  %8s %16s %20s %14s\n','N_h','ln N + ln|2sin|','ln N + g + ln|2sin|','Bg1 (T)');
W=nan(numel(nul),4);
for i=1:numel(nul)
    N=nul(i);
    R=cogging_mec(M,Msf,N,3,M.muI,kfr,1e-6);
    Br=R.AG.field(R.U(1:Msf,1),0,thu); Br=Br(:).';
    br  = log(R.numax_eff)+log(abs(2*sin(d/2)));
    brg = br + gam;
    W(i,:)=[R.numax_eff br brg amp(Br,thu,p)];
    fprintf('  %8d %16.12f %20.12f %14.6f\n',W(i,1),br,brg,W(i,4));
end

fprintf('\n  ---- 2. CROISSANCE SUR LE BALAYAGE, PLEINE PRECISION ----\n');
gs = 100*(W(end,2)-W(1,2))/W(1,2);      % sans gamma
gg = 100*(W(end,3)-W(1,3))/W(1,3);      % avec gamma
fprintf('    SANS gamma : %.12f -> %.12f\n',W(1,2),W(end,2));
fprintf('                 croissance = %.10f %%   -> arrondi : %+.0f %%\n',gs,gs);
fprintf('    AVEC gamma : %.12f -> %.12f\n',W(1,3),W(end,3));
fprintf('                 croissance = %.10f %%   -> arrondi : %+.1f %%\n',gg,gg);

fprintf('\n  ---- 3. PENTE PAR DECADE (la garde) ----\n');
ps=polyfit(log10(W(:,1)),W(:,2),1);
pg=polyfit(log10(W(:,1)),W(:,3),1);
fprintf('    pente SANS gamma : %.12f par decade\n',ps(1));
fprintf('    pente AVEC gamma : %.12f par decade\n',pg(1));
fprintf('    ln 10            : %.12f\n',log(10));
fprintf('    ecart SANS gamma a ln10 : %+.4f %%\n',100*(ps(1)-log(10))/log(10));
fprintf('    ecart AVEC gamma a ln10 : %+.4f %%\n',100*(pg(1)-log(10))/log(10));
dp=abs(ps(1)-pg(1));
fprintf('    |pente_sans - pente_avec| = %.3e\n',dp);
if dp<1e-12
    fprintf('    GARDE 1 PASSEE : gamma est additif, la pente ne bouge pas.\n');
else
    fprintf('    GARDE 1 ECHOUEE : les deux pentes different, calcul a rejeter.\n');
end

fprintf('\n  ---- 4. GARDE 2 : verrou contre deverrouillage a N_h = M_s/2 ----\n');
Rl=cogging_mec(M,Msf,0,3,M.muI,kfr,1e-6);
Ru=cogging_mec(M,Msf,nu_lock,3,M.muI,kfr,1e-6);
Bl=Rl.AG.field(Rl.U(1:Msf,1),0,thu); Bl=Bl(:).';
Bu=Ru.AG.field(Ru.U(1:Msf,1),0,thu); Bu=Bu(:).';
fprintf('    Bg1 verrouille   : %.15f T\n',amp(Bl,thu,p));
fprintf('    Bg1 deverrouille : %.15f T\n',amp(Bu,thu,p));
fprintf('    ecart max sur Br : %.3e T\n',max(abs(Bl-Bu)));
if max(abs(Bl-Bu))<1e-12
    fprintf('    GARDE 2 PASSEE : identite au bit pres.\n');
else
    fprintf('    GARDE 2 ECHOUEE.\n');
end

fprintf('\n  ---- 5. VALEUR LIMITE DU CROCHET SOUS VERROU ----\n');
fprintf('    ln(pi)         = %.12f\n',log(pi));
fprintf('    ln(pi) + gamma = %.12f\n',log(pi)+gam);

fprintf('\n  ---- 6. PANNEAU LaTeX, Table 6 panneau bas, AVEC gamma ----\n\n');
fprintf('$N_h$ & \\multicolumn{2}{c}{$\\ln N+\\gamma+\\ln|2\\sin(d/2)|$}\n');
fprintf('      & \\multicolumn{2}{c}{$B_{g1}$ (\\si{\\tesla}), dev.}\\\\\n');
for i=1:numel(nul)
    if i==1
        fprintf('\\num{%d}  & \\multicolumn{2}{c}{\\num{%.4f}} & \\multicolumn{2}{c}{\\num{%.5f}\\quad ---}\\\\\n', ...
            W(i,1),W(i,3),W(i,4));
    else
        fprintf('\\num{%d} & \\multicolumn{2}{c}{\\num{%.4f}} & \\multicolumn{2}{c}{\\num{%.5f}\\quad $%+.2f\\%%$}\\\\\n', ...
            W(i,1),W(i,3),W(i,4),100*(W(i,4)-W(1,4))/W(1,4));
    end
end

save('G1_gamma.mat','W','gam','gs','gg','ps','pg');
fprintf('\n  duree %.0f s\n=== G1 termine ===\n',toc(t0));
diary off;
