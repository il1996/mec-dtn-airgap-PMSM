%% RUN_R4_UNLOCK - deverrouiller la troncature du pavage (§3.5, Article I)
%
%  BLOC R4 de SPEC_CLAUDE_CODE_v4. Le §3.5 enonce une prediction
%  FALSIFIABLE -- "Unlock the truncation from the tiling, and the drift
%  returns" -- qu'il ne teste pas. C'est la seule prise franche laissee sur
%  l'Article I. Ce bloc la teste.
%
%  MECANISME, a ne pas re-expliquer autrement. La chaine impose
%  N = M_s/2 sur un pavage uniforme d = 2.pi/M_s, donc le terme divergent
%  de la forme fermee (T17)
%        ln N + ln|2 sin(d/2)|  ->  ln(M_s/2) + ln(2.pi/M_s) = ln(pi)
%  est INDEPENDANT de M_s : stationnaire PAR CONSTRUCTION.
%
%  /!\ NE PAS re-introduire l'explication par l'epaisseur d'anneau
%  (X_g = 0,139 contre 0,0030, facteur 47). Elle est FAUSSE et l'Article I
%  la refute deja en §3.5 : la queue s'installe des n ~ 7 sur le PMSM
%  contre n ~ 337 sur la MAS, donc PLUS TOT, pas plus tard. Si l'epaisseur
%  etait le mecanisme, cette machine deriverait DAVANTAGE.
%
%  DEVERROUILLAGE RETENU. Plutot que de porter la condensation de Schur --
%  qui exigerait de reconstruire le reseau de reluctances sur Ns noeuds au
%  lieu de Nsurf, donc de REECRIRE la chaine -- on exploite le fait que la
%  contrainte de cogging_mec porte sur numax INFERIEUR a Nsurf/2. Monter
%  numax a pavage FIXE suffit a liberer la troncature, et respecte la regle
%  du dossier : piloter la chaine existante, ne pas la refaire.
%
%  GARDE (v4 §R4). A numax = M_s/2 exactement, la chaine deverrouillee doit
%  redonner la valeur verrouillee AU BIT PRES. Sinon le portage a change
%  autre chose que la troncature.
clear; clc; t0=tic;
diary('R4_unlock_out.txt'); diary on;
M=machine_bldc(); p=M.p; Ns=M.Ns;
if isfile('kfringe_ident.mat'), S0=load('kfringe_ident.mat'); kfr=S0.kbest; else, kfr=0.325; end
thu=linspace(0,2*pi,2001); thu(end)=[];
amp=@(y,th,k)abs((2/numel(y))*sum(y(:).'.*exp(-1i*k*th(:).')));
mu0=4*pi*1e-7; L=M.ls;

fprintf('=== R4 : deverrouillage troncature-pavage ===\n\n');
fprintf('  CONFIGURATION\n');
fprintf('    machine     : PMSM 15/14, 750 W (machine_bldc)\n');
fprintf('    chaine      : cogging_mec + airgap_magnet\n');
fprintf('    base        : P0 (constante par morceaux), celle du §3.5\n');
fprintf('    solveur     : lineaire, mu_i = %g\n',M.muI);
fprintf('    k_fringe    : %.4f\n',kfr);
fprintf('    pavage      : FIXE, uniforme\n');
fprintf('    reference EF: %s\n',M.FEA.dir);

%% ======== 1. LE VERROU : ce que fait la chaine publiee ================
fprintf('\n  ---- 1. CHAINE VERROUILLEE (numax = M_s/2) ----\n');
fprintf('  %8s %8s %14s %14s %16s\n', ...
    'M_s','numax','Bg1 (T)','lam_a (Wb)','ln N + ln|2sin|');
Msl=[540 1080 2160 4320 8640];
V=nan(numel(Msl),4);
for i=1:numel(Msl)
    Msv=Msl(i);
    R=cogging_mec(M,Msv,0,3,M.muI,kfr,1e-6);      % 0 -> verrou impose
    Br=R.AG.field(R.U(1:Msv,1),0,thu); Br=Br(:).';
    d=2*pi/Msv; N=R.numax_eff;
    br=log(N)+log(abs(2*sin(d/2)));
    V(i,:)=[Msv N amp(Br,thu,p) br];
    fprintf('  %8d %8d %14.6f %14s %16.6f\n',Msv,N,V(i,3),'--',br);
end
fprintf('  crochet : %.6f -> %.6f, soit %+.4f %% sur un facteur %d en M_s\n', ...
    V(1,4),V(end,4),100*(V(end,4)-V(1,4))/V(1,4),Msl(end)/Msl(1));
fprintf('  ln(pi) = %.6f  |  ln(pi)+gamma = %.6f\n',log(pi),log(pi)+0.5772156649);
fprintf('  => STATIONNAIRE. Le terme divergent ne croit pas : c''est le verrou.\n');

%% ======== 2. LE DEVERROUILLAGE : pavage FIXE, numax libre =============
Msf=1080;                                   % pavage FIXE
nu_lock=floor(Msf/2);
fprintf('\n  ---- 2. CHAINE DEVERROUILLEE (M_s = %d FIXE, numax libre) ----\n',Msf);
fprintf('  %10s %8s %14s %16s %14s\n', ...
    'numax','n/verrou','Bg1 (T)','ln N + ln|2sin|','ecart Bg1');
nul=nu_lock*[1 2 4 8 16];
W=nan(numel(nul),4);
d=2*pi/Msf;
for i=1:numel(nul)
    N=nul(i);
    R=cogging_mec(M,Msf,N,3,M.muI,kfr,1e-6);
    Br=R.AG.field(R.U(1:Msf,1),0,thu); Br=Br(:).';
    br=log(R.numax_eff)+log(abs(2*sin(d/2)));
    W(i,:)=[R.numax_eff N/nu_lock amp(Br,thu,p) br];
    fprintf('  %10d %8.0f %14.6f %16.6f %13.4f %%\n', ...
        W(i,1),W(i,2),W(i,3),W(i,4),100*(W(i,3)-W(1,3))/W(1,3));
end
fprintf('  crochet : %.6f -> %.6f, soit %+.4f %%\n', ...
    W(1,4),W(end,4),100*(W(end,4)-W(1,4))/W(1,4));

%% ---- pente mesuree contre pente predite -----------------------------
pp=polyfit(log10(W(:,1)),W(:,4),1);
fprintf('\n  pente du crochet : %.6f par decade de N (predit : ln10 = %.6f)\n', ...
    pp(1),log(10));
fprintf('  ecart a la prediction : %+.4f %%\n',100*(pp(1)-log(10))/log(10));
qq=polyfit(log10(W(:,1)),W(:,3),1);
fprintf('  pente de Bg1     : %.6e T par decade de N\n',qq(1));
fprintf('  coefficient predit (2*mu0*L/pi)*ln10 = %.6e\n',(2*mu0*L/pi)*log(10));

%% ======== 3. GARDE ====================================================
fprintf('\n  ---- 3. GARDE : la chaine deverrouillee redonne-t-elle le verrou ? ----\n');
Rl=cogging_mec(M,Msf,0,3,M.muI,kfr,1e-6);          % verrou impose
Ru=cogging_mec(M,Msf,nu_lock,3,M.muI,kfr,1e-6);    % deverrouille A LA MEME valeur
Bl=Rl.AG.field(Rl.U(1:Msf,1),0,thu); Bl=Bl(:).';
Bu=Ru.AG.field(Ru.U(1:Msf,1),0,thu); Bu=Bu(:).';
fprintf('    numax verrouille   : %d\n',Rl.numax_eff);
fprintf('    numax deverrouille : %d\n',Ru.numax_eff);
fprintf('    Bg1 verrouille     : %.15f T\n',amp(Bl,thu,p));
fprintf('    Bg1 deverrouille   : %.15f T\n',amp(Bu,thu,p));
fprintf('    ecart              : %.3e T\n',abs(amp(Bl,thu,p)-amp(Bu,thu,p)));
fprintf('    ecart max sur Br   : %.3e T\n',max(abs(Bl-Bu)));
if max(abs(Bl-Bu))<1e-12
    fprintf('    GARDE PASSEE : identite au bit pres.\n');
else
    fprintf('    GARDE ECHOUEE : le portage a change autre chose.\n');
end

%% ======== 4. PANNEAU POUR LA TABLE 8 =================================
fprintf('\n  ---- 4. PANNEAU A AJOUTER A LA TABLE 8 ----\n\n');
fprintf('\\multicolumn{4}{l}{\\emph{Truncation unlocked from the tiling}, $M_s = %d$}\\\\\n',Msf);
fprintf('$N_h$ & $\\ln N + \\ln|2\\sin(d/2)|$ & $B_{g1}$ (T) & dev. \\\\\n');
for i=1:numel(nul)
    fprintf('%d & %.4f & %.5f & $%+.2f\\%%$ \\\\\n', ...
        W(i,1),W(i,4),W(i,3),100*(W(i,3)-W(1,3))/W(1,3));
end
save('R4_unlock.mat','V','W','Msl','Msf','nul');
fprintf('\n  duree %.0f s\n=== R4 termine ===\n',toc(t0));
diary off;
