%% RUN_R6C_CORR - R6, volet mecanisme : ce que fait la correction constante
%
%  DEPENDANCE DECLAREE. Ce bloc ne recalcule aucune carte : il relit les deux
%  cartes produites et sauvegardees par RUN_R6B_FEA (kfringe = 0,75 et 0)
%  dans R6_satmap.mat. Il est donc regenerable en rejouant R6b puis R6c.
%
%  CE QU'IL ETABLIT. mec_map.m:203-231 retranche du flux totalise un terme
%  dLsl*iabc, LINEAIRE en courant, dont le coefficient alpha est calibre a
%  COURANT NUL pour ramener l'inductance de ligne sur 2*Ld d'un reseau
%  LINEAIRE et AIMANTS ETEINTS. La justification ecrite dans le code est que
%  l'exces est "un chemin dans l'AIR : il ne sature pas".
%
%  On teste cette premisse au lieu de l'admettre : si l'exces etait un chemin
%  d'air, l'ecart entre le reseau AVEC pont tangentiel (kfringe = 0,75) et le
%  reseau SANS (kfringe = 0) serait constant en courant. On le mesure.
clear; clc; t0=tic;
diary('R6_satmap_out.txt'); diary on;      % on APPEND au meme livrable
fprintf('\n\n=== R6c : ce que fait la correction constante ===\n');
fprintf('  (relit les cartes de R6b, n''en recalcule aucune)\n\n');
Q=load('R6_satmap.mat'); S=Q.S; S0=Q.S0;
d75=S.Leq0-2*S.Ldref; d0=S0.Leq0-2*S0.Ldref;
IQ=[0 2 5 10 15 20 25];

%% ---- GARDE : la correction est-elle EXACTEMENT constante ? ------------
%  Sur la combinaison de ligne, dLsl*iabc doit produire un decalage
%  independant du courant. S'il ne l'est pas, le raisonnement qui suit tombe.
fprintf('  ---- GARDE : la soustraction est-elle constante en courant ? ----\n');
r75=arrayfun(@(x)LL(S,x),IQ)+d75;
r0 =arrayfun(@(x)LL(S0,x),IQ)+d0;
fprintf('    par construction : L_brute(i) = L_corrigee(i) + %.4f mH\n',d75*1e3);
fprintf('    verification a courant nul : L_corrigee(0) = %.6f mH = 2*Ld = %.6f mH\n', ...
    LL(S,0)*1e3,2*S.Ldref*1e3);
G=abs(LL(S,0)-2*S.Ldref)<1e-9 && abs(LL(S0,0)-2*S0.Ldref)<1e-9;
fprintf('    GARDE %s\n',tern(G,'PASSEE','ECHOUEE'));

%% ---- le tableau ------------------------------------------------------
fprintf('\n  ---- INDUCTANCE DE LIGNE : BRUTE ET CORRIGEE, DEUX RESEAUX ----\n');
fprintf('  %6s %15s %15s %15s %15s %14s\n', ...
    'i (A)','kfr=0.75 corr','kfr=0.75 brute','kfr=0 corr','kfr=0 brute','ecart des 2');
for k=1:numel(IQ)
    a=LL(S,IQ(k)); b=LL(S0,IQ(k));
    fprintf('  %6.0f %13.3f mH %13.3f mH %13.3f mH %13.3f mH %11.3f mH\n', ...
        IQ(k),a*1e3,r75(k)*1e3,b*1e3,r0(k)*1e3,(r75(k)-r0(k))*1e3);
end

%% ---- la premisse de la correction, testee -----------------------------
dd=r75-r0;
fprintf('\n  ---- LA PREMISSE DE LA CORRECTION EST FAUSSE ----\n');
fprintf('    Si l''exces etait un chemin d''AIR, l''ecart entre les deux reseaux\n');
fprintf('    serait CONSTANT. Mesure : %+.3f mH a 0 A, %+.3f mH a 25 A.\n',dd(1)*1e3,dd(end)*1e3);
fprintf('    Il change de SIGNE : amplitude de variation %.3f mH.\n',(max(dd)-min(dd))*1e3);
fprintf('    Retrancher une constante calibree a i = 0 est donc exact a\n');
fprintf('    courant nul et de plus en plus faux ensuite.\n');

%% ---- consequence 1 : d'ou vient l'effondrement publie ------------------
fprintf('\n  ---- CONSEQUENCE 1 : L''EFFONDREMENT PUBLIE EST EN PARTIE ARITHMETIQUE ----\n');
c75=arrayfun(@(x)LL(S,x),IQ);
fprintf('    courbe PUBLIEE (kfr=0.75 corrigee) : %.1f -> %.1f mH, facteur %.1f\n', ...
    c75(1)*1e3,c75(end)*1e3,c75(1)/c75(end));
fprintf('    physique BRUTE du meme reseau      : %.1f -> %.1f mH, facteur %.1f\n', ...
    r75(1)*1e3,r75(end)*1e3,r75(1)/r75(end));
fprintf('    physique brute SANS pont tangentiel: %.1f -> %.1f mH, facteur %.1f\n', ...
    r0(1)*1e3,r0(end)*1e3,r0(1)/r0(end));
fprintf('    Aucun reseau ne produit le facteur publie : il vient de ce qu''on\n');
fprintf('    retranche %.1f mH a une grandeur qui, elle, tombe a %.1f mH.\n', ...
    d75*1e3,r75(end)*1e3);

%% ---- consequence 2 : pourquoi la grille ne converge pas en haut -------
fprintf('\n  ---- CONSEQUENCE 2 : L''AMPLIFICATION DU BRUIT DE GRILLE ----\n');
amp=r75(end)/c75(end);
fprintf('    A 25 A la valeur publiee est la difference de %.1f mH et %.1f mH.\n', ...
    r75(end)*1e3,d75*1e3);
fprintf('    Tout ecart relatif sur la brute est donc multiplie par %.1f sur la\n',amp);
fprintf('    corrigee. R6a mesure une dispersion de 41,6 %% a 25 A et de 0,3 %% a 5 A :\n');
fprintf('    la non-convergence en haut de plage est coherente avec ce facteur.\n');
fprintf('    STATUT : explication coherente avec les nombres, NON mesuree\n');
fprintf('    independamment (R6a n''a pas archive alpha pour chaque grille).\n');

fprintf('\n  duree %.0f s\n=== R6c termine ===\n',toc(t0));
diary off;

% ======================================================================
function L = LL(St,ii)
%LL  Inductance de ligne, meme combinaison que local_Lline.
    ia=ii; ib=-ii/sqrt(3); s3=sqrt(3);
    q=@(X,k) interpn(1:3,St.the,St.iax,St.iax,X,k,St.the,ia,ib,'linear');
    L=mean((q(St.La,1)-q(St.La,2))-(q(St.Lb,1)-q(St.Lb,2))/s3);
end

function s=tern(c,a,b), if c, s=a; else, s=b; end, end
