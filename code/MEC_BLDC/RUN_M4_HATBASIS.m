%% RUN_M4_HATBASIS - M-4 : troncature liberee EN BASE CHAPEAU sur le PMSM
%
%  CE QUE CE BLOC DECIDE. R4 a libere la troncature du pavage en base P0 et
%  montre que le crochet reprend sa croissance a la pente ln 10. La question
%  que le manuscrit laisse ouverte est l'autre moitie : EN BASE CHAPEAU, la
%  serie converge-t-elle sur CETTE machine ? C'est la seule experience qui
%  permettrait a l'Article I de demontrer son resultat central sur sa propre
%  machine, au lieu de l'importer de l'article compagnon.
%
%  LES DEUX ISSUES SONT PUBLIABLES, ET IL FAUT LE DIRE AVANT DE LANCER.
%    - Si B_g1 se stabilise quand N croit a pavage fixe, la base chapeau
%      tient sa promesse ici aussi et le §10 y gagne son premier item.
%    - Si elle ne se stabilise PAS, c'est un resultat plus important que
%      tout le reste du papier et il se publie tel quel.
%
%  GARDE. Au point de recouvrement N = M_s/2, la chaine liberee doit
%  redonner la chaine verrouillee AU BIT PRES, dans chaque base. C'est la
%  garde de R4, reprise ici : elle prouve que seule la troncature a change.
clear; clc; t0=tic;
diary('M4_hatbasis_out.txt'); diary on;
M=machine_bldc(); p=M.p;
if isfile('kfringe_ident.mat'), S0=load('kfringe_ident.mat'); kfr=S0.kbest; else, kfr=0.325; end
thu=linspace(0,2*pi,2001); thu(end)=[];
amp=@(y,th,k)abs((2/numel(y))*sum(y(:).'.*exp(-1i*k*th(:).')));
Ms=1080; nu_lock=floor(Ms/2); MUL=[1 2 4 8];

fprintf('=== M-4 : troncature liberee, base P0 CONTRE base chapeau ===\n\n');
fprintf('  CONFIGURATION\n');
fprintf('    machine   : PMSM 15/14, 750 W (machine_bldc)\n');
fprintf('    chaine    : cogging_mec + airgap_magnet\n');
fprintf('    pavage    : M_s = %d, FIXE\n',Ms);
fprintf('    verrou    : N = M_s/2 = %d\n',nu_lock);
fprintf('    N balaye  : verrou x %s\n',mat2str(MUL));
fprintf('    solveur   : lineaire, mu_r = %g | k_fringe = %.4f\n',M.muI,kfr);

BAS={'p0','p1'}; NM={'P0 (constante)','P1 (chapeau)'};
R=nan(numel(MUL),2);
for b=1:2
    fprintf('\n  ---- base %s ----\n',NM{b});
    fprintf('  %8s %10s %16s %14s\n','N','N/verrou','Bg1 (T)','ecart au verrou');
    for k=1:numel(MUL)
        N=nu_lock*MUL(k);
        S=cogging_mec(M,Ms,N,3,M.muI,kfr,1e-6,BAS{b});
        Br=S.AG.field(S.U(1:Ms,1),0,thu); Br=Br(:).';
        R(k,b)=amp(Br,thu,p);
        fprintf('  %8d %10d %16.6f %13.4f %%\n', ...
            N,MUL(k),R(k,b),100*(R(k,b)-R(1,b))/R(1,b));
    end
end

%% ---- diagnostic : increments et leurs rapports ------------------------
fprintf('\n  ---- INCREMENTS DE Bg1 ET LEURS RAPPORTS ----\n');
fprintf('  Rapport -> 1 : queue logarithmique, pas de limite.\n');
fprintf('  Rapport < 1 et decroissant : convergence.\n\n');
fprintf('  %10s %14s %10s %14s %10s\n','N','dP0','rap. P0','dP1','rap. P1');
d0=diff(R(:,1)); d1=diff(R(:,2));
for k=1:numel(d0)
    r0=NaN; r1=NaN;
    if k<numel(d0), r0=d0(k+1)/d0(k); r1=d1(k+1)/d1(k); end
    fprintf('  %10d %14.3e %10.3f %14.3e %10.3f\n',nu_lock*MUL(k+1),d0(k),r0,d1(k),r1);
end
dsp=@(x)100*(max(x)-min(x))/mean(x);
fprintf('\n  dispersion de Bg1 sur le balayage :\n');
fprintf('    base P0      %8.4f %%\n',dsp(R(:,1)));
fprintf('    base chapeau %8.4f %%\n',dsp(R(:,2)));

%% ---- GARDE -------------------------------------------------------------
fprintf('\n  ---- GARDE : recouvrement au verrou, dans chaque base ----\n');
G=true;
for b=1:2
    Sl=cogging_mec(M,Ms,0,3,M.muI,kfr,1e-6,BAS{b});        % verrou impose
    Su=cogging_mec(M,Ms,nu_lock,3,M.muI,kfr,1e-6,BAS{b});  % libere au verrou
    Bl=Sl.AG.field(Sl.U(1:Ms,1),0,thu); Bl=Bl(:).';
    Bu=Su.AG.field(Su.U(1:Ms,1),0,thu); Bu=Bu(:).';
    e=max(abs(Bl-Bu));
    fprintf('    base %-14s ecart max sur Br : %.3e T  %s\n', ...
        NM{b},e,ternstr(e<1e-12,'identite au bit pres','ECHEC'));
    if e>=1e-12, G=false; end
end
fprintf('    GARDE %s\n',ternstr(G,'PASSEE','ECHOUEE'));

%% ---- verdict -----------------------------------------------------------
fprintf('\n  ---- VERDICT ----\n');
rp1=d1(end)/d1(end-1);
if abs(rp1)<0.7
    fprintf('    La base chapeau CONVERGE sur cette machine : rapport\n');
    fprintf('    d''increments %.3f, dispersion %.4f %% contre %.4f %% en P0.\n', ...
        rp1,dsp(R(:,2)),dsp(R(:,1)));
    fprintf('    L''Article I peut demontrer son resultat central sur sa\n');
    fprintf('    PROPRE machine, et le premier item de << further work >>\n');
    fprintf('    est acquitte.\n');
else
    fprintf('    La base chapeau NE CONVERGE PAS de facon nette ici :\n');
    fprintf('    rapport d''increments %.3f. C''est un resultat plus\n',rp1);
    fprintf('    important que le reste et il doit etre publie tel quel,\n');
    fprintf('    non enfoui. Voir la note de ce bloc.\n');
end

%% ---- panneau LaTeX -----------------------------------------------------
fprintf('\n  ---- PANNEAU POUR LE PANNEAU BAS DE LA TABLE 8 ----\n\n');
fprintf('$N_h$ & $B_{g1}$ p.c.\\ (\\si{\\tesla}) & $B_{g1}$ hat (\\si{\\tesla}) & dev.\\ hat\\\\\n');
for k=1:numel(MUL)
    fprintf('%d & \\num{%.5f} & \\num{%.5f} & $%+.3f\\%%$\\\\\n', ...
        nu_lock*MUL(k),R(k,1),R(k,2),100*(R(k,2)-R(1,2))/R(1,2));
end

save('M4_hatbasis.mat','R','MUL','Ms','nu_lock','G');
fprintf('\n  duree %.0f s\n=== M-4 termine ===\n',toc(t0));
diary off;

% ======================================================================
function s=ternstr(c,a,b), if c, s=a; else, s=b; end, end
