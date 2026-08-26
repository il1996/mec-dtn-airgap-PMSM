%% RUN_V1_PMSM_BASE - V1 : la base P1 sur le bore du PMSM
%
%  VERROU. T16/T17 n'ont ete menes que sur la MAS. Hypothese a tester : les
%  bons resultats PMSM tiennent en partie a une troncature grossiere qui
%  MASQUE la queue logarithmique.
%
%  ECART A LA SPECIFICATION, ET SA RAISON. V1 demande de faire varier N SEUL
%  a pavage fixe. C'est IMPOSSIBLE sur cette chaine : cogging_mec impose
%  numax = Nsurf/2 (ligne 27), car en dessous l'operateur de couronne est de
%  RANG DEFICIENT et au-dessus les harmoniques replient sur une base a Nsurf
%  colonnes. La troncature est donc VERROUILLEE sur la discretisation de
%  surface -- difference structurelle avec la MAS, ou airgap_dtn_tooth
%  condense une grille fine et laisse N libre.
%
%  On fait donc varier Nsurf, avec N = Nsurf/2 qui suit, et on le DECLARE :
%  le balayage mesure l'effet conjoint pavage+troncature, non la troncature
%  seule. C'est une limite de la chaine PMSM, et c'est en soi un resultat.
clear; clc; t0=tic;
diary('V1_pmsm_basis_out.txt'); diary on;
M=machine_bldc(); p=M.p;
if isfile('kfringe_ident.mat'), S0=load('kfringe_ident.mat'); kfr=S0.kbest; else, kfr=0.325; end
thu=linspace(0,2*pi,2001); thu(end)=[];
amp=@(y,th,k)abs((2/numel(y))*sum(y(:).'.*exp(-1i*k*th(:).')));

fprintf('=== V1 : base P0 vs P1 sur le bore du PMSM ===\n');
fprintf('  machine 15/14 750 W | kfringe = %.4f | solveur lineaire\n',kfr);
fprintf('  reference EF : Bg1 1.0746 | rang8 0.01966 T\n');
fprintf('  NOTE : numax = Nsurf/2 impose par cogging_mec (rang de Y).\n');
fprintf('         Le balayage varie donc pavage ET troncature ensemble.\n\n');

Nsurf=[540 1080 1260 2160 4320 8640];
lab={'Bg1 (T)','|Br| moyen (T)','Br crete (T)','Bt rms (T)','rang 8 (T)'};
V=nan(numel(Nsurf),5,2);

for b=1:2
    bas={'p0','p1'}; bs=bas{b};
    fprintf('  --- base %s ---\n',upper(bs));
    fprintf('  %7s %6s %10s %12s %11s %10s %11s\n', ...
        'Nsurf','N','Bg1','|Br| moyen','Br crete','Bt rms','rang 8');
    for k=1:numel(Nsurf)
        Nk=Nsurf(k);
        R=cogging_mec(M,Nk,0,3,M.muI,kfr,1e-6,bs);
        [Br,Bt]=R.AG.field(R.U(1:Nk,1),0,thu);
        Br=Br(:).'; Bt=Bt(:).';
        V(k,:,b)=[amp(Br,thu,p) mean(abs(Br)) max(Br) sqrt(mean(Bt.^2)) amp(Br,thu,8)];
        fprintf('  %7d %6d %10.5f %12.5f %11.5f %10.5f %11.6f\n', ...
            Nk,R.numax_eff,V(k,:,b));
    end
    fprintf('\n');
end

%% ---- diagnostic : increments et leurs rapports -----------------------
fprintf('  --- increments de Bg1 et leurs rapports ---\n');
fprintf('  %8s %13s %9s %13s %9s\n','Nsurf','dP0','rap.P0','dP1','rap.P1');
d0=diff(V(:,1,1)); d1=diff(V(:,1,2));
for k=1:numel(d0)
    r0=NaN; r1=NaN;
    if k<numel(d0), r0=d0(k+1)/d0(k); r1=d1(k+1)/d1(k); end
    fprintf('  %8d %13.3e %9.3f %13.3e %9.3f\n',Nsurf(k+1),d0(k),r0,d1(k),r1);
end
dsp=@(x)100*(max(x)-min(x))/mean(x);
fprintf('\n  dispersion sur les %d discretisations :\n',numel(Nsurf));
for j=1:5
    fprintf('  %-18s P0 %8.3f %%   P1 %8.3f %%\n',lab{j},dsp(V(:,j,1)),dsp(V(:,j,2)));
end
fprintf('\n  valeurs a la plus fine discretisation (Nsurf = %d) :\n',Nsurf(end));
for j=1:5
    fprintf('  %-18s P0 %10.5f   P1 %10.5f   ecart %+7.2f %%\n', ...
        lab{j},V(end,j,1),V(end,j,2),100*(V(end,j,2)-V(end,j,1))/V(end,j,1));
end
fprintf(['\n  LECTURE. Rapports -> 1 : queue logarithmique. Rapports ~0,3 :\n' ...
         '  convergence geometrique. Si P0 disperse peu, c''est que la chaine\n' ...
         '  PMSM opere SOUS le seuil ou le defaut se manifeste.\n']);
save('V1_pmsm.mat','V','Nsurf');
fprintf('\n  duree %.0f s\n=== V1 termine ===\n',toc(t0));
diary off;
