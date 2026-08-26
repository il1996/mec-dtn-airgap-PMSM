%% RUN_IND_MESH - inductances sur le maillage raffine
%
%  *** ETAT : VALIDE. Ce script ALIMENTE la colonne Mesh du Tableau 7. ***
%
%  HISTORIQUE, conserve parce qu'il explique un en-tete qui a longtemps dit
%  le contraire. Une version anterieure donnait La = 12,25 mH contre 50,21
%  mesures par la FEA (-75,6 %) : la FMM etait injectee dans les branches
%  TANGENTIELLES avec la meme valeur a TOUTES les couches, alors que la
%  fraction d'ampere-tours enlacee doit croitre avec la profondeur dans
%  l'encoche -- c'est ce qui donne le h/(3b) classique. Le motif Fcol a
%  depuis ete gradue en profondeur et le defaut est CORRIGE.
%
%  ETAT ACTUEL, mesure : La = 50,38 a 51,55 mH selon Ms et n_sh, soit +0,3 %
%  a +2,7 % contre la reference. La ligne Ms=540 n_sh=1 (Ld = 52,395 mH,
%  +0,1 %) est celle publiee en Table 7.
%
%  L'en-tete "NE PAS UTILISER" a survecu a la correction pendant plusieurs
%  mois sur un script qui produisait des donnees publiees. C'est le pire des
%  deux mondes : soit un relecteur ecarte des resultats valides, soit il
%  conclut que le manuscrit s'appuie sur du code que ses auteurs declarent
%  invalide. Tache C2 de SPEC_CLAUDE_CODE_v2.
%
%  INTERET. Dans le reseau a une dent, l'inductance est la somme de deux
%  morceaux calcules separement : une part d'entrefer issue du reseau, et
%  une fuite d'encoche AJOUTEE analytiquement (permeance de trapeze +
%  lambda_tip), avec kfringe force a 0 pour ne pas compter deux fois la
%  fuite de bec. Trois conventions a tenir, documentees mais fragiles.
%  Sur le maillage, la fuite d'encoche est RESOLUE : les colonnes d'air de
%  l'encoche et leurs branches tangentielles portent explicitement le flux
%  qui traverse l'encoche d'un flanc a l'autre. Plus de formule ajoutee,
%  plus de kfringe a annuler.
%
%  FLUX TOTALISE. La FMM de bobinage est injectee dans les branches
%  TANGENTIELLES (Fcol distribue sur les colonnes d'encoche). Pour un tel
%  reseau, le flux totalise d'une phase s'obtient EXACTEMENT par
%      lambda_ph = sum_k F_k^(ph) * Phi_k / i
%  ou F_k^(ph) est le motif de FMM de la phase a 1 A. C'est la forme duale
%  de W = 1/2 sum(F.Phi) : elle englobe toutes les fuites, y compris celles
%  qui ne traversent aucune dent.
clear; clc; t0=tic;
M=machine_bldc(); Ns=M.Ns;
if isfile('kfringe_ident.mat'), S0=load('kfringe_ident.mat'); kfr=S0.kbest; else, kfr=0.325; end
rd=@(f)readmatrix(f,'FileType','text','Delimiter','\t','NumHeaderLines',1);
dL=rd(fullfile(M.FEA.dir,'magnetostique(Armature-Field)','Output Variables Table 1.tab'));
LaF=dL(1,2); MF=(dL(1,3)+dL(1,4))/2; LdF=LaF-MF;
fprintf('=== Inductances : maillage raffine vs reseau a une dent vs FEA ===\n');
fprintf('  FEA : La = %.3f mH | M = %.4f mH | Ld = %.3f mH\n\n',LaF*1e3,MF*1e3,LdF*1e3);

%  motifs de FMM de branche des trois phases, a 1 A, sur le MEME maillage
mkF=@(Ms,nsh,ii)mesh_bldc(M,Ms,4,3,ii,0,kfr,M.muI,nsh);
fprintf('  %5s %5s %6s | %9s %10s %9s | %8s %8s %8s\n', ...
    'Ms','nsh','noeuds','La (mH)','M (mH)','Ld (mH)','ec. La','ec. M','ec. Ld');
for Ms=[360 540 900]
    for nsh=[1 2]
        MEa=mkF(Ms,nsh,[1;0;0]); MEa.Isrc(:)=0;   % aimants eteints
        Sa=solve_bldc_mesh(MEa);
        %  motifs de FMM des phases b et c sur la meme topologie
        Fb=mkF(Ms,nsh,[0;1;0]).E;  Fc=mkF(Ms,nsh,[0;0;1]).E;
        la=sum(MEa.E.*Sa.Phi);                     % lambda_a a 1 A
        lb=sum(Fb.*Sa.Phi); lc=sum(Fc.*Sa.Phi);
        La=la; Mm=0.5*(lb+lc); Ld=La-Mm;
        fprintf('  %5d %5d %6d | %9.3f %10.4f %9.3f | %7.1f %% %7.1f %% %7.1f %%\n', ...
            Ms,nsh,MEa.N,La*1e3,Mm*1e3,Ld*1e3, ...
            100*(La-LaF)/LaF,100*(Mm-MF)/abs(MF),100*(Ld-LdF)/LdF);
    end
end

%% ---- reference interne : reseau a une dent ----------------------------
RI=inductance_mec(M,1260,M.muI,0);
fprintf('\n  %-22s %9.3f %10.4f %9.3f | %7.1f %% %7.1f %% %7.1f %%\n', ...
    'reseau a une dent',RI.La*1e3,RI.M*1e3,RI.Ld*1e3, ...
    100*(RI.La-LaF)/LaF,100*(RI.M-MF)/abs(MF),100*(RI.Ld-LdF)/LdF);
fprintf('     (dont fuite d''encoche ANALYTIQUE : %.2f mH sur %.2f, soit %.0f %%)\n', ...
    RI.L_slot*1e3,RI.La*1e3,100*RI.L_slot/RI.La);

%% ---- decomposition sur le maillage ------------------------------------
MEa=mkF(540,2,[1;0;0]); MEa.Isrc(:)=0; Sa=solve_bldc_mesh(MEa);
La=sum(MEa.E.*Sa.Phi);
%  C2. Les deux lignes precedentes se terminaient par "*0" : le diagnostic
%  valait ZERO PAR CONSTRUCTION, alors que la conclusion imprimee juste en
%  dessous porte l'argument central de la §3.5 (fuite d'encoche RESOLUE et
%  non ajoutee). Une affirmation etayee par un calcul desactive est pire
%  qu'une affirmation nue. Le calcul est retabli ci-dessous.
%  QUELLE GRANDEUR DECOMPOSER ? Pas lambda. lambda = somme(F.Phi) est une
%  forme DUALE, ponderee par l'endroit ou la FMM est INJECTEE et non par
%  celui ou le flux PASSE. La FMM de bobinage etant injectee sur des
%  branches de fer, E est nul sur toute branche d'air et la part d'air
%  vaudrait 0 PAR CONSTRUCTION -- un calcul juste qui ne mesure pas ce que
%  la ligne suivante annonce. La grandeur qui repond est l'ENERGIE
%  MAGNETIQUE, W = 1/2 * Phi * dU, sommee par type de branche.
mAir=~MEa.iron;
%  ATTENTION. La chute aux bornes de la RELUCTANCE d'une branche a source
%  vaut U(a) - U(b) + E, et non U(a) - U(b) : omettre le terme source donne
%  2W = -19,45 mH contre lambda = 50,48. Le controle 2W = lambda ci-dessous
%  est precisement la garde qui attrape cette erreur -- il est conserve.
dU=Sa.U(MEa.a)-Sa.U(MEa.b)+MEa.E(:);
W=0.5*Sa.Phi(:).*dU(:);
%  L'ENTREFER N'EST PAS UNE BRANCHE. L'operateur de couronne AG.Y est un
%  bloc DENSE reliant tous les noeuds de surface : son energie 1/2 U'YU
%  n'apparait dans aucune somme sur les branches. L'omettre laisse
%  2W = 31,03 mH contre lambda = 50,48 -- il manque exactement 19,45 mH,
%  qui est l'energie de la couronne. Le controle 2W = lambda ne ferme
%  qu'une fois ce terme ajoute, et c'est lui qui le prouve.
%  SIGNE. solve_bldc_mesh assemble le bloc de couronne par
%  A(ids,ids) = A(ids,ids) - Y (l.84), c'est-a-dire avec le signe OPPOSE
%  a celui des conductances de branche. L'energie de la couronne dans la
%  convention du reseau est donc -1/2 U'YU. Avec le signe brut on obtient
%  une energie NEGATIVE pour un operateur semi-defini POSITIF -- absurdite
%  qui designe la convention, et que le controle 2W = lambda confirme.
Wann=-0.5*(Sa.Usurf(:).'*(MEa.AG.Y*Sa.Usurf(:)));
Wbr_air=sum(W(mAir)); Wfer=sum(W(~mAir));
Wair=Wbr_air+Wann; Wtot=Wair+Wfer;
fprintf('\n  --- ou passe l''ENERGIE d''induit (maillage, Ms=540, nsh=2) ---\n');
fprintf('  energie TOTALE                          : %12.6e J/A^2\n',Wtot);
fprintf('  AIR : branches (encoche + becs)         : %12.6e  (%.1f %%)\n', ...
    Wbr_air,100*Wbr_air/Wtot);
fprintf('  AIR : couronne d''entrefer (1/2 U''YU)    : %12.6e  (%.1f %%)\n', ...
    Wann,100*Wann/Wtot);
fprintf('  FER : dents + culasses                  : %12.6e  (%.1f %%)\n', ...
    Wfer,100*Wfer/Wtot);
fprintf('  ---------------------------------------------------------------\n');
fprintf('  AIR TOTAL                               : %12.6e  (%.1f %%)\n', ...
    Wair,100*Wair/Wtot);
fprintf('  CONTROLE energetique 2W vs lambda       : %12.4f vs %.4f mH  (%+.2e)\n', ...
    2*Wtot*1e3,La*1e3,2*Wtot-La);
fprintf('\n  (pour memoire, la decomposition de lambda = somme(F.Phi) donne\n');
fprintf('   0 %% d''air : c''est une propriete de la forme duale, PAS un\n');
fprintf('   resultat physique. Elle a ete retiree.)\n');
%  Comparaison au modele localise, qui AJOUTE une permeance d'encoche
%  analytique : c'est la grandeur homologue, et le rapprochement est le
%  seul qui etaye reellement la §3.5.
fprintf('  pour memoire, reseau a une dent : fuite d''encoche ANALYTIQUE\n');
fprintf('    %.2f mH sur %.2f, soit %.0f %% -- une FORMULE AJOUTEE.\n', ...
    RI.L_slot*1e3,RI.La*1e3,100*RI.L_slot/RI.La);
fprintf('  Sur le maillage la meme fuite est portee par des CELLULES D''AIR\n');
fprintf('  resolues, sans permeance ajoutee ni kfringe a annuler.\n');
%  NOTE. La part "traversant la surface d'entrefer" a ete RETIREE et non
%  reparee : les branches de surface portent E = 0, donc le produit F.Phi
%  y est nul par definition et ne mesure rien. La quantifier demanderait le
%  flux d'entrefer par l'operateur (AG.Y*Usurf), grandeur d'une autre
%  nature qui ne se compare pas a une part de lambda.
fprintf('\n  (%.0f s)\n',toc(t0));
