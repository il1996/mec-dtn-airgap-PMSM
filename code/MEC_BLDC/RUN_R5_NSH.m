%% RUN_R5_NSH - trancher n_sh : deux couches ou quatre (bloc R5, v4)
%
%  DEFAUT SIGNALE. Le texte de l'Article I annonce quelque part n_sh = 4
%  tandis que la Table 12 rapporte n_sh = 2. La v4 fournit la garde qui
%  tranche : "le nombre de couches radiales declare dans le texte doit etre
%  coherent avec le nombre d'inconnues declare (5400 a n_sh=1, 6480 a la
%  seconde configuration). Ce dernier chiffre tranche a lui seul."
%
%  METHODE. On ne lit pas le manuscrit pour trancher, on le CONFRONTE :
%  chaque chiffre de comparaison est RELU du .tex par expression reguliere
%  et chaque chiffre de reference est RECALCULE par la chaine. Aucun nombre
%  n'est transcrit a la main, dans un sens ni dans l'autre.
%
%  GARDE. Le compte d'inconnues produit par mesh_bldc doit reproduire les
%  comptes declares dans le manuscrit, ET le B_g1 des colonnes publiees doit
%  sortir de la configuration que le compte designe. Si les deux gardes
%  passent, le chiffre d'inconnues designe la configuration sans ambiguite
%  et l'occurrence divergente est le defaut.
clear; clc; t0=tic;
diary('R5_nsh_out.txt'); diary on;
M=machine_bldc(); Ms=540;
if isfile('kfringe_ident.mat'), S0=load('kfringe_ident.mat'); kfr=S0.kbest; else, kfr=0.325; end
TEX=fullfile('..','article','ArticleI_DtN_PMSM.tex');
thu=linspace(0,2*pi,3601); thu(end)=[];
amp=@(B,th,k) 2*abs(mean(B(:).'.*exp(-1i*k*th(:).')));

fprintf('=== R5 : trancher n_sh (2 ou 4) ===\n\n');
fprintf('  CONFIGURATION\n');
fprintf('    machine     : PMSM 15/14, 750 W (machine_bldc)\n');
fprintf('    chaine      : mesh_bldc + solve_bldc_mesh\n');
fprintf('    pavage      : Ms = %d, nst = 4, nys = 3, phi = 0\n',Ms);
fprintf('    solveur     : lineaire, mu_r = %g\n',M.muI);
fprintf('    k_fringe    : %.4f\n',kfr);
fprintf('    base        : P0\n');
fprintf('    manuscrit   : %s\n',TEX);
if ~isfile(TEX), error('R5:tex','manuscrit introuvable : %s',TEX); end
tx=fileread(TEX); tl=regexp(tx,'\r\n|\n|\r','split');

%% ======== 1. CE QUE LA CHAINE PRODUIT ================================
fprintf('\n  ---- 1. INCONNUES PRODUITES PAR mesh_bldc ----\n');
fprintf('  Ls = 2*n_sh + nst + nys ; rangees de noeuds = Ls + 1\n\n');
fprintf('  %6s %6s %10s %12s %14s\n','n_sh','Ls','rangees','inconnues','Bg1 (T)');
NSH=[1 2 3 4]; K=numel(NSH); Nu=nan(1,K); Lsv=nan(1,K); Bg=nan(1,K);
for k=1:K
    ME=mesh_bldc(M,Ms,4,3,[],0,kfr,M.muI,NSH(k));
    S=solve_bldc_mesh(ME,1e-9,30);
    Br=ME.AG.field(S.Usurf,0,thu);
    Nu(k)=ME.N; Lsv(k)=ME.Ls; Bg(k)=amp(Br(:).',thu,M.p);
    fprintf('  %6d %6d %10d %12d %14.5f\n',NSH(k),ME.Ls,ME.N/Ms,ME.N,Bg(k));
end

%% ======== 2. CE QUE LE MANUSCRIT DECLARE =============================
fprintf('\n  ---- 2. COMPTES D''INCONNUES RELUS DU MANUSCRIT ----\n');
%  (a) la note de la Table 12 : "Unknowns: 5400 at n_sh=1, 6480 at n_sh=2"
t=regexp(tx,['Unknowns:\s*\\num\{(\d+)\}\s*at\s*\$n_\{\\mathrm\{sh\}\}=(\d+)\$' ...
             '\s*,\s*\\num\{(\d+)\}\s*at\s*\$n_\{\\mathrm\{sh\}\}=(\d+)\$'],'tokens','once');
if isempty(t), error('R5:parse','note de la Table 12 non reconnue.'); end
dec=[str2double(t{2}) str2double(t{1}); str2double(t{4}) str2double(t{3})];
fprintf('  note de la Table 12 (tab:pmsm) :\n');
for i=1:2
    j=find(NSH==dec(i,1));
    ok=~isempty(j) && Nu(j)==dec(i,2);
    fprintf('    n_sh = %d -> %d inconnues declarees | chaine : %d  [%s]\n', ...
        dec(i,1),dec(i,2),Nu(j),ternstr(ok));
end
%  (b) les lignes des deux tables de convergence
rw=regexp(tx,'(?<Ms>\d+)\s*&\s*\\num\{(?<U>\d+)\}','names');
fprintf('\n  lignes "Ms & inconnues" trouvees dans le manuscrit :\n');
fprintf('  %8s %12s %10s %8s\n','Ms','inconnues','rangees','n_sh');
nshRow=nan(1,numel(rw));
for i=1:numel(rw)
    Msi=str2double(rw(i).Ms); Ui=str2double(rw(i).U);
    rows=Ui/Msi; ns=(rows-1-4-3)/2;                 % Ls+1 = 2*n_sh+7+1
    nshRow(i)=ns;
    fprintf('  %8d %12d %10g %8g\n',Msi,Ui,rows,ns);
end
fprintf('  => n_sh distincts impliques par les lignes : %s\n',mat2str(unique(nshRow)));

%% ======== 3. LA GARDE : le compte designe-t-il la bonne colonne ? =====
fprintf('\n  ---- 3. GARDE : B_g1 des colonnes publiees ----\n');
b=regexp(tx,'\$B_\{g1\}\$\s*\(T\)\s*&\s*([\d.]+)\s*&\s*([\d.]+)\s*&\s*([\d.]+)\s*&\s*([\d.]+)', ...
    'tokens','once');
if isempty(b), error('R5:parse','ligne B_g1 de la Table 12 non reconnue.'); end
pub=cellfun(@str2double,b(1:2));
fprintf('  B_g1 publie, colonnes "Mesh" de la Table 12 : %.5f | %.5f T\n',pub);
fprintf('  %6s %12s %14s %12s\n','n_sh','Bg1 chaine','ecart col.1','ecart col.2');
for k=1:K
    fprintf('  %6d %12.5f %13.2e %12.2e\n',NSH(k),Bg(k),Bg(k)-pub(1),Bg(k)-pub(2));
end
tol=5e-6;                                    % demi-unite du dernier chiffre
g1=abs(Bg(NSH==1)-pub(1))<tol;
g2=abs(Bg(NSH==2)-pub(2))<tol;
g4=abs(Bg(NSH==4)-pub(2))<tol;
fprintf('\n  colonne 1 <- n_sh = 1 : %s\n',ternstr(g1));
fprintf('  colonne 2 <- n_sh = 2 : %s\n',ternstr(g2));
fprintf('  colonne 2 <- n_sh = 4 : %s  (doit etre NON)\n',ternstr(g4));
GARDE = g1 && g2 && ~g4 && Nu(NSH==1)==5400 && Nu(NSH==2)==6480;
fprintf('\n  GARDE %s\n',ternstr(GARDE,'PASSEE','ECHOUEE'));

%% ======== 4. TOUTES LES OCCURRENCES DE n_sh DANS LE MANUSCRIT =========
fprintf('\n  ---- 4. OCCURRENCES DE n_sh, LIGNE PAR LIGNE ----\n');
hit=find(~cellfun(@isempty,regexp(tl,'n_\{\\mathrm\{sh\}\}\s*=\s*\d','once')));
for i=hit(:).'
    v=regexp(tl{i},'n_\{\\mathrm\{sh\}\}\s*=\s*(\d)','tokens');
    v=cellfun(@(c)str2double(c{1}),v);
    fprintf('  L%-5d n_sh = %-8s | %s\n',i,mat2str(v),strtrim(tl{i}));
end

%% ======== 5. VERDICT ==================================================
fprintf('\n  ---- 5. VERDICT ----\n');
bad=find(~cellfun(@isempty,regexp(tl, ...
    'shoe layering is varied|\$n_\{\\mathrm\{sh\}\}=1\$ and \$n_\{\\mathrm\{sh\}\}=4\$','once')));
fprintf('  phrase de sensibilite (§Methods) :\n');
for i=bad(:).', fprintf('    L%-5d %s\n',i,strtrim(tl{i})); end
fprintf('\n  Correction A APPLIQUER PAR L''AUTEUR (aucun .tex modifie ici) :\n');
fprintf('    fichier : article/ArticleI_DtN_PMSM.tex\n');
for i=bad(:).'
    s=tl{i};
    if contains(s,'=4$')
        fprintf('    ligne %d\n',i);
        fprintf('      avant : %s\n',strtrim(s));
        fprintf('      apres : %s\n',strtrim(strrep(s,'=4$','=2$')));
    end
end
save('R5_nsh.mat','NSH','Nu','Lsv','Bg','pub','dec','GARDE');
fprintf('\n  duree %.0f s\n=== R5 termine ===\n',toc(t0));
diary off;

% ======================================================================
function s=ternstr(c,a,b)
if nargin<2, a='OUI'; end
if nargin<3, b='NON'; end
if c, s=a; else, s=b; end
end
