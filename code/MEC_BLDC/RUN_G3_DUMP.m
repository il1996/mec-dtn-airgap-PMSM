%% RUN_G3_DUMP - vidage cible de deux .mat, pour clore G2-A et G2-B
clear; clc;
diary('G3_dump_out.txt'); diary on;
fprintf('=== G3 : vidage cible ===\n');

fprintf('\n---- A1_table7.mat : les colonnes de la Table 13 ----\n');
A=load('A1_table7.mat');
fn=fieldnames(A);
for i=1:numel(fn)
    v=A.(fn{i});
    fprintf('  %-6s class=%-10s size=%s\n',fn{i},class(v),mat2str(size(v)));
end
if isfield(A,'lab'), fprintf('  lab :\n'); for i=1:numel(A.lab), fprintf('    %2d  %s\n',i,A.lab{i}); end, end
fprintf('\n  Rm (colonnes de modele maille) :\n');
if isfield(A,'Rm'), disp(A.Rm); end
fprintf('\n  VL (colonne lumped) :\n'); if isfield(A,'VL'), fprintf('    %2d  %.15f\n',[1:numel(A.VL);A.VL(:).']); end
fprintf('\n  VF (colonne FEA) :\n');    if isfield(A,'VF'), fprintf('    %2d  %.15f\n',[1:numel(A.VF);A.VF(:).']); end

fprintf('\n---- la matrice 2 x 2 de k_r ----\n');
M=machine_bldc();
rd=@(f)readmatrix(f,'FileType','text','Delimiter','\t','NumHeaderLines',1);
fea=M.FEA.dir;
dmm=rd(fullfile(fea,'magnetostique(Magnetic_loading)','Calculator Expressions Table 1.tab'));
mmf_m_F=dmm(1,2:15); mmf_g_F=dmm(1,16:29);
krF=rd(fullfile(fea,'magnetostique(Magnetic_loading)','Output Variables Table 2.tab')); krF=krF(1,2);
krC=mean(mmf_m_F)/mean(mmf_g_F(2:end));
ikr=16;                                   % k_r est la 16e ligne de la Table 13
if isfield(A,'lab')
    j=find(contains(lower(A.lab),'k_r')); if ~isempty(j), ikr=j(1); end
end
fprintf('  indice de k_r retenu : %d  (%s)\n',ikr,A.lab{ikr});
kr_m1=A.Rm(1).v(ikr); kr_m2=A.Rm(2).v(ikr); kr_lu=A.VL(ikr); kr_fe=A.VF(ikr);
fprintf('  inconnues : n_sh=1 -> %d | n_sh=2 -> %d\n',A.Rm(1).N,A.Rm(2).N);
fprintf('  k_r mesh n_sh=1 : %.15f\n',kr_m1);
fprintf('  k_r mesh n_sh=2 : %.15f\n',kr_m2);
fprintf('  k_r lumped      : %.15f\n',kr_lu);
fprintf('  k_r FEA (table) : %.15f   (relu du projet : %.15f)\n',kr_fe,krF);
fprintf('  reference hors point aberrant : %.15f\n',krC);
fprintf('\n  %-16s %20s %20s\n','colonne','/ %.6f (tous)','/ hors aberrant');
fprintf('  %-16s %19.4f %% %19.4f %%\n','mesh n_sh=1',100*(kr_m1-krF)/krF,100*(kr_m1-krC)/krC);
fprintf('  %-16s %19.4f %% %19.4f %%\n','mesh n_sh=2',100*(kr_m2-krF)/krF,100*(kr_m2-krC)/krC);
fprintf('  %-16s %19.4f %% %19.4f %%\n','lumped',100*(kr_lu-krF)/krF,100*(kr_lu-krC)/krC);

fprintf('\n---- R7_scorecard.mat : lignes de rendement en points ----\n');
R=load('R7_scorecard.mat');
fn=fieldnames(R);
for i=1:numel(fn)
    v=R.(fn{i});
    fprintf('  %-6s class=%-10s size=%s\n',fn{i},class(v),mat2str(size(v)));
end
SC=R.SC;
if iscell(SC)
    fprintf('\n  SC est un cell. Lignes 25 et 29 :\n');
    for r=[25 29]
        fprintf('    ligne %d :',r);
        for c=1:size(SC,2)
            x=SC{r,c};
            if isnumeric(x), fprintf('  %.10f',x); else, fprintf('  %s',string(x)); end
        end
        fprintf('\n');
    end
end
if istable(SC)
    fprintf('\n  SC est une table. Colonnes : %s\n',strjoin(SC.Properties.VariableNames,', '));
    disp(SC([25 29],:));
end
E=R.E(:).';
fprintf('\n  ecarts relatifs #25 et #29 : %+.10f  %+.10f\n',E(25),E(29));

%  recomposer les deux lignes en points depuis SC si possible
mv=[]; rv=[];
if iscell(SC)
    for r=[25 29]
        num=[];
        for c=1:size(SC,2), if isnumeric(SC{r,c}) && isscalar(SC{r,c}), num(end+1)=SC{r,c}; end, end %#ok<AGROW>
        fprintf('    ligne %d, scalaires : ',r); fprintf('%.10f ',num); fprintf('\n');
        if numel(num)>=2, mv(end+1)=num(1); rv(end+1)=num(2); end %#ok<AGROW>
    end
end
if numel(mv)==2
    fprintf('\n  CONVERSION EN POINTS, pleine precision :\n');
    fprintf('    ligne 25 : %.10f - %.10f = %+.10f p.p.\n',mv(1),rv(1),mv(1)-rv(1));
    fprintf('    ligne 29 : %.10f - %.10f = %+.10f p.p.\n',mv(2),rv(2),mv(2)-rv(2));
    Ep=E; Ep(25)=mv(1)-rv(1); Ep(29)=mv(2)-rv(2);
    fprintf('\n  RECOMPTAGE\n');
    fprintf('    dans les 5 %%, tout en relatif        : %d / 29\n',sum(abs(E)<5));
    fprintf('    dans les 5 %%, 25 et 29 en points     : %d / 29\n',sum(abs(Ep)<5));
    fprintf('    a moins d un point du seuil, relatif : ');
    fprintf('#%d(%+.4f) ',[find(abs(E)<5 & abs(E)>4);E(abs(E)<5 & abs(E)>4)]); fprintf('\n');
    fprintf('    a moins d un point du seuil, points  : ');
    fprintf('#%d(%+.4f) ',[find(abs(Ep)<5 & abs(Ep)>4);Ep(abs(Ep)<5 & abs(Ep)>4)]); fprintf('\n');
end

fprintf('\n---- largeurs d encadrement Newton, pleine precision ----\n');
X=load('X1_table5b.mat');
w=X.EN(:,3)-X.EL(:,3);
for i=1:4
    fprintf('  n_sh = %d : lineaire %+.6f %%  Newton %+.6f %%  largeur %.6f points\n', ...
        X.NSH(i),X.EL(i,3),X.EN(i,3),w(i));
end
fprintf('  M_s de la colonne retenue : %d\n',X.Mss(3));
fprintf('\n=== G3 termine ===\n');
diary off;
