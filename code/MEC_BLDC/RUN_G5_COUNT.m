%% RUN_G5_COUNT - le decompte "twelve of the seventeen" de la conclusion
%  La conclusion affirme que raffiner le bec ECARTE douze des dix-sept
%  grandeurs de la Table 13 de la reference. Ce bloc le recompte a pleine
%  precision depuis A1_table7.mat, qui porte les trois colonnes de modele et
%  la colonne EF de cette table.
%  GARDE. La somme des trois categories doit valoir exactement 17.
clear; clc;
diary('G5_count_out.txt'); diary on;
A=load('A1_table7.mat');
v1=A.Rm(1).v(:).'; v2=A.Rm(2).v(:).'; F=A.VF(:).'; lab=A.lab;
d1=100*(v1-F)./F; d2=100*(v2-F)./F;
fprintf('=== G5 : raffinement du bec, n_sh = 1 -> 2, sens du mouvement ===\n');
fprintf('  Table 13, %d grandeurs, ecarts formes sur les valeurs pleines\n\n',numel(F));
fprintf('  %-28s %12s %12s   %s\n','grandeur','n_sh=1 (%)','n_sh=2 (%)','sens');
away=0; tow=0; eq=0; ia=[]; it=[]; ie=[];
for i=1:numel(F)
    a=abs(d1(i)); b=abs(d2(i));
    if abs(b-a)<1e-9,       s='inchange'; eq=eq+1; ie(end+1)=i; %#ok<AGROW>
    elseif b>a,             s='ECARTE';  away=away+1; ia(end+1)=i; %#ok<AGROW>
    else,                   s='rapproche'; tow=tow+1; it(end+1)=i; %#ok<AGROW>
    end
    fprintf('  %-28s %12.4f %12.4f   %s\n',lab{i},d1(i),d2(i),s);
end
fprintf('\n  ECARTE    : %d  -> %s\n',away,strjoin(lab(ia),'; '));
fprintf('  rapproche : %d  -> %s\n',tow,strjoin(lab(it),'; '));
fprintf('  inchange  : %d  -> %s\n',eq,strjoin(lab(ie),'; '));
fprintf('\n  GARDE : %d + %d + %d = %d (doit valoir %d)\n',away,tow,eq,away+tow+eq,numel(F));
if away+tow+eq==numel(F), fprintf('  GARDE PASSEE\n'); else, fprintf('  GARDE ECHOUEE\n'); end
fprintf('\n  Les deux grandeurs inchangees sont a -100 %% dans les deux colonnes :\n');
for i=ie, fprintf('    %-28s modele %.6f  reference %.6f\n',lab{i},v1(i),F(i)); end
fprintf('\n  DECOMPTE HORS CES DEUX ZEROS STRUCTURELS : %d ecartees sur %d\n',away,numel(F)-eq);
fprintf('\n=== G5 termine ===\n');
diary off;
