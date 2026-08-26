%% RUN_G4_POINTS - lignes de rendement du Tableau 20 converties en POINTS
%  Le §5.3 du manuscrit impose que les grandeurs qui SONT des pourcentages
%  soient comparees par une difference en POINTS. Deux lignes du Tableau 20
%  y contreviennent. Ce bloc les convertit a pleine precision et recompte
%  tout ce qui en depend.
%  GARDE. Les ecarts relatifs recomposes depuis les deux valeurs doivent
%  redonner exactement le vecteur E sauvegarde par R7.
clear; clc;
diary('G4_points_out.txt'); diary on;
fprintf('=== G4 : rendements en points, et recomptage ===\n');
R=load('R7_scorecard.mat'); S=R.SC; E=R.E(:).';
fprintf('\n  champs de SC : %s\n',strjoin(fieldnames(S).',', '));
n=numel(S);
mv=nan(1,n); rv=nan(1,n); nm=strings(1,n);
fn=fieldnames(S);
for i=1:n
    for k=1:numel(fn)
        x=S(i).(fn{k});
        if ischar(x)||isstring(x), nm(i)=string(x); end
    end
    num=[];
    for k=1:numel(fn)
        x=S(i).(fn{k});
        if isnumeric(x)&&isscalar(x), num(end+1)=x; end %#ok<AGROW>
    end
    if numel(num)>=2, mv(i)=num(1); rv(i)=num(2); end
end
fprintf('\n  GARDE : ecart relatif recompose contre E sauvegarde\n');
Erec=100*(mv-rv)./rv;
g=max(abs(Erec-E));
fprintf('    ecart maximal sur les 29 : %.3e\n',g);
if g<1e-9, fprintf('    GARDE PASSEE\n'); else, fprintf('    GARDE ECHOUEE -- resultat a rejeter\n'); end

fprintf('\n  ---- les deux lignes de rendement ----\n');
for r=[25 29]
    fprintf('    #%d  %-32s\n',r,nm(r));
    fprintf('        modele %.10f   reference %.10f\n',mv(r),rv(r));
    fprintf('        relatif %+.6f %%    POINTS %+.6f p.p.\n',E(r),mv(r)-rv(r));
end

Ep=E; Ep(25)=mv(25)-rv(25); Ep(29)=mv(29)-rv(29);
fprintf('\n  ---- RECOMPTAGE ----\n');
fprintf('    dans les 5, tout en relatif    : %d / 29\n',sum(abs(E)<5));
fprintf('    dans les 5, 25 et 29 en points : %d / 29\n',sum(abs(Ep)<5));
ex=[20 24 27];                                  % les trois declarees non validees
kk=true(1,29); kk(ex)=false;
fprintf('    hors les trois non validees, relatif : %d / 26\n',sum(abs(E(kk))<5));
fprintf('    hors les trois non validees, points  : %d / 26\n',sum(abs(Ep(kk))<5));
fprintf('\n    a moins d un point du seuil de 5 --\n');
fprintf('      convention actuelle : ');
i1=find(abs(E)<5 & abs(E)>4);  for i=i1, fprintf('#%d %s (%+.4f)  ',i,nm(i),E(i)); end
fprintf('\n      lignes 25 et 29 en points : ');
i2=find(abs(Ep)<5 & abs(Ep)>4); for i=i2, fprintf('#%d %s (%+.4f)  ',i,nm(i),Ep(i)); end
fprintf('\n');
fprintf('\n    nombre a moins d un point : %d en relatif, %d en points\n',numel(i1),numel(i2));
fprintf('\n=== G4 termine ===\n');
diary off;
