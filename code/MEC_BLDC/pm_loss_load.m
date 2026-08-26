function [P,D] = pm_loss_load(M, Rpm, Uarm, ts, iabc, thm, nr, nth)
%PM_LOSS_LOAD  Pertes par courants de Foucault dans les aimants EN CHARGE.
%
%   A vide, seule la modulation de denture cree un dB/dt dans l'aimant
%   (pm_loss). En charge s'ajoute le champ d'INDUIT : avec une alimentation
%   six-step la FMM d'induit est riche en harmoniques d'espace ET de temps
%   qui ne tournent pas avec le rotor et qui balayent donc l'aimant. C'est
%   le terme dominant : la FEA donne 3.18 W en charge contre 0.33 W a vide.
%   Controle prealable : dans BLDC.aedt, 'Eddy Effect'=true porte sur 14
%   objets exactement -- les 14 aimants. SolidLoss est donc bien la perte
%   aimant seule, la comparaison est homogene.
%
%   PRINCIPE. Le potentiel de surface stator est la somme de deux
%   contributions, toutes deux fournies par le meme reseau MEC :
%       U_s(theta_s,t) = U_pm(theta_s ; theta_m(t)) + sum_k i_k(t) U_k(theta_s)
%   U_pm vient du balayage a vide (cogging_mec), U_k du probleme d'induit a
%   1 A (inductance_mec) ; U_k ne depend pas de la position rotor, la
%   geometrie de la couronne etant invariante. Le potentiel vecteur dans
%   l'aimant s'en deduit par le DtN ETENDU, puis
%       J_z = -sigma (dA/dt - <dA/dt>)   (courant net nul par aimant)
%       P   = sigma L INT (dA/dt - <dA/dt>)^2 dS
%
%   DERIVEE ANALYTIQUE (et non par differences finies). A point materiel
%   fixe du rotor, avec Theta = theta_r + theta_m(t),
%       A = r * sum_n c_n(r) [ Pc_n(t) sin(n Theta) - Ps_n(t) cos(n Theta) ]
%       dA/dt = r * sum_n c_n [ (Pc_n' + n om Ps_n) sin(n Theta)
%                             + (n om Pc_n - Ps_n') cos(n Theta) ]
%   Le terme n*omega -- la ROTATION du motif d'induit devant l'aimant --
%   est ainsi exact ; une difference finie sur une grille temporelle grossiere
%   le sous-estimait fortement (les harmoniques d'encoche defilent a
%   Ns*n/60 = 335 Hz, ceux de la commutation dix fois plus vite).
%   En developpant sin(n(theta_r+theta_m)) le facteur de position sort de
%   l'integrale spatiale : le calcul se ramene a DEUX produits matriciels
%   par rayon, ce qui autorise un pas de temps fin (~5 us).

AG=Rpm.AG; nu=AG.nu(:); phis=Rpm.phis; mu0=AG.mu0; mur=AG.mur; Um=AG.Um;
Nt=numel(ts); om=[diff(thm(:).')./diff(ts(:).') 0]; om(end)=om(end-1);

%% ---------- potentiel de surface total a chaque instant -----------------
%  PERIODICITE. Le potentiel de surface vu du STATOR ne se repete qu'au
%  TOUR COMPLET : avec gcd(Ns,Nm) = gcd(15,14) = 1, une rotation d'un pas
%  polaire ou d'un pas d'encoche ne ramene jamais la meme configuration
%  rotor/denture. Replier la position sur un balayage plus court introduit
%  une discontinuite de U_pm, dont la derivee temporelle produit une
%  impulsion parasite d'amplitude 1/dt : la perte calculee devient alors
%  proportionnelle au nombre de pas de temps (erreur constatee x400).
%  On n'accepte donc le repliement QUE si le balayage couvre un tour.
sp=phis(end)-phis(1);
if abs(sp-2*pi)<1e-9
    tq=mod(thm(:).',2*pi);                                  % repliement exact
else
    tq=thm(:).';                                            % doit deja couvrir
    if min(tq)<phis(1)-1e-9 || max(tq)>phis(end)+1e-9
        warning('pm_loss_load:span', ...
            'le balayage rotor ne couvre pas la fenetre simulee');
    end
end
Upm=interp1(phis.',Rpm.Usurf.',tq.','linear','extrap').';   % Nsurf x Nt
%  CHAMP D'INDUIT. Deux formes acceptees :
%   - Uarm matrice Nsurf x 3 : potentiel par ampere, INDEPENDANT de la
%     position (reseau lineaire) ;
%   - Uarm structure .U (Nsurf x 3 x Np) et .phis : potentiel par ampere
%     calcule A CHAQUE POSITION ROTOR avec la permeabilite de dent reelle
%     (arm_field_sat). Les dents etant polarisees par les aimants, cette
%     permeabilite TOURNE avec le rotor : le champ d'induit en est module,
%     ce qui cree dans le repere rotor des harmoniques asynchrones absentes
%     de chacun des deux problemes pris separement.
if isstruct(Uarm)
    tqa=mod(thm(:).',2*pi); Ua=zeros(size(Upm));
    for k=1:3
        Uk=reshape(Uarm.U(:,k,:),[],numel(Uarm.phis));      % Nsurf x Np
        Ua=Ua+interp1(Uarm.phis.',Uk.',tqa.','linear','extrap').'.*iabc(k,:);
    end
    Us=Upm+Ua;
else
    Us=Upm+Uarm*iabc;                                       % superposition
end
Usc=AG.Wc*Us; Uss=AG.Ws*Us;                                 % numax x Nt
alu=AG.alphau(:);
Pc=alu.*Usc; Ps=alu.*Uss;

%% ---------- grille rotor ------------------------------------------------
%  Le champ d'induit et de denture decroit en cosh(n ln(r/rmi))/sinh(n Um)
%  a travers l'aimant : pour les rangs eleves il est confine sur une COUCHE
%  LIMITE contre la surface exterieure. La quadrature radiale doit donc etre
%  raffinee de ce cote (repartition en cosinus) et suffisamment dense.
if nargin<7||isempty(nr),  nr=25;  end
if nargin<8||isempty(nth), nth=720; end
uq=(1-cos(linspace(0,pi/2,nr)));
rq=AG.rmi+(AG.Rro-AG.rmi)*uq/uq(end);
thr=linspace(0,2*pi,nth+1); thr(end)=[]; dth=thr(2)-thr(1);
SINr=sin(nu*thr); COSr=cos(nu*thr);                         % numax x nth
cnm=cos(nu*thm(:).'); snm=sin(nu*thm(:).');                 % numax x Nt

%% ---------- passage AU REPERE ROTOR AVANT de deriver -------------------
%  A = r sum_n c_n [ Pc sin(n(th_r+th_m)) - Ps cos(n(th_r+th_m)) ]
%    = r sum_n c_n [ ac_n(t) sin(n th_r) + bc_n(t) cos(n th_r) ]
%  avec ac = Pc cos(n th_m) - Ps sin(n th_m), bc = Pc sin(n th_m) + Ps cos(n th_m).
%  C'est sur ac et bc -- les coefficients VUS PAR L'AIMANT -- qu'il faut
%  deriver. Deriver Pc, Ps et ajouter separement le terme de rotation n*omega
%  reviendrait a soustraire deux quantites presque egales : la composante
%  synchrone du champ des aimants est ici ~400 fois plus grande que la
%  modulation de denture qui subsiste dans le repere rotor, et la difference
%  de deux termes de cet ordre perd toute precision (erreur constatee x430).
%  Dans le repere rotor cette composante est STATIQUE : elle disparait de la
%  derivee d'elle-meme, sans annulation numerique.
ac=Pc.*cnm-Ps.*snm; bc=Pc.*snm+Ps.*cnm;                     % numax x Nt
dac=gradient(ac,ts(:).',1); dbc=gradient(bc,ts(:).',1);

%% ---------- masques des aimants ----------------------------------------
arc=M.embrace*2*pi/M.Nm; msk=false(M.Nm,nth);
for k=0:M.Nm-1
    msk(k+1,:)=abs(angle(exp(1i*(thr-k*2*pi/M.Nm))))<=arc/2;
end

%% ---------- dA/dt sur toute la couronne ---------------------------------
Gall=zeros(nr,Nt,nth);
for j=1:nr
    r=rq(j); v=log(r/AG.rmi);
    cn=-mu0*mur*(1/r)*cosh(nu*v)./sinh(nu*Um);              % = kk_n/n
    Gall(j,:,:)=r*((cn.*dac).'*SINr + (cn.*dbc).'*COSr);    % Nt x nth = dA/dt
end

%% ---------- courant net nul SUR LA SECTION ENTIERE DE L'AIMANT ---------
%  La condition physique est INT_S J dS = 0 sur toute la section de
%  l'aimant, donc une SEULE constante K par aimant et par instant :
%      K = (INT_S dA/dt dS)/S,   J = -sigma (dA/dt - K).
%  L'imposer rayon par rayon (une constante par couronne elementaire)
%  reviendrait a supposer chaque couronne isolee des autres : cela retire
%  aussi la VARIATION RADIALE de dA/dt, qui est forte ici (le champ de
%  denture decroit en cosh(n ln(r/rmi)) a travers l'aimant), et sous-estime
%  la perte d'un facteur 2 a 4.
wr=zeros(nr,1);                                             % poids de trapz
wr(1)=(rq(2)-rq(1))/2; wr(end)=(rq(end)-rq(end-1))/2;
if nr>2, wr(2:end-1)=(rq(3:end)-rq(1:end-2))/2; end
wS=wr(:).*rq(:)*dth;                                        % poids d'aire
Pt=zeros(1,Nt);
for k=1:M.Nm
    mk=msk(k,:); nk=sum(mk); Sk=sum(wS)*nk;
    num=zeros(1,Nt);
    for j=1:nr, num=num+wS(j)*sum(Gall(j,:,mk),3); end
    K=num/Sk;                                               % 1 x Nt
    for j=1:nr
        g=reshape(Gall(j,:,mk),Nt,nk);
        Pt=Pt+wS(j)*sum((g-K.').^2,2).';
    end
end
Pt=M.sigma_pm*M.ls*Pt;
P=mean(Pt);
D.Pt=Pt; D.ts=ts; D.rq=rq; D.thr=thr;
%  part du seul champ d'induit (aimants eteints dans le potentiel de surface)
D.Uarm_rms=sqrt(mean(sum((Us-Upm).^2,1)));
D.Upm_rms =sqrt(mean(sum(Upm.^2,1)));
end
