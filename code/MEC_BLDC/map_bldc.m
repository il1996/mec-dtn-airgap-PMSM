function map_bldc(ME, S, ax, Bmax, ttl, showdir)
%MAP_BLDC  Carte 2D de |B| sur le maillage raffine, style ANSYS.
%
%   Chaque cellule du maillage porte quatre branches ; on reconstruit en
%   chaque NOEUD un vecteur B en accumulant les branches qui y aboutissent,
%   en separant les RADIALES des TANGENTIELLES -- sans cette separation, le
%   moyennage sur les quatre branches divise chaque composante par deux et
%   la carte sort deux fois trop faible (piege releve dans map_field de
%   MEC_IM). L'entrefer et la couronne d'aimant sont peints avec le champ
%   ANALYTIQUE du DtN etendu, continu en (r,theta).
%
%   Palette et echelle identiques a la carte ANSYS du .docx : jet a 11
%   bandes, 0 a 1.75 T.

if nargin<4||isempty(Bmax), Bmax=1.75; end
if nargin<5, ttl=''; end
if nargin<6||isempty(showdir), showdir=true; end
M=machine_bldc(); Rso=M.Rso; Rsh=31e-3/2;
Ms=ME.Ms; Ls=ME.Ls; re=ME.re; rc=ME.rc; th=ME.ths; dth=ME.dth;
axes(ax); hold on; axis equal off;
cm=jet(11); colormap(ax,cm); caxis(ax,[0 Bmax]);

%% ---------- vecteur B nodal --------------------------------------------
xn=zeros(ME.N,1); yn=xn;
for la=1:Ls
    id=(la-1)*Ms+(1:Ms);
    xn(id)=rc(la)*cos(th); yn(id)=rc(la)*sin(th);
end
id=ME.Ng0+(1:Ms); xn(id)=M.Rsi*cos(th); yn(id)=M.Rsi*sin(th);
dx=xn(ME.b)-xn(ME.a); dy=yn(ME.b)-yn(ME.a); dl=hypot(dx,dy); dl(dl==0)=1;
ux=dx./dl; uy=dy./dl;
rn=hypot(xn(ME.a),yn(ME.a)); rn(rn==0)=1;
isRad=abs(ux.*xn(ME.a)./rn+uy.*yn(ME.a)./rn)>0.7;
Bs=S.B;
acc=@(m)deal(accumarray([ME.a(m);ME.b(m)],[Bs(m).*ux(m);Bs(m).*ux(m)],[ME.N 1]), ...
             accumarray([ME.a(m);ME.b(m)],[Bs(m).*uy(m);Bs(m).*uy(m)],[ME.N 1]), ...
             accumarray([ME.a(m);ME.b(m)],1,[ME.N 1]));
[vxR,vyR,cR]=acc(isRad); [vxT,vyT,cT]=acc(~isRad);
cR(cR==0)=1; cT(cT==0)=1;
vx=vxR./cR+vxT./cT; vy=vyR./cR+vyT./cT; Bn=hypot(vx,vy);

%% ---------- entrefer + aimant : champ analytique du DtN ----------------
nq=1440; tq=linspace(0,2*pi,nq+1); tq(end)=[];
[Brg,Btg]=ME.AG.field(S.Usurf,ME.phi,tq); Bg=hypot(Brg,Btg);
nr=14; rr=linspace(ME.AG.Rro,M.Rsi,nr);
[TT,RR]=meshgrid([tq 2*pi],rr);
surf(RR.*cos(TT),RR.*sin(TT),-ones(size(RR)),repmat([Bg Bg(1)],nr,1),'EdgeColor','none');
%% ---------- ROTOR : aimants et culasse, resolus en (r,theta) ------------
%  AIMANTS. Le flux qui traverse l'aimant se conserve : a travers son
%  epaisseur la section varie en r, donc B_m(r,theta) = B_r(theta)*Rro/r.
%  Entre deux aimants (arc inter-polaire) il ne reste que la frange.
arc=M.embrace*2*pi/M.Nm;
inmag=false(size(tq));
for k=0:M.Nm-1
    inmag=inmag|abs(angle(exp(1i*(tq-(k*2*pi/M.Nm+ME.phi)))))<=arc/2;
end
nm=10; rm=linspace(ME.AG.rmi,ME.AG.Rro,nm);
[TTm,RRm]=meshgrid([tq 2*pi],rm);
Bmg=(abs([Brg Brg(1)]).*[inmag inmag(1)])*ME.AG.Rro./RRm;
Bmg(:,~[inmag inmag(1)])=Bmg(:,~[inmag inmag(1)])*0.25;   % arc inter-polaire
surf(RRm.*cos(TTm),RRm.*sin(TTm),-ones(size(RRm)),Bmg,'EdgeColor','none');

%  CULASSE ROTORIQUE. Elle ne porte PAS un champ uniforme : elle conduit le
%  flux d'un pole vers le suivant, donc son induction suit le FLUX CUMULE
%  entrant par les aimants -- maximale entre deux poles, nulle sous l'axe
%  d'un pole, et alternant 14 fois par tour. C'est cette alternance qui
%  fixe les pertes fer rotoriques et la contrainte de dimensionnement de
%  la culasse ; la peindre d'une seule couleur la faisait disparaitre.
Phry=cumtrapz([tq 2*pi],[Brg Brg(1)])*ME.AG.Rro*M.ls;
Phry=Phry-mean(Phry);
hry=ME.AG.rmi-Rsh;
Bry=abs(Phry)/(hry*M.ls*M.Ki);
nry=12; rry=linspace(Rsh,ME.AG.rmi,nry);
[TTr,RRr]=meshgrid([tq 2*pi],rry);
%  la section offerte au flux tangentiel croit avec r : B decroit vers
%  l'arbre comme le rapport des epaisseurs restantes
wgt=(ME.AG.rmi-RRr)/hry; wgt=max(wgt,0.05);
surf(RRr.*cos(TTr),RRr.*sin(TTr),-ones(size(RRr)),repmat(Bry,nry,1).*wgt, ...
    'EdgeColor','none');
P=sect(0,Rsh,0,2*pi); patch(P(:,1),P(:,2),[.75 .75 .75],'EdgeColor','none');
%  contours des aimants
for k=0:M.Nm-1
    a0=k*2*pi/M.Nm+ME.phi; P=sect(ME.AG.rmi,ME.AG.Rro,a0-arc/2,a0+arc/2);
    plot([P(:,1);P(1,1)],[P(:,2);P(1,2)],'k','LineWidth',.3);
end

%% ---------- cellules du maillage statorique ----------------------------
X=[]; Y=[]; C=[];
for la=1:Ls
    t1=th-dth/2; t2=th+dth/2; r1=re(la); r2=re(la+1);
    X=[X, [r1*cos(t1); r1*cos(t2); r2*cos(t2); r2*cos(t1)]]; %#ok<AGROW>
    Y=[Y, [r1*sin(t1); r1*sin(t2); r2*sin(t2); r2*sin(t1)]]; %#ok<AGROW>
    C=[C, Bn((la-1)*Ms+(1:Ms)).'];                           %#ok<AGROW>
end
patch(X,Y,C,'EdgeColor','none');
tt=linspace(0,2*pi,600);
plot(Rso*cos(tt),Rso*sin(tt),'k','LineWidth',.6);
plot(M.Rsi*cos(tt),M.Rsi*sin(tt),'k','LineWidth',.4);
plot(ME.AG.rmi*cos(tt),ME.AG.rmi*sin(tt),'k','LineWidth',.4);

%% ---------- direction du flux ------------------------------------------
if showdir
    ss=max(1,round(Ms/40));
    idq=[];
    for la=1:Ls, idq=[idq, (la-1)*Ms+(1:ss:Ms)]; end %#ok<AGROW>
    q=0.004./max(Bn(idq),0.4);
    quiver(xn(idq),yn(idq),vx(idq).*q,vy(idq).*q,0,'k','LineWidth',0.3, ...
        'MaxHeadSize',0.5);
end
xlim([-Rso Rso]*1.03); ylim([-Rso Rso]*1.03); view(2);
if ~isempty(ttl), title(ttl,'FontSize',9); end
end

function i=band(v,Bmax), i=max(1,min(11,floor(abs(v)/Bmax*11)+1)); end
function P=sect(ri,ro,t1,t2)
    n=48; ta=linspace(t1,t2,n);
    P=[[ri*cos(ta).' ri*sin(ta).']; [ro*cos(fliplr(ta)).' ro*sin(fliplr(ta)).']];
end
