function D = drive_mec(M, the, lam, Leff, opt)
%DRIVE_MEC  Regime TRANSITOIRE en charge de l'entrainement BLDC, pilote MEC.
%
%   Le reseau de reluctances (avec DtN etendu a la couronne d'aimant)
%   fournit la TOTALITE des parametres electromagnetiques. Deux niveaux :
%
%   (a) MODELE LINEAIRE  -- arguments the, lam, Leff :
%         lam(k,:) flux totalise a vide [Wb] contre l'angle ELECTRIQUE the,
%         Leff = L_a - M inductance cyclique. Valide tant que le flux
%         d'induit reste petit devant celui des aimants.
%
%   (b) MODELE SATURABLE -- opt.map (sortie de mec_map) :
%         lambda_k(theta, i_alpha, i_beta) et T(theta, i_alpha, i_beta)
%         tabules par le MEME reseau resolu avec la courbe B(H) reelle.
%       Ce niveau est INDISPENSABLE au demarrage : sous 500 V le courant
%       monte a 23 A, la FMM de bobine atteint 3500 At, la dent depasse
%       2 T et le flux des aimants est CHASSE de la dent -- la FEM
%       s'effondre de 90 % et le courant peut monter deux fois plus haut
%       que ne le prevoit un modele lineaire.
%
%   Le convertisseur et la mecanique sont recopies de BLDC.aedt :
%       netlist 'vf-15-14.ckt' : bus fractionne +/-Vdc/2, pont triphase,
%       6 interrupteurs commandes par la POSITION (PULSE de periode
%       720/Nm deg mec, largeur 720/(3*Nm) -> conduction 120 deg elec),
%       diodes de recuperation, etoile a NEUTRE ISOLE, resistance Rph ;
%       MotionSetup1 : J, frottement visqueux, couple resistant constant.
%   Loi de commutation verifiee sur les sondes de grille FEA (RUN_DIAG_DRIVE):
%   periode mesuree 51.44 deg mec contre 720/Nm = 51.4286 theoriques, calage
%   A+ 0 / C- 60 / B+ 120 / A- 180 / C+ 240 / B- 300 deg elec a 1.6 deg pres.
%
%   MODELE ELECTRIQUE. Etoile a neutre isole :
%       v_k - v_n = R i_k + om * dlam_k/dth_m + sum_j (dlam_k/di_j) di_j/dt
%   avec sum(i_k) = 0. L'etat de chaque bras (+Vdc/2, -Vdc/2 ou flottant)
%   est deduit des grilles ET du signe du courant : une phase non commandee
%   dont le courant n'est pas nul conduit par sa diode de recuperation, ce
%   qui reproduit les pics de commutation.
%
%   COUPLE. Modele lineaire : T = sum_k i_k dlam_k/dth_m + T_detente.
%   Modele saturable : T lu dans la cartographie, ou il est calcule par le
%   TENSEUR DE MAXWELL sur les harmoniques de surface -- seule forme valable
%   quand lambda n'est plus lineaire en i.

p = M.Nm/2;
g = @(f,d) getfield_or(opt,f,d);
dt    = g('dt',1e-6);       tend = g('tend',34.26e-3);
J     = g('J',1e-3);        Bf   = g('Bf',6.07927101854027e-4);
Tload = g('Tload',4.8);     Vdc  = g('Vdc',M.Vdc);
Rph   = g('Rph',M.Rph);     dphi = g('dphi',0);
om0   = g('om0',0);         th0  = g('th0',0);
itol  = g('itol',1e-9);     MAP  = g('map',[]);
s3    = sqrt(3);

%% ---------- tables MEC sur une periode electrique -----------------------
Ntab = 4096; dta = 2*pi/Ntab; thtab = (0:Ntab-1)*dta;
if ~isempty(the)
    the=the(:).';
    if abs(the(end)-2*pi) < 1e-9, the(end)=[]; lam(:,end)=[]; end
    LAM = zeros(3,Ntab);
    for k=1:3
        LAM(k,:) = interp1([the 2*pi],[lam(k,:) lam(k,1)],thtab,'pchip');
    end
    F = fft(LAM,[],2); kk = [0:Ntab/2-1, 0, -Ntab/2+1:-1];
    PSIP = p*real(ifft(1i*kk.*F,[],2));      % dlam/dth_m  [Wb/rad mec]
else
    LAM=zeros(3,Ntab); PSIP=zeros(3,Ntab);
end
if isfield(opt,'Tcog') && ~isempty(opt.Tcog)
    tc = opt.Tcog(:).'; thc = opt.thc(:).';
    if abs(thc(end)-2*pi)<1e-9, thc(end)=[]; tc(end)=[]; end
    TCOG = interp1([thc 2*pi],[tc tc(1)],thtab,'pchip');
else
    TCOG = zeros(1,Ntab);
end

%% ---------- interpolants de la cartographie saturable -------------------
useMap = ~isempty(MAP);
if useMap
    th4=[MAP.the 2*pi]; ia4=MAP.iax;          % periodisation en theta
    wrap=@(X) cat(2,X,X(:,1,:,:));
    Fpsi=cell(1,3); FLa=cell(1,3); Flb=cell(1,3); Flam=cell(1,3);
    for k=1:3
        Fpsi{k}=griddedInterpolant({th4,ia4,ia4},squeeze(wrap(MAP.psi(k,:,:,:))),'linear','nearest');
        FLa{k} =griddedInterpolant({th4,ia4,ia4},squeeze(wrap(MAP.La (k,:,:,:))),'linear','nearest');
        Flb{k} =griddedInterpolant({th4,ia4,ia4},squeeze(wrap(MAP.Lb (k,:,:,:))),'linear','nearest');
        Flam{k}=griddedInterpolant({th4,ia4,ia4},squeeze(wrap(MAP.lam(k,:,:,:))),'linear','nearest');
    end
    Tw=cat(1,MAP.T,MAP.T(1,:,:));
    FT=griddedInterpolant({th4,ia4,ia4},Tw,'linear','nearest');
end

%% ---------- table des grilles (120 deg, calage netlist) -----------------
GU = false(3,Ntab); GD = false(3,Ntab);
de = thtab*180/pi;
inw = @(x,a,b) mod(x-a,360) < mod(b-a,360);
GU(1,:) = inw(de,  0,120);  GD(1,:) = inw(de,180,300);   % A+ / A-
GU(2,:) = inw(de,120,240);  GD(2,:) = inw(de,300, 60);   % B+ / B-
GU(3,:) = inw(de,240,360);  GD(3,:) = inw(de, 60,180);   % C+ / C-

%% ---------- integration --------------------------------------------------
Nt = round(tend/dt)+1;
t  = (0:Nt-1)*dt;
I  = zeros(3,Nt); OM = zeros(1,Nt); TH = zeros(1,Nt);
TE = zeros(1,Nt); IDC = zeros(1,Nt); VT = zeros(3,Nt);
E  = zeros(3,Nt); LAMt = zeros(3,Nt); LEQ = zeros(1,Nt);
UN = zeros(3,Nt); VN = zeros(1,Nt);
i  = zeros(3,1); om = om0; th = th0; half = Vdc/2;

for n = 1:Nt
    thm = mod(p*th,2*pi);                    % angle elec, convention ANSYS
    ja  = idxof(thm,dta,Ntab);
    tq  = mod(thm+dphi,2*pi);                % angle dans le repere MEC

    ial = i(1); ibe = (i(2)-i(3))/s3;
    if useMap
        psi=zeros(3,1); La=zeros(3,1); Lb=zeros(3,1); lamk=zeros(3,1);
        for k=1:3
            psi(k) =Fpsi{k}(tq,ial,ibe); La(k)=FLa{k}(tq,ial,ibe);
            Lb(k)  =Flb{k}(tq,ial,ibe);  lamk(k)=Flam{k}(tq,ial,ibe);
        end
        Tmap=FT(tq,ial,ibe);
    else
        jm=idxof(tq,dta,Ntab);
        psi=PSIP(:,jm); lamk=LAM(:,jm);
        La=[Leff;0;0]; Lb=[0;0;0];           % place tenue, non utilise
    end
    e = om*psi;

    % --- etat des trois bras : interrupteurs puis diodes de roue libre ---
    s = zeros(3,1);
    for k=1:3
        if     GU(k,ja),      s(k)= 1;
        elseif GD(k,ja),      s(k)=-1;
        elseif i(k)> itol,    s(k)=-1;
        elseif i(k)<-itol,    s(k)= 1;
        else,  s(k)=0; i(k)=0;
        end
    end
    C = s~=0; nc = sum(C); v = s*half;
    di = zeros(3,1); Leq = NaN; vn = 0;

    if nc==2
        d=zeros(3,1); d(C)=sign(v(C));       % +1 vers le rail haut
        if ~useMap
            Leq=2*Leff; La=[Leff;Leff;Leff]; Lb=[0;0;0];
        end
        da=d(1); db=(d(2)-d(3))/s3;
        if useMap, Leq=sum(d.*(La*da+Lb*db)); end
        x=(sum(d.*v)-Rph*sum(d.*i)-om*sum(d.*psi))/Leq;
        di=d*x;
        kp=find(d>0);                        % phase reliee au rail haut
        vn=v(kp)-Rph*i(kp)-om*psi(kp)-(La(kp)*da+Lb(kp)*db)*x;
    elseif nc==3
        if ~useMap
            rhs=v-Rph*i-e; vn=mean(rhs); di=(rhs-vn)/Leff;
        else
            P=La+Lb/s3; Q=2*Lb/s3;
            b=v-Rph*i-e;
            z=[P Q ones(3,1)]\b;
            di=[z(1); z(2); -z(1)-z(2)];
            Leq=P(1); vn=-z(3);
        end
    end
    %  tension au NOEUD D'ENROULEMENT (equivalent de NodeVoltage(Ie_x)
    %  d'ANSYS : apres la resistance de phase). Phase conductrice :
    %  v_k - R i_k ; phase flottante : v_n + FEM propre.
    u=v-Rph*i; u(~C)=vn+om*psi(~C);

    if useMap, T=Tmap; else, T=sum(i.*psi)+TCOG(idxof(tq,dta,Ntab)); end
    dom = (T - Bf*om - Tload)/J;

    I(:,n)=i; OM(n)=om; TH(n)=th; TE(n)=T; VT(:,n)=v.*C;
    E(:,n)=e; LAMt(:,n)=lamk; LEQ(n)=Leq; UN(:,n)=u; VN(n)=vn;
    IDC(n)=sum(i(s==1));

    iold = i;
    i  = i + dt*di;
    om = om + dt*dom;
    th = th + dt*om;
    for k=1:3
        if ~GU(k,ja) && ~GD(k,ja) && iold(k)*i(k)<0, i(k)=0; end
    end
    r = sum(i); a = i~=0;
    if any(a), i(a) = i(a) - r/sum(a); end
end

%% ---------- sorties ------------------------------------------------------
D.t=t; D.i=I; D.om=OM; D.n=OM*60/(2*pi); D.th=TH; D.the=mod(p*TH,2*pi);
D.T=TE; D.idc=IDC; D.v=VT; D.e=E; D.lam=LAMt; D.Leq=LEQ;
D.u=UN; D.vn=VN;
%  enveloppe de tension aux bornes, definition ANSYS de 'BEMF' :
%  max(Ie_a,Ie_b,Ie_c) - min(Ie_a,Ie_b,Ie_c)
D.env=max(UN,[],1)-min(UN,[],1);
D.Pcu = Rph*sum(I.^2,1);
D.Pconv = sum(E.*I,1);
D.Pmec  = TE.*OM;
D.Pdc   = Vdc*IDC;
D.Pfric = Bf*OM.^2;
D.dt=dt; D.dphi=dphi; D.useMap=useMap;
end

%% ------------------------------------------------------------------------
function j = idxof(x,dta,Ntab)
j = floor(x/dta)+1; if j>Ntab, j=Ntab; elseif j<1, j=1; end
end

function v = getfield_or(s,f,d)
if isstruct(s) && isfield(s,f) && ~isempty(s.(f)), v=s.(f); else, v=d; end
end
