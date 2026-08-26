function F = magnet_subdomain(M, theta, Nh)
%MAGNET_SUBDOMAIN  Champ a vide d'un SPM par sous-domaines analytiques (DtN).
%
%   F = magnet_subdomain(M, theta, Nh) resout le champ des aimants dans les
%   DEUX couronnes concentriques aimant (rmi<=r<=Rro) + entrefer (Rro<=r<=Rsi),
%   les culasses stator et rotor etant traitees en EQUIPOTENTIELLES (phi=0),
%   harmonique par harmonique (systeme 4x4 exact). C'est l'operateur DtN
%   ETENDU a la couronne d'aimant : la magnetisation radiale est la source,
%   contrairement a un DtN d'entrefer nu qui ne couvre que Rro->Rsi.
%
%   Physique (potentiel scalaire phi, H=-grad phi) :
%     aimant :  B = mu0*mu_r*H + mu0*M_r        -> Laplace: lap(phi)=div(M)/mu_r
%     air    :  B = mu0*H                        -> lap(phi)=0
%   Magnetisation radiale, aimants alternes, arc alpha_p (embrace) :
%     M_r(theta) = sum_{n impair} M_n cos(n p theta),
%     M_n = (4 Br/(mu0 n pi)) sin(n pi alpha_p/2)
%   Ordre spatial mecanique nu = n*p ; rayons normalises par Rsi (evite r^nu).
%
%   Entrees : M (machine_bldc), theta (angles mec, rad), Nh (n impairs max).
%   Sorties : F.Br,F.Bt (au mi-entrefer M.rmid), F.Br_bore (a Rsi^-),
%             F.nu, F.Brn (amplitudes radiales par harmonique au mi-entrefer),
%             F.eval(r) -> [Br,Bt] a un rayon quelconque de l'entrefer.

mu0 = 4*pi*1e-7;
Rs = M.Rsi; Rro = M.Rro; rmi = M.rmi; mur = M.mu_r; p = M.p;
xo = Rro/Rs;  xi = rmi/Rs;               % rayons normalises (<1)
theta = theta(:).';

nn = 1:2:(2*Nh-1);                        % harmoniques impairs de pole
nu = nn*p;                                % ordre spatial mecanique

% amplitudes de magnetisation radiale
Mn = (4*M.Br./(mu0*nn*pi)).*sin(nn*pi*M.embrace/2);

% stockage des coefficients de la couronne d'AIR : phi_A = (a1 x^nu + a2 x^-nu)
A1 = zeros(size(nu)); A2 = A1;
for k = 1:numel(nu)
    v  = nu(k);
    if v*log(1/xi) > 60, continue; end    % harmonique trop haut : x^-v deborde,
                                          % contribution au champ d'entrefer ~0
    Pn = Mn(k)/(mur*(1-v^2));             % coeff. particulier phi_p = Pn*r
    % inconnues X = [a1 a2 c1 c2] ; potentiels normalises en x=r/Rs
    % phi_A = a1 x^v + a2 x^-v                          (air)
    % phi_M = c1 x^v + c2 x^-v + Pn*Rs*x                (aimant)  [Pn*r = Pn*Rs*x]
    S = zeros(4); rhs = zeros(4,1);
    % (1) phi_A(x=1)=0  (culasse stator)
    S(1,:) = [1, 1, 0, 0];                        rhs(1)=0;
    % (2) phi_M(x=xi)=0 (culasse rotor)
    S(2,:) = [0, 0, xi^v, xi^-v];                 rhs(2)= -Pn*Rs*xi;
    % (3) continuite de phi en x=xo
    S(3,:) = [xo^v, xo^-v, -xo^v, -xo^-v];        rhs(3)= Pn*Rs*xo;
    % (4) continuite de Br en x=xo :
    %     -dphiA/dr = -mu_r dphiM/dr + M_n     (Br/mu0 par harmonique)
    %     dphi/dr = (v/Rs)(a1 x^{v-1} - a2 x^{-v-1}) [+ Pn dans l'aimant]
    cA =  (v/Rs)*[xo^(v-1), -xo^(-v-1)];
    cM =  (v/Rs)*[xo^(v-1), -xo^(-v-1)];
    S(4,:) = [-cA(1), -cA(2), mur*cM(1), mur*cM(2)];
    rhs(4) = Mn(k) - mur*Pn;                       % from -(-mur*Pn) moved: see note
    X = S\rhs;
    A1(k)=X(1); A2(k)=X(2);
end

% --- champ dans l'entrefer au rayon r : Br=-mu0 dphi/dr, Bt=-mu0/r dphi/dth
    function [Br,Bt] = eval_air(r)
        x = r/Rs;
        Br = zeros(size(theta)); Bt = Br;
        for kk = 1:numel(nu)
            v = nu(kk);
            dphidr = (v/Rs)*(A1(kk)*x^(v-1) - A2(kk)*x^(-v-1));  % d/dr (coeff)
            phicoef= (A1(kk)*x^v + A2(kk)*x^-v);
            Br = Br - mu0*dphidr*cos(v*theta);
            Bt = Bt + mu0*(v/r)*phicoef*sin(v*theta);
        end
    end

[F.Br, F.Bt]           = eval_air(M.rmid);
[F.Br_bore, ~]         = eval_air(Rs*0.99999);
F.nu = nu; F.theta = theta;
% amplitude radiale par harmonique au mi-entrefer (pour spectre)
x = M.rmid/Rs;
F.Brn = arrayfun(@(kk) -mu0*(nu(kk)/Rs)*(A1(kk)*x^(nu(kk)-1)-A2(kk)*x^(-nu(kk)-1)), 1:numel(nu));
F.eval = @eval_air;
end
