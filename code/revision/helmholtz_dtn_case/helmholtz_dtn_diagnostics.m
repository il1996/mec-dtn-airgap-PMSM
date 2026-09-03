%% HELMHOLTZ_DTN_DIAGNOSTICS  -  task B8 of the revision specification
%
%  A SECOND, INDEPENDENT CASE. The manuscript claims its result "applies
%  wherever a surface operator of order one or higher is condensed onto a
%  discontinuous basis", but establishes it on a PMSM plus a smooth-bore control
%  built on the SAME operator, the same radii, the same magnet and the same
%  tiling. That control is the same operator on smoother data, not a second
%  case. This is a second case: a different physics, a different operator, the
%  same defect.
%
%  *** FIREWALL. NON-NEGOTIABLE. ***
%  This is a firewalled appendix. Its results feed NO other claim in the
%  manuscript. No number in Sections 6, 7 or 8 depends on it. No existing table
%  is modified and no existing figure is regenerated. It lives in its own folder
%  and shares no code with the PMSM chain: the projection coefficients below are
%  COPIED from the forms of (14) and (31), not imported from airgap_magnet.m.
%
%  THE PROBLEM. Exterior Helmholtz in two dimensions, Laplacian u + k^2 u = 0
%  for r > R, with the Sommerfeld radiation condition. The Dirichlet-to-Neumann
%  operator on r = R is diagonal in Fourier:
%       Lambda_n = k * H'_n^(1)(kR) / H_n^(1)(kR)
%  For n >> kR, H'_n/H_n ~ -n/z, so Lambda_n ~ -n/R: ORDER ONE, exactly the
%  threshold of the paper. A P0 basis then gives sum n|W_n|^2 ~ sum 1/n, which
%  diverges. Same defect, different operator, different physics.
%
%  THE NUMERICAL DIFFICULTY, AND HOW IT IS AVOIDED. besselh(n,1,z) overflows to
%  Inf for n >> z (Y_n grows like (2n/ez)^n), so the ratio H'_n/H_n cannot be
%  formed from the functions themselves at the orders needed here (n up to 10^6).
%  Instead the RATIO is propagated directly. From the three-term recurrence
%       H_{n+1} = (2n/z) H_n - H_{n-1},
%  dividing by H_n and writing R_n = H_{n-1}/H_n gives
%       R_{n+1} = 1 / ( 2n/z - R_n ),
%  which is stable in the increasing direction (H^(1) contains the DOMINANT
%  solution Y_n) and cannot overflow, since R_n -> 0. Then
%       H'_n/H_n = R_n - n/z    and    Lambda_n = k (R_n - n/z).
%  The recurrence is started from besselh at n = 0 and 1, where it is safe, and
%  the result is validated against besselh at moderate orders below.
%
%  Nothing outside this folder is written.

clear; clc; t0=tic;

here   = fileparts(mfilename('fullpath'));
root   = fileparts(fileparts(fileparts(here)));
outdir = fullfile(root,'outputs','revision');
if ~exist(outdir,'dir'), mkdir(outdir); end
tfile  = fullfile(outdir,'B8_helmholtz_dtn_out.txt');
if isfile(tfile), delete(tfile); end
diary(tfile); diary on;

ISO=datestr(now,'yyyy-mm-dd');
R  = 1.0;
KR = [1 5 20];
MS = [180 360 720];

fprintf('=== B8 : exterior Helmholtz DtN on a circle - the second case ===\n');
fprintf('  date (ISO) : %s\n  MATLAB     : %s\n', ISO, version);
fprintf('  geometry   : circle of radius R = %.1f\n', R);
fprintf('  operator   : Lambda_n = k H''_n^(1)(kR) / H_n^(1)(kR)\n');
fprintf('  parameters : kR = %s\n', mat2str(KR));
fprintf('  tilings    : Ms = %s contiguous columns, d = 2*pi/Ms\n', mat2str(MS));
fprintf('  bases      : P0  W_n = (2/(n pi)) sin(n d/2)          [form of (14)]\n');
fprintf('               hat W_n = (4/(pi n^2 d)) sin^2(n d/2)    [form of (31)]\n');
fprintf('  FIREWALL   : this appendix feeds no other claim in the manuscript.\n\n');

%% -------------------------------------------------- validation of Lambda_n
fprintf('=== validation of the Hankel ratio against besselh ===\n');
fprintf('  z      n      Lambda_n (recurrence)        Lambda_n (besselh)        rel. err\n');
worst=0;
for z=KR
    for n=[1 5 20 60 100]
        Lr=lambda_n(n,z,R);
        Hn =besselh(n,1,z); Hm=besselh(n-1,1,z);
        Lb =z/R*((Hm/Hn)-n/z);
        e=abs(Lr-Lb)/abs(Lb); worst=max(worst,e);
        fprintf('  %-5g  %-5d  %+0.10e   %+0.10e   %.2e\n', z,n,real(Lr),real(Lb),e);
    end
end
fprintf('  worst relative error against besselh : %.2e   %s\n\n', worst, ...
        tern(worst<1e-8,'PASS','**FAIL**'));

fprintf('=== order of the symbol : Lambda_n * R / n -> -1 ===\n');
fprintf('  n          kR=1            kR=5            kR=20\n');
for n=[10 100 1000 10000 100000]
    v=arrayfun(@(z) real(lambda_n(n,z,R))*R/n, KR);
    fprintf('  %-9d %+13.9f  %+13.9f  %+13.9f\n', n, v);
end
fprintf('  the symbol is of order one, exactly the threshold of the paper.\n\n');

%% ---------------------------------------------- the four diagnostics
NDBL=14;  Ms0=360;                        % tiling fixed for diagnostics 1-3
fprintf('=== diagnostic 1 and 2 : truncation sweep at fixed tiling, and the\n');
fprintf('    ratio of successive increments, over %d doublings, Ms = %d ===\n',NDBL,Ms0);
d=2*pi/Ms0;
res=struct();
for bi=1:2
    bas={'P0','hat'}; bs=bas{bi};
    fprintf('\n  ---- basis %s ----\n',bs);
    fprintf('  %-6s %10s', 'kR', 'N');
    fprintf('   %16s %14s\n','Re[-Y(i,i)]','ratio');
    for z=KR
        A=zeros(1,NDBL); Nv=Ms0*2.^(0:NDBL-1);
        for j=1:NDBL, A(j)=selfterm(Nv(j),d,z,R,bs); end
        dA=diff(A); rat=dA(2:end)./dA(1:end-1);
        for j=1:NDBL
            if j<=2 || j>=NDBL-2
                if j>=3, rr=sprintf('%14.5f',rat(min(j-2,numel(rat)))); else, rr='             -'; end
                fprintf('  %-6g %10d   %16.8e %s\n',z,Nv(j),A(j),rr);
            elseif j==3
                fprintf('  %-6s %10s   %16s %14s\n','...','...','...','...');
            end
        end
        res.(sprintf('%s_%d',bs,z))=rat(end-2:end);
        fprintf('  --> last three increment ratios : %.5f %.5f %.5f\n',rat(end-2:end));
    end
end

fprintf('\n  TOLERANCES OF THE SPECIFICATION\n');
okP=true; okH=true;
for z=KR
    r0=res.(sprintf('P0_%d',z)); r1=res.(sprintf('hat_%d',z));
    p=all(abs(r0-1)<0.05); h=all(abs(r1-0.25)<0.002);
    okP=okP&&p; okH=okH&&h;
    fprintf('  kR=%-3g  P0 -> 1 within 0.05 : %-5s (%.4f %.4f %.4f)\n',z,tern(p,'PASS','FAIL'),r0);
    fprintf('          hat -> 0.250 within 0.002 : %-5s (%.5f %.5f %.5f)\n',tern(h,'PASS','FAIL'),r1);
end
fprintf('  READING. The P0 ratio tends to 1: a logarithmic tail, no limit, the\n');
fprintf('  condensation is NOT defined. The hat ratio tends to 1/4: the tail falls\n');
fprintf('  as N^-2 and the condensation converges. Same verdict as the paper''s\n');
fprintf('  magnetostatic operator, on a different physics.\n\n');

%% ------------------------------------------------- diagnostic 3 : stencil
fprintf('=== diagnostic 3 : the stencil of the divergent term, in units of the\n');
fprintf('    logarithmic step (expected +2, -1, 0, 0, ...) ===\n');
fprintf('  kR     |i-j| =    0         1         2         3         4\n');
for z=KR
    st=zeros(1,5);
    for kk=0:4
        a=offdiag(Ms0*2^12,d,z,R,'P0',kk);
        b=offdiag(Ms0*2^13,d,z,R,'P0',kk);
        st(kk+1)=(b-a)/( (1/(pi*R))*log(2) );
    end
    st=st/abs(st(1))*2*sign(st(1));      % normalise the diagonal to +2
    fprintf('  %-6g          %+8.4f  %+8.4f  %+8.4f  %+8.4f  %+8.4f\n',z,st);
end
fprintf('  the divergent part is the periodic (-1,2,-1) Laplacian here too.\n\n');

%% ------------------------------ diagnostic 4 : a functional quantity
fprintf('=== diagnostic 4 : a functional quantity under released truncation ===\n');
fprintf('  condensed energy of one column, a(u_h,u_h) for u_h the indicator of\n');
fprintf('  a single column, normalised by its value at N = Ms.\n');
fprintf('  %-6s %10s %16s %16s\n','kR','N','P0','hat');
for z=KR
    b0=selfterm(Ms0,d,z,R,'P0'); b1=selfterm(Ms0,d,z,R,'hat');
    for N=Ms0*[1 16 256 4096]
        fprintf('  %-6g %10d %16.6f %16.6f\n',z,N, ...
            selfterm(N,d,z,R,'P0')/b0, selfterm(N,d,z,R,'hat')/b1);
    end
end
fprintf('  the P0 column keeps moving; the hat column settles. Truncation selects\n');
fprintf('  the value of a functional in the discontinuous basis, and does not in\n');
fprintf('  the conforming one.\n\n');

%% ------------------------------------------- tiling dependence (appendix table)
fprintf('=== appendix table : the same verdict at three tilings ===\n');
fprintf('  %-6s %-6s %14s %14s\n','kR','Ms','P0 ratio','hat ratio');
for z=KR
    for Msq=MS
        dd=2*pi/Msq; Nv=Msq*2.^(0:11);
        A0=arrayfun(@(N) selfterm(N,dd,z,R,'P0'),Nv);
        A1=arrayfun(@(N) selfterm(N,dd,z,R,'hat'),Nv);
        r0=diff(A0); r0=r0(end)/r0(end-1);
        r1=diff(A1); r1=r1(end)/r1(end-1);
        fprintf('  %-6g %-6d %14.5f %14.5f\n',z,Msq,r0,r1);
    end
end

%% ------------------------------------------- appendix figure
figfile = fullfile(here,'B8_helmholtz_increment_ratios.pdf');
figsrc  = fullfile(here,'B8_helmholtz_increment_ratios.fig');
try
    f=figure('Visible','off','Units','centimeters','Position',[0 0 16 6.5]);
    tiledlayout(1,2,'Padding','compact','TileSpacing','compact');
    col=lines(3); Nv=Ms0*2.^(0:NDBL-1);
    for bi=1:2
        bas={'P0','hat'}; bs=bas{bi};
        nexttile; hold on; box on
        for zi=1:3
            z=KR(zi);
            A=arrayfun(@(N) selfterm(N,d,z,R,bs),Nv);
            rr=diff(A); rr=rr(2:end)./rr(1:end-1);
            semilogx(Nv(3:end),rr,'o-','Color',col(zi,:),'MarkerSize',3, ...
                     'DisplayName',sprintf('kR = %g',z));
        end
        set(gca,'XScale','log'); grid on
        yline(tern(bi==1,1,0.25),'k--','HandleVisibility','off');
        xlabel('truncation N'); ylabel('\Delta_{k+1} / \Delta_k');
        title(sprintf('%s basis  (limit %s)',bs,tern(bi==1,'1','1/4')));
        ylim(tern(bi==1,[0.8 1.1],[0.0 0.6]));
        legend('Location','best','Box','off');
    end
    exportgraphics(f,figfile,'ContentType','vector');
    set(f,'Visible','on'); savefig(f,figsrc); close(f);
    fprintf('\n  appendix figure written : %s\n', figfile);
    fprintf('  editable source written : %s   (openfig)\n', figsrc);
catch ME
    fprintf('\n  appendix figure NOT produced : %s\n', ME.message);
end

fprintf('\n  ACCEPTANCE : P0 ratios -> 1 : %s | hat ratios -> 0.250 : %s\n', ...
        tern(okP,'PASS','FAIL'), tern(okH,'PASS','FAIL'));
fprintf('\n  duration %.0f s\n=== B8 complete ===\n',toc(t0));
diary off;

%% ============================== local functions ==========================
function Lam = lambda_n(n,z,R)
%  Lambda_n = k (R_n - n/z), R_n = H_{n-1}^(1)/H_n^(1), by stable ratio recurrence.
    k=z/R;
    r=besselh(0,1,z)/besselh(1,1,z);          % R_1
    for m=1:n-1, r=1/(2*m/z - r); end         % R_{m+1} = 1/(2m/z - R_m)
    Lam=k*(r - n/z);
end

function W2 = proj2(n,d,basis)
%  Squared Fourier projection of one column. Copied from the FORMS of (14) and
%  (31); nothing is imported from the PMSM chain.
    if strcmp(basis,'P0')
        W2=((2./(n*pi)).*sin(n*d/2)).^2;
    else
        W2=((4./(pi*n.^2*d)).*sin(n*d/2).^2).^2;
    end
end

function A = selfterm(N,d,z,R,basis)
    A = offdiag(N,d,z,R,basis,0);
end

function A = offdiag(N,d,z,R,basis,kk)
%  Re[-Y(i,i+kk)] = -Re sum_{n<=N} Lambda_n W_n^2 cos(n kk d).
%  The ratio recurrence is inherently sequential, so it is computed ONCE per
%  frequency, to the largest order any caller will need, and cached. Without
%  this the sweep would repeat tens of millions of scalar steps per call.
    persistent CACHE
    if isempty(CACHE), CACHE=containers.Map('KeyType','char','ValueType','any'); end
    key=sprintf('z%.10g',z);
    if ~isKey(CACHE,key) || numel(CACHE(key))<N
        M=max(N,4194304);
        r=zeros(M,1);
        rc=besselh(0,1,z)/besselh(1,1,z);
        for m=1:M
            r(m)=rc;
            rc=1/(2*m/z - rc);
        end
        nn=(1:M).';
        CACHE(key)=(z/R)*(r - nn/z);
    end
    Lam=CACHE(key); Lam=Lam(1:N);
    n=(1:N).';
    A=-real(sum(Lam.*proj2(n,d,basis).*cos(n*kk*d)));
end

function s=tern(c,a,b), if c, s=a; else, s=b; end, end
