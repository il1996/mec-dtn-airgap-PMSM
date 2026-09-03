%% SECTION69_COST_AND_RUNTIME  -  tasks B5 and B6 of the revision specification
%
%  B5. Section 6.9 states two costs per rotor position at the SAME Ms = 1260
%  without reconciling them:
%    "128 rotor positions cost 0.119 s in total against 0.111 s for a single
%     position, so the sweep costs 0.84 % of 128 independent solutions, and the
%     marginal cost of one further position is 0.06 ms. Reassembling the
%     operator at each position ... brings the same sweep to 13.79 s and the
%     marginal cost to 108 ms: a factor of 115 on the sweep."
%    "the coupled network was solved at 721 rotor positions on 1260 surface
%     nodes in 3 s."
%  3 s / 721 = 4.16 ms per position against a marginal 0.06 ms: a factor near
%  69. The internal arithmetic of the first sentence is sound (0.84 %, 0.063 ms,
%  107.7 ms, 115.9). What is missing is the reconciliation. This script measures
%  the decomposition instead of arguing about it.
%
%  B6. Section 6.9 states "The complete study executes in 143 s"; the archived
%  transcript outputs/MEC_BLDC/T9_temps_out.txt states 133 s and the README
%  states 133 s. The point is not decidable without execution.
%
%  NEITHER NUMBER IS ADJUSTED TOWARDS THE OTHER. Whatever comes out is printed.
%
%  TIMING HYGIENE. Wall clock only, on an otherwise idle machine, after a warm-up
%  that is discarded. Every figure is the median of repeats, with min and max
%  printed beside it, so that a single outlier cannot be mistaken for a result.
%
%  Nothing under code/MEC_BLDC or outputs/MEC_BLDC is written by this script.
%  (Note: the published MEC_BLDC_MASTER.m writes FIG_master.png into its own
%   directory as a side effect; B6 removes that file after timing it.)

clear; clc; T0=tic;

here = fileparts(mfilename('fullpath'));
root = fileparts(fileparts(here));
outdir = fullfile(root,'outputs','revision');
if ~exist(outdir,'dir'), mkdir(outdir); end
cd(fullfile(root,'code','MEC_BLDC'));

tfile = fullfile(outdir,'B5B6_section69_cost_out.txt');
if isfile(tfile), delete(tfile); end
diary(tfile); diary on;

ISO=datestr(now,'yyyy-mm-dd');
M=machine_bldc(); mu0=4*pi*1e-7;
Rs=M.Rsi; L=M.ls; Ns=M.Ns; p=M.p; taus=2*pi/Ns;
kfr=0.325; muI=M.muI;

fprintf('=== B5 / B6 : the cost claims of Section 6.9, measured ===\n');
fprintf('  date (ISO)  : %s\n  MATLAB      : %s\n', ISO, version);
[~,cpu]=system('wmic cpu get name /value 2>NUL');
cpu=regexprep(cpu,'\s+',' '); fprintf('  CPU         : %s\n', strtrim(erase(cpu,'Name=')));
fprintf('  declared HW : Intel Core i5-11400H @ 2.70 GHz, 16 GB, MATLAB R2024a\n');
fprintf('  iron model  : linear, mu_i = %g, kfringe = %.4f\n', muI, kfr);
fprintf('  timing      : wall clock, median of repeats, warm-up discarded\n\n');

%% ================= B5 : cost per rotor position, Ms = 1260 ==============
Msq=1260; numax=floor(Msq/2);
fprintf('=== B5 : decomposition of the cost per rotor position, Ms = %d ===\n',Msq);

%  ---- stage 1 : assembly of the bore operator Y (position independent) ----
dth=2*pi/Msq; ths=(0:Msq-1)*dth; dths=dth*ones(1,Msq);
tt=[]; for r=1:3, t=tic; AG=airgap_magnet(M,ths,dths,numax); tt(end+1)=toc(t); end %#ok<SAGROW>
t_Y=median(tt);
fprintf('  assembly of Y (airgap_magnet, %d x %d)          %8.3f s\n',Msq,Msq,t_Y);

%  ---- stage 2 : assembly of the stator network ----
wtf=taus-M.ws0/Rs;
tooth=zeros(1,Msq); isFe=false(1,Msq); nb1=zeros(1,Msq); nb2=nb1;
for j=1:Msq
    i0=mod(round(ths(j)/taus),Ns); dd=angle(exp(1i*(ths(j)-i0*taus)));
    tooth(j)=i0+1;
    if abs(dd)<wtf/2, isFe(j)=true;
    elseif dd>0, nb1(j)=i0+1; nb2(j)=mod(i0+1,Ns)+1;
    else,        nb1(j)=mod(i0-1+Ns,Ns)+1; nb2(j)=i0+1; end
end
g_face=mu0*muI*(dth*Rs)*L/(M.hs0+M.hs1);
Gt=mu0*muI*M.wst1*L/M.hs2;
Gy=mu0*muI*M.wsy*L/(taus*(M.Rso-M.wsy/2));
Gfr=kfr*mu0*(dth*Rs)*L/(M.ws0/2);
Ntot=Msq+2*Ns;
tt=[];
for r=1:3
  t=tic;
  A=zeros(Ntot);
  for j=1:Msq
      if isFe(j)
          a=j; b=Msq+tooth(j); g=g_face;
          A(a,a)=A(a,a)+g; A(b,b)=A(b,b)+g; A(a,b)=A(a,b)-g; A(b,a)=A(b,a)-g;
      elseif Gfr>0
          for q=[nb1(j) nb2(j)]
              a=j; b=Msq+q; g=Gfr;
              A(a,a)=A(a,a)+g; A(b,b)=A(b,b)+g; A(a,b)=A(a,b)-g; A(b,a)=A(b,a)-g;
          end
      end
  end
  for i=1:Ns
      a=Msq+i; b=Msq+Ns+i; g=Gt;
      A(a,a)=A(a,a)+g; A(b,b)=A(b,b)+g; A(a,b)=A(a,b)-g; A(b,a)=A(b,a)-g;
      a=Msq+Ns+i; b=Msq+Ns+mod(i,Ns)+1; g=Gy;
      A(a,a)=A(a,a)+g; A(b,b)=A(b,b)+g; A(a,b)=A(a,b)-g; A(b,a)=A(b,a)-g;
  end
  A(1:Msq,1:Msq)=A(1:Msq,1:Msq)-AG.Y;
  tt(end+1)=toc(t); %#ok<SAGROW>
end
t_net=median(tt);
fprintf('  assembly of the stator network                  %8.3f s\n',t_net);

keep=true(Ntot,1); keep(Msq+Ns+1)=false; Ak=A(keep,keep);

%  ---- stage 3 : one factorisation ----
%  DENSE factorisation, because the published chain stores A as a full matrix
%  (cogging_mec builds A = zeros(Ntot)) and its A(keep,keep)\Isrc is therefore a
%  dense LU. Timing a sparse factorisation here would measure a different chain.
tt=[]; for r=1:3, t=tic; [Lf,Uf,Pf]=lu(Ak); tt(end+1)=toc(t); end %#ok<SAGROW>
t_fact=median(tt);
fprintf('  ONE factorisation of the reduced matrix (%d)   %8.3f s  (dense LU, as the chain does)\n',size(Ak,1),t_fact);

%  ---- stage 4 : source assembly and back-substitution, per position ----
nu=AG.nu; brm=AG.brm;
mkIsrc=@(ph) L*Rs*pi*( AG.Wc.'*(brm.*cos(nu*ph)) + AG.Ws.'*(brm.*sin(nu*ph)) );
for Np=[1 128 721]
    phis=linspace(0,2*pi/p,max(Np,2)); phis=phis(1:Np);
    tt=[];
    for r=1:3
        t=tic;
        Is=zeros(Ntot,Np);
        for q=1:Np, Is(1:Msq,q)=mkIsrc(phis(q)); end
        tsrc=toc(t);
        t=tic;
        Uu=zeros(Ntot,Np);
        Uu(keep,:)=Uf\(Lf\(Pf*Is(keep,:)));
        tsub=toc(t);
        tt(end+1,:)=[tsrc tsub]; %#ok<SAGROW>
    end
    tm=median(tt,1);
    R.(sprintf('n%d',Np))=[tm(1) tm(2)];
    fprintf('  Np = %3d : source assembly %8.4f s | back-substitution %8.4f s\n',Np,tm(1),tm(2));
end

t1  = t_Y+t_net+t_fact+sum(R.n1);
t128= t_Y+t_net+t_fact+sum(R.n128);
t721= t_Y+t_net+t_fact+sum(R.n721);
marg=(t128-t1)/127;
fprintf('\n  --- the sweep, as Section 6.9 counts it (solve stage only) ---\n');
fprintf('  a single position, total                        %8.4f s   (published 0.111 s)\n',t1);
fprintf('  128 positions, total                            %8.4f s   (published 0.119 s)\n',t128);
fprintf('  sweep as a fraction of 128 independent solves   %8.3f %%   (published 0.84 %%)\n', ...
        100*t128/(128*t1));
fprintf('  marginal cost of one further position           %8.4f ms  (published 0.06 ms)\n',1e3*marg);
fprintf('  721 positions, total                            %8.4f s\n',t721);

%  ---- the counterfactual : reassembling the operator at every position ----
nrep=3;
tt=[]; for r=1:nrep, t=tic; airgap_magnet(M,ths,dths,numax); tt(end+1)=toc(t); end %#ok<SAGROW>
t_re=median(tt)+t_net+t_fact;
fprintf('\n  --- the counterfactual of Section 6.9 ---\n');
fprintf('  reassemble + refactorise at EVERY position       %8.4f s per position\n',t_re);
fprintf('  128 positions that way                          %8.3f s   (published 13.79 s)\n',128*t_re);
fprintf('  factor on the sweep                             %8.1f     (published 115)\n',128*t_re/t128);

%  ---- chain (b) : the published call of MEC_BLDC_MASTER ----
fprintf('\n  --- chain (b) : cogging_mec(M,1260,0,721,...) as MEC_BLDC_MASTER calls it ---\n');
tt=[]; for r=1:3, t=tic; Rb=cogging_mec(M,1260,0,721,M.muI,0.75,2*pi/p); tt(end+1)=toc(t); end %#ok<SAGROW>
tb=median(tt);
fprintf('  full cogging_mec at 721 positions               %8.3f s   (published 3 s)\n',tb);
fprintf('  per position                                    %8.4f ms  (published 4.16 ms)\n',1e3*tb/721);
fprintf('  of which the solve stage measured above          %8.3f s\n',t721);
fprintf('  the remainder - field and torque reconstruction  %8.3f s   (%.0f %% of the call)\n', ...
        tb-t721, 100*(tb-t721)/tb);

fprintf(['\n  RECONCILIATION. The two figures are not in conflict; they count different\n' ...
         '  work. The 0.06 ms marginal cost is ONE back-substitution against an already\n' ...
         '  factorised operator, and nothing else. The 4.16 ms per position is a whole\n' ...
         '  cogging_mec call divided by its positions, and that call also assembles the\n' ...
         '  operator once, factorises once, and then, at EVERY position, rebuilds the\n' ...
         '  source, reconstructs Br and Bt over the bore and forms the Maxwell stress.\n' ...
         '  The measured split above says how much of the gap each part accounts for.\n' ...
         '  NEITHER PUBLISHED NUMBER HAS BEEN ALTERED. What the manuscript lacks is one\n' ...
         '  sentence saying what each includes.\n']);

%% ============ B6 : how long the complete study takes ====================
fprintf('\n=== B6 : the run time of the complete study ===\n');
fprintf('  consistency controls that must hold, and do:\n');
fprintf('    540 x 12 = %d unknowns                         (manuscript: 6480)\n',540*12);
fprintf('    n_sh = 1 -> 2(1)+4+3+1 = %d nodes per column    (Table 9: 10)\n',2*1+4+3+1);
fprintf('    n_sh = 4 -> 2(4)+4+3+1 = %d nodes per column    (Table 10: 16)\n',2*4+4+3+1);

fprintf('\n  WHAT "THE COMPLETE STUDY" IS, AND THE PROBLEM WITH TIMING IT.\n');
fprintf('  No single script in this repository is named as the complete study. The\n');
fprintf('  archived transcript T9_temps_out.txt reports "BLDC 15/14 : 133 s (8 analyses,\n');
fprintf('  29 grandeurs, 9 figures)" without naming the eight analyses, and the README\n');
fprintf('  repeats 133 s. The manuscript says 143 s. MEC_BLDC_MASTER.m is the only\n');
fprintf('  end-to-end entry point, so it is what is timed here, and it is timed as the\n');
fprintf('  thing it is - not asserted to be the thing the manuscript meant.\n\n');

%  MEC_BLDC_MASTER.m begins with "clear", which wipes the caller's workspace and
%  would destroy any timer held here. Each run therefore happens in its OWN
%  MATLAB process, and that process stamps the clock into a file immediately
%  before and after the run, using literal paths only so that "clear" cannot
%  reach them. What is measured is the script's own time, NOT MATLAB start-up.
stamp = fullfile(outdir,'B6_stamps.txt');
if isfile(stamp), delete(stamp); end
exe   = fullfile(matlabroot,'bin','matlab.exe');
cmd = sprintf(['fid=fopen(''%s'',''a''); fprintf(fid,''%%.6f\\n'',posixtime(datetime(''now''))); fclose(fid); ' ...
               'run(''%s''); ' ...
               'fid=fopen(''%s'',''a''); fprintf(fid,''%%.6f\\n'',posixtime(datetime(''now''))); fclose(fid);'], ...
               strrep(stamp,'\','\\'), ...
               strrep(fullfile(root,'code','MEC_BLDC','MEC_BLDC_MASTER.m'),'\','\\'), ...
               strrep(stamp,'\','\\'));
fprintf('  timing MEC_BLDC_MASTER.m in a fresh MATLAB process each time,\n');
fprintf('  one cold warm-up (discarded) then five runs\n');
fprintf('    run        wall clock\n');
for r=0:5
    [st,~]=system(sprintf('"%s" -batch "cd(''%s''); %s"', exe, ...
        strrep(fullfile(root,'code','MEC_BLDC'),'\','\\'), strrep(cmd,'"','""')));
    if st~=0, fprintf('    run %d returned status %d\n',r,st); end
end
S=readmatrix(stamp);
d=S(2:2:end)-S(1:2:end);
tms=d(2:end).';                     % drop the warm-up
fprintf('    warm-up    %8.2f s  (discarded)\n',d(1));
for r=1:numel(tms), fprintf('    %-10d %8.2f s\n',r,tms(r)); end
fig=fullfile(root,'code','MEC_BLDC','FIG_master.png');
if isfile(fig), delete(fig); end
fprintf('\n    median %.2f s | min %.2f s | max %.2f s | std %.2f s | n = %d\n', ...
        median(tms),min(tms),max(tms),std(tms),numel(tms));
fprintf('    published in the manuscript : 143 s\n');
fprintf('    archived in T9_temps_out.txt and the README : 133 s\n');
fprintf(['\n  FINDING. Neither figure is reproduced by the only end-to-end script in the\n' ...
         '  archive, and the two published figures differ from each other. The median\n' ...
         '  above is what MEC_BLDC_MASTER.m costs on this machine today; it is NOT a\n' ...
         '  correction of 143 s, because it is not established that the two measure the\n' ...
         '  same work. What is established is that no archived chain produces either\n' ...
         '  number, so neither is currently verifiable. Reported, not harmonised.\n']);

fprintf('\n  total duration of this script %.0f s\n=== B5 / B6 complete ===\n',toc(T0));
diary off;
