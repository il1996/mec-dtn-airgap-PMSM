%% RUN_COG_MESH - couple de detente sur le maillage raffine
%
%  CONTRAINTES D'ECHANTILLONNAGE. La detente d'un 15/14 est d'ordre
%  LCM(15,14) = 210. Deux resolutions doivent suivre :
%    - en ESPACE : la grille de surface porte l'operateur DtN, il lui faut
%      Ms >= 840 pour que l'ordre 210 ne soit pas replie (meme critere que
%      Nsurf dans cogging_mec) ;
%    - en POSITION : le balayage couvre UN PAS POLAIRE (2*pi/Nm), ou la
%      detente accomplit 210/14 = 15 periodes ; 841 positions donnent 56
%      points par periode.
%  Cout : 841 resolutions de Newton. Le demarrage a chaud sur la position
%  precedente les ramene a 2-3 iterations chacune.
clear; clc; t0=tic;
M=machine_bldc(); Ns=M.Ns; Nm=M.Nm; p=M.p; L=M.ls; mu0=4*pi*1e-7;
if isfile('kfringe_ident.mat'), S0=load('kfringe_ident.mat'); kfr=S0.kbest; else, kfr=0.325; end
Ms=900; Np=841; nsh=2;
phis=linspace(0,2*pi/Nm,Np);
T=zeros(1,Np); it=zeros(1,Np); U0=[];
fprintf('=== Detente sur maillage raffine : Ms=%d, %d positions ===\n',Ms,Np);
for q=1:Np
    ME=mesh_bldc(M,Ms,4,3,[],phis(q),kfr,'nl',nsh);
    S=solve_bldc_mesh(ME,1e-10,25,U0); U0=S.U; it(q)=S.iter;
    AG=ME.AG; nu=AG.nu;
    Usc=AG.Wc*S.Usurf; Uss=AG.Ws*S.Usurf;
    cq=cos(nu*phis(q)); sq=sin(nu*phis(q));
    Brc=AG.bru.*Usc+AG.brmq.*cq;  Brs=AG.bru.*Uss+AG.brmq.*sq;
    Btc=-AG.btu.*Uss-AG.btmq.*sq; Bts=AG.btu.*Usc+AG.btmq.*cq;
    T(q)=(L*M.rmid^2/mu0)*pi*sum(Brc.*Btc+Brs.*Bts);
    if mod(q,120)==0
        fprintf('  %4d/%d  (%.0f s, %.1f it/pos)\n',q,Np,toc(t0),mean(it(1:q)));
    end
end
T=T-mean(T);
Y=abs(fft(T(1:end-1))); Y=Y(2:floor((Np-1)/2));
[~,km]=max(Y); ordre=round(km*Nm);
Tpp=(max(T)-min(T))*1e3;

%% ---------- reference : reseau a une dent + FEA transitoire ------------
Rc=cogging_mec(M,1260,0,841,M.muI,kfr);
rd=@(f)readmatrix(f,'FileType','text','Delimiter','\t','NumHeaderLines',1);
dt=rd(fullfile(M.FEA.dir,'transitoire (Back_emf)','Torque Plot.tab'));
th=dt(:,1)*1e-3*(M.speed/60)*2*pi; Tt=dt(:,2);
LCMv=lcm(Ns,Nm);
X=[cos(LCMv*th) sin(LCMv*th) ones(size(th)) th th.^2];
c=X\Tt; r=Tt-X*c; A210=hypot(c(1),c(2)); s210=std(r)*sqrt(2/numel(th));

fprintf('\n  %-26s %11s %11s\n','','maillage','1 dent');
fprintf('  %-26s %11.3f %11.3f  mN.m\n','detente c-c',Tpp,Rc.Tpp);
fprintf('  %-26s %11d %11d\n','ordre dominant',ordre,Rc.order);
fprintf('  %-26s %11d %11s\n','ordre attendu LCM(15,14)',LCMv,'210');
fprintf('\n  reference FEA (transitoire, maillage preserve) :\n');
fprintf('    %.3f +/- %.3f mN.m a l''ordre %d\n',2*A210,4*s210,LCMv);
fprintf('    ecart maillage : %+.3f mN.m | ecart 1 dent : %+.3f mN.m\n', ...
    Tpp-2*A210,Rc.Tpp-2*A210);
fprintf('\n  Newton : %.1f iterations par position en moyenne\n',mean(it));
fprintf('  (%.0f s)\n',toc(t0));
save('cog_mesh.mat','T','phis','Tpp','ordre','it','-v7.3');
