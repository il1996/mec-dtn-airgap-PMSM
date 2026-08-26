%% RUN_MAP2D_MESH - cartes 2D sur le MAILLAGE RAFFINE (Newton + DtN etendu)
clear; clc; close all; t0=tic;
M=machine_bldc(); Ns=M.Ns; p=M.p;
if isfile('kfringe_ident.mat'), S0=load('kfringe_ident.mat'); kfr=S0.kbest; else, kfr=0.325; end
rd=@(f)readmatrix(f,'FileType','text','Delimiter','\t','NumHeaderLines',1);
d4=rd(fullfile(M.FEA.dir,'magnetostique(Magnetic_loading)','Calculator Expressions Plot 4.tab'));
thu=linspace(0,2*pi,3601); thu(end)=[];
BrF=interp1([d4(:,2)*pi/180;2*pi],[d4(:,3);d4(1,3)],thu,'linear','extrap');
amp=@(B,th,k) 2*abs(mean(B(:).'.*exp(-1i*k*th(:).')));

fprintf('=== Cartes 2D sur maillage raffine ===\n');
Ms=540;
%  (1) a vide, non lineaire
ME0=mesh_bldc(M,Ms,4,3,[],0,kfr,'nl');
S0v=solve_bldc_mesh(ME0,1e-8,40);
[Br0,~]=ME0.AG.field(S0v.Usurf,0,thu);
fprintf('  a vide   : Newton %d it, residu %.2e | Bg1 %+.1f %% | rang 8 %+.1f %% | B fer max %.2f T\n', ...
    S0v.iter,S0v.err,100*(amp(Br0,thu,p)-amp(BrF,thu,p))/amp(BrF,thu,p), ...
    100*(amp(Br0,thu,8)-amp(BrF,thu,8))/amp(BrF,thu,8),max(abs(S0v.B(ME0.iron))));
%  (2) en charge : paire a+/b- du six-step
iabc=[1.90;-1.90;0];
ME1=mesh_bldc(M,Ms,4,3,iabc,0,kfr,'nl');
S1v=solve_bldc_mesh(ME1,1e-8,40);
fprintf('  en charge: Newton %d it, residu %.2e | B fer max %.2f T\n', ...
    S1v.iter,S1v.err,max(abs(S1v.B(ME1.iron))));

%% ---------- FIGURE ------------------------------------------------------
f=figure('Color','w','Position',[20 20 1500 560],'Visible','off');
ax1=subplot(1,3,1);
map_bldc(ME0,S0v,ax1,1.75, ...
    sprintf('(a) |B| a vide - %d cellules, Newton',ME0.Ms*ME0.Ls),true);
c=colorbar(ax1,'southoutside'); c.Label.String='|B| (T) - echelle de la carte ANSYS';
ax2=subplot(1,3,2);
map_bldc(ME1,S1v,ax2,1.75, ...
    sprintf('(b) |B| en charge - i_a=%.2f, i_b=%.2f A',iabc(1),iabc(2)),true);
c2=colorbar(ax2,'southoutside'); c2.Label.String='|B| (T)';
%  (c) densite de courant
ax3=subplot(1,3,3); axes(ax3); hold on; axis equal off;
Ssl=0.5*(M.ws1+M.ws2)*M.hs2; taus=2*pi/Ns;
AT=zeros(Ns,3); PH={[1 -2 -15 3 14],[6 -7 -5 8 4],[11 -12 -10 13 9]};
for k=1:3
    for cc=PH{k}
        i=abs(cc); s=sign(cc);
        AT(mod(i-2,Ns)+1,k)=AT(mod(i-2,Ns)+1,k)+s*M.Ntc;
        AT(i,k)            =AT(i,k)            -s*M.Ntc;
    end
end
J=(AT*iabc)/Ssl/1e6; Jm=max(abs(J));
cmJ=[flipud(autumn(128)); winter(128)];
r1=M.Rsi+M.hs0+M.hs1; r2=r1+M.hs2;
for la=1:ME1.Ls
    t1=ME1.ths-ME1.dth/2; t2=ME1.ths+ME1.dth/2;
    ra=ME1.re(la); rb=ME1.re(la+1);
    for cc=1:ME1.Ms
        if ME1.isFe(la,cc), col=[.88 .88 .88];
        else
            k=mod(round(ME1.ths(cc)/taus-0.5),Ns)+1;
            col=cmJ(max(1,min(256,round(128+127*J(k)/Jm))),:);
        end
        P=[ra*cos([t1(cc) t2(cc)]) rb*cos([t2(cc) t1(cc)]); ...
           ra*sin([t1(cc) t2(cc)]) rb*sin([t2(cc) t1(cc)])].';
        patch(P(:,1),P(:,2),col,'EdgeColor','none');
    end
end
P=sect(ME1.AG.rmi,ME1.AG.Rro,0,2*pi); patch(P(:,1),P(:,2),[.75 .85 1],'EdgeColor','none');
P=sect(0,ME1.AG.rmi,0,2*pi); patch(P(:,1),P(:,2),[.8 .8 .8],'EdgeColor','none');
tt=linspace(0,2*pi,600); plot(M.Rso*cos(tt),M.Rso*sin(tt),'k','LineWidth',.6);
colormap(ax3,cmJ); caxis(ax3,[-Jm Jm]); c3=colorbar(ax3,'southoutside');
c3.Label.String='J (A/mm^2)';
xlim([-M.Rso M.Rso]*1.03); ylim([-M.Rso M.Rso]*1.03);
title(sprintf('(c) Densite de courant : %.2f A/mm^2',Jm),'FontSize',9);
sgtitle('BLDC 15/14 - cartes 2D sur MAILLAGE RAFFINE (Newton exact + DtN etendu)');
OUT=fileparts(mfilename('fullpath')); if isempty(OUT), OUT=pwd; end
exportgraphics(f,fullfile(OUT,'BLDC_FIG9_cartes2D_maillage.png'),'Resolution',150);
set(f,'Visible','on'); savefig(f,fullfile(OUT,'BLDC_FIG9_cartes2D_maillage.fig'));
fprintf('\n-> BLDC_FIG9_cartes2D_maillage.fig / .png  (%.0f s)\n',toc(t0));

function P=sect(ri,ro,t1,t2)
    n=200; ta=linspace(t1,t2,n);
    P=[[ri*cos(ta).' ri*sin(ta).']; [ro*cos(fliplr(ta)).' ro*sin(fliplr(ta)).']];
end
