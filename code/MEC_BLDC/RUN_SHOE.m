%% RUN_SHOE  -  La FACE DE DENT SEGMENTEE ferme-t-elle le deficit de -28 % ?
%  Hypothese : le fer du bec, suppose infiniment permeable tangentiellement,
%  ne l'est pas ; les cornes, plus chargees, satureraient et gradueraient la
%  frontiere fer/air -> plus de 1ere bande de denture, moins de 2e.
%  Test : Nf = 1 (equipotentiel) -> 9 segments, en lineaire ET en mu(B).
clear; clc;
M=machine_bldc();
D=readmatrix(fullfile(M.FEA.dir,'magnetostique(Magnetic_loading)', ...
    'Calculator Expressions Plot 4.tab'),'FileType','text','Delimiter','\t','NumHeaderLines',1);
thF=D(:,2)*pi/180; BrF=D(:,3); [thF,ix]=sort(thF); BrF=BrF(ix);
fu=@(y,t,k)hypot(2*mean(y(:).*cos(k*t(:))),2*mean(y(:).*sin(k*t(:))));
F.Bg1=fu(BrF,thF,M.p); F.a8=fu(BrF,thF,8); F.a22=fu(BrF,thF,22);
F.a23=fu(BrF,thF,23); F.a37=fu(BrF,thF,37);

K=8; nm=ceil(4*K*pi/(M.ws0/M.Rsi));
fprintf('=== Face de dent SEGMENTEE : bec distribue + saturation locale ===\n');
fprintf('porte-a-faux du bec = %.2f mm de chaque cote ; hauteur de bec %.2f mm\n\n', ...
    (2*pi*M.Rsi/M.Ns - M.ws0 - M.wst1)/2*1e3,(M.hs0+M.hs1)*1e3);
fprintf('%-22s %8s %8s %8s %8s | %7s %7s\n','modele','Bg1','a8','a22','a23','a22/a8','a37/a22');
for Nf=[1 3 5 9]
    for mo={M.muI,'nl'}
        R=subdomain_mec2(M,K,nm,Nf,mo{1},0);
        if ischar(mo{1}), lab=sprintf('Nf=%d  mu(B)',Nf); else, lab=sprintf('Nf=%d  mu=%g',Nf,mo{1}); end
        fprintf('%-22s %8.4f %8.5f %8.5f %8.5f | %7.3f %7.3f\n', ...
            lab,R.Bg1,R.a8,R.a22,R.a23,R.a22/R.a8,R.a37/R.a22);
    end
end
fprintf('%-22s %8.4f %8.5f %8.5f %8.5f | %7.3f %7.3f\n','FEA (ANSYS)', ...
    F.Bg1,F.a8,F.a22,F.a23,F.a22/F.a8,F.a37/F.a22);

% --- etat du bec : y a-t-il vraiment saturation ? ---
R=subdomain_mec2(M,K,nm,9,'nl',0);
fprintf('\n--- Etat magnetique du bec (Nf=9, mu(B)) ---\n');
fprintf('mu_r des branches de bec : min %.0f  max %.0f\n',min(R.muShoe(:,1)),max(R.muShoe(:,1)));
dphi=R.phiF(1,:)-mean(R.phiF(1,:));
fprintf('variation du potentiel le long de la face (dent 1) : %.4g A crete-a-crete\n', ...
    max(dphi)-min(dphi));
fprintf('a comparer a la FMM d''entrefer (~720 A) : %.3f %%\n', ...
    100*(max(dphi)-min(dphi))/720);
fprintf('\n=> si cette variation est <<1 %%, la face est de fait EQUIPOTENTIELLE\n');
fprintf('   et la segmentation ne peut pas modifier le spectre de denture.\n');
