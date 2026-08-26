%% RUN_SAT - courbe de saturation L_ll(i) par MEC non lineaire
clear; clc;
M=machine_bldc();
if isfile('kfringe_ident.mat'), S=load('kfringe_ident.mat'); kfr=S.kbest; else, kfr=0.75; end
t0=tic;
S=inductance_sat(M,1260,0);
RI=inductance_mec(M,1260,5000,0);
fprintf('=== Inductance de ligne saturee (MEC non lineaire) ===\n');
fprintf('  reference lineaire : 2*L_d = %.2f mH (mu_r = %g)\n',2*RI.Ld*1e3,M.muI);
fprintf('  %8s %10s %10s %10s %10s %10s\n','i (A)','lam (Wb)','L_inc(mH)','L_app(mH)','B_t max','tau=L/2R');
for q=1:numel(S.i)
    fprintf('  %8.3f %10.4f %10.2f %10.2f %10.3f %8.2f ms\n', ...
        S.i(q),S.lam(q),S.L(q)*1e3,S.Lapp(q)*1e3,S.Bt_max(q),1e3*S.L(q)/(2*M.Rph));
end
fprintf('\n  L a 1 A   : %.2f mH  (2*L_d lineaire = %.2f mH)\n', ...
    interp1(S.i,S.L,1)*1e3,2*RI.Ld*1e3);
fprintf('  L a 23 A  : %.2f mH  -> tau = %.2f ms\n', ...
    interp1(S.i,S.L,23)*1e3,1e3*interp1(S.i,S.L,23)/(2*M.Rph));
fprintf('  courant de court-circuit Vdc/2R = %.2f A\n',M.Vdc/(2*M.Rph));
fprintf('  (%.1f s)\n',toc(t0));
