%% TABLE12_FLUXLINEA_REPRINT  -  task B7 of the revision specification
%
%  PROBLEM. The 'flux line A' row of Table 12 prints
%      0.00586 | 0.00587 | 0.00586 | 0.00585
%  to three significant figures, and draws from it the deviations
%      +0.12 % | +0.22 % | +0.09 %.
%  The rounding step of the printed row is 1e-5 / 5.85e-3 = 0.17 %, so the
%  published deviations are FINER THAN THE GRAIN OF THE ROW, and two columns
%  carry the same printed value with two different deviations. A reader cannot
%  check the row, and the row decides the 9-against-8 count between the meshed
%  and lumped iron models (Deliverable 1, item X1).
%
%  WHAT THIS SCRIPT DOES. It reprints the row at six significant figures from
%  the archived full-precision data of the run that produced Table 12, and
%  re-forms the deviations at full precision. THE MODEL IS NOT RE-RUN: the
%  values are read from code/MEC_BLDC/A1_table7.mat, saved by RUN_A1_TABLE7.m
%  (line 134) in the same execution that wrote outputs/MEC_BLDC/A1_table7_out.txt.
%
%  TWO CONTROLS, kept apart because they answer different questions.
%    (1) PROVENANCE. The full-precision deviations must reproduce the figures
%        the declared transcript prints, at that transcript's own precision
%        (0.115 / 0.221 / 0.087 %), to within 0.005 point. If they do not, the
%        row does not come from the declared transcript: REPORT, do not adjust.
%    (2) TRANSCRIPTION. Each two-decimal figure in Table 12 must be the correct
%        rounding of its own full-precision value. Comparing against the
%        manuscript's two-decimal row is NOT a provenance test: a disagreement
%        in the last digit is a rounding fault, not a provenance fault, and
%        conflating the two would misreport a sound number as an unsourced one.
%
%  Nothing under code/MEC_BLDC or outputs/MEC_BLDC is written.

clear; clc;

here = fileparts(mfilename('fullpath'));
root = fileparts(fileparts(here));
outdir = fullfile(root,'outputs','revision');
if ~exist(outdir,'dir'), mkdir(outdir); end
cd(fullfile(root,'code','MEC_BLDC'));

tfile = fullfile(outdir,'B7_fluxlineA_reprint_out.txt');
if isfile(tfile)
    movefile(tfile, fullfile(outdir, sprintf('B7_fluxlineA_reprint_out.superseded_%s.txt', ...
        datestr(now,'yyyymmdd_HHMMSS'))));
end
diary(tfile); diary on;

ISO_DATE = datestr(now,'yyyy-mm-dd');
src = fullfile(root,'code','MEC_BLDC','A1_table7.mat');
d   = dir(src);

fprintf('=== B7 : Table 12, row "flux line A", reprinted at six significant figures ===\n');
fprintf('  date (ISO)   : %s\n', ISO_DATE);
fprintf('  MATLAB       : %s\n', version);
fprintf('  source       : %s\n', 'code/MEC_BLDC/A1_table7.mat');
fprintf('  saved        : %s   (by RUN_A1_TABLE7.m, line 134)\n', datestr(d.datenum,'yyyy-mm-dd HH:MM:SS'));
fprintf('  transcript of the same run : outputs/MEC_BLDC/A1_table7_out.txt\n');
fprintf('  THE MODEL IS NOT RE-RUN. Values are read, not recomputed.\n\n');

S = load(src);                       % Rm, VL, VF, lab
j = find(strcmp(S.lab,'ligne de flux A (Wb/m)'));
if isempty(j), error('row label not found in the archived lab list'); end
fprintf('  row index %d of %d, label "%s"\n\n', j, numel(S.lab), S.lab{j});

v = [S.Rm(1).v(j) S.Rm(2).v(j) S.VL(j) S.VF(j)];
colname = {'mesh n_sh=1','mesh n_sh=2','lumped','reference (FE)'};

fprintf('  ---- as printed in Table 12 (three significant figures) ----\n');
fprintf('    %-15s %-15s %-15s %-15s\n', colname{:});
fprintf('    %-15.5f %-15.5f %-15.5f %-15.5f\n', v);
fprintf('    rounding step of that display : %.3g Wb/m  =  %.3f %% of the row\n\n', ...
        1e-5, 100*1e-5/v(4));

fprintf('  ---- reprinted at six significant figures ----\n');
fprintf('    %-15s %-15s %-15s %-15s\n', colname{:});
fprintf('    %-15.6g %-15.6g %-15.6g %-15.6g\n', v);
fprintf('    full double precision :\n');
for i=1:4, fprintf('      %-16s %.17g\n', colname{i}, v(i)); end
fprintf('\n');

dev  = 100*(v(1:3)-v(4))/v(4);
PUB  = [0.12 0.22 0.09];             % as printed in Table 12
PUBT = [0.115 0.221 0.087];          % as printed in the A1 transcript, 3 decimals

fprintf('  ---- deviations, re-formed at full precision ----\n');
fprintf('    %-15s %12s %12s %12s\n','column','recomputed','Table 12','transcript');
for i=1:3
    fprintf('    %-15s %11.4f %% %10.2f %% %10.3f %%\n', colname{i}, dev(i), PUB(i), PUBT(i));
end
%  TWO SEPARATE CONTROLS. They answer different questions and must not be mixed.
%    (1) PROVENANCE: do the archived full-precision values reproduce the figures
%        the declared transcript prints, at the transcript's own precision?
%    (2) TRANSCRIPTION: does the manuscript's two-decimal row agree with the
%        correctly rounded full-precision value?
gapT = abs(dev-PUBT);
fprintf('\n    (1) PROVENANCE - against the declared transcript (3 decimals)\n');
fprintf('        largest gap : %.4f point (tolerance 0.005)  -> %s\n', max(gapT), ...
        tern(max(gapT)<=0.005,'PASSED','FAILED'));
if max(gapT) <= 0.005
    fprintf('        The row does come from outputs/MEC_BLDC/A1_table7_out.txt. The printed\n');
    fprintf('        row was simply too coarse to let a reader check it.\n');
else
    fprintf('        *** The row does NOT come from the declared transcript. Reported, not adjusted.\n');
end

fprintf('\n    (2) TRANSCRIPTION - against the manuscript row (2 decimals)\n');
fprintf('        %-15s %10s %10s %10s\n','column','full prec.','correct 2dp','printed');
bad = false;
for i=1:3
    corr = round(dev(i),2);
    ok = abs(corr-PUB(i)) < 1e-9;
    bad = bad || ~ok;
    fprintf('        %-15s %9.4f %% %9.2f %% %9.2f %%  %s\n', colname{i}, dev(i), corr, PUB(i), ...
            tern(ok,'ok','MISMATCH'));
end
if bad
    fprintf(['        *** FINDING: at least one printed deviation is not the correct rounding\n' ...
             '        of its own full-precision value. The cause is DOUBLE ROUNDING: the\n' ...
             '        transcript prints 0.115 %%, and 0.115 rounded again to two decimals gives\n' ...
             '        0.12 %%, whereas the underlying 0.1149 %% rounds to 0.11 %%. The value is\n' ...
             '        sound; the last printed digit is not. Corrected below. No value adjusted.\n']);
else
    fprintf('        All printed deviations are the correct rounding of their full-precision value.\n');
end

fprintf('\n  ---- why the coarse row could not carry them ----\n');
vr = round(v,5);
devr = 100*(vr(1:3)-vr(4))/vr(4);
fprintf('    deviations re-formed from the ROUNDED values that Table 12 prints:\n');
for i=1:3
    fprintf('      %-15s %8.3f %%   (published %.2f %%)\n', colname{i}, devr(i), PUB(i));
end
fprintf(['    These do not agree with the published figures, which is the point: the\n' ...
         '    row as printed cannot be used to check its own deviations, and columns 1\n' ...
         '    and 3 share a printed value while carrying different deviations.\n']);

fprintf('\n  ---- replacement row for Table 12 ----\n');
fprintf('    flux line A (Wb/m)   %.6f   %.6f   %.6f   %.6f\n', v);
fprintf('                         %+.2f %%    %+.2f %%    %+.2f %%        --\n', dev);

fprintf('\n=== B7 complete ===\n');
diary off;

function s=tern(c,a,b), if c, s=a; else, s=b; end, end
