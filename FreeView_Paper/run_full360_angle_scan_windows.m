% 2026-03-27 by GitHub Copilot
% Batch-plot non-folded 360-degree angle scan curves for sliding time windows.

clear variables; clear global; fclose('all'); clc

rootDir = fileparts(mfilename('fullpath'));
addpath(fullfile(rootDir, 'Methods', 'function_library_cus', 'ANA'));
addpath(fullfile(rootDir, '..', 'function_library'));

%% ========================= User Config =========================
cfg = struct();
cfg.version = 'v2';
cfg.verDir = '';  % leave empty to use default path pattern: E:\Desktop\临时文件\AttenSamp\FreeView_<version>
cfg.step_ms = 100;
cfg.window_sizes_ms = [500, 1000];
cfg.start_ms = 0;
cfg.end_ms = 5000;
cfg.R_min = 1;
cfg.R_max = 7.5;
cfg.ang_bin_size = 11.25;
cfg.start_angle = -22.5/2;
cfg.norm_mode = 'sum1';
cfg.norm_label = 'Proportion';
cfg.skip_corr = false;
cfg.get_thr_data = true;

% Plot colors (same logic as stat_Nfix_time.m)
CardColor = [61 124 76] / 255;
ObliColor = [69 87 191] / 255;
GapColor  = [121 137 170] / 255;

saveRoot = fullfile(rootDir, 'Results', 'Full360_time_windows', cfg.version);
if ~exist(saveRoot, 'dir')
    mkdir(saveRoot);
end

%% ========================= Prepare Data =========================
[dirs, exclude_sub, learn_stage_n, last_trial] = local_get_version_settings(cfg.version, cfg.verDir, cfg.get_thr_data);

if ~exist(dirs.fix, 'dir')
    error('Fixation directory does not exist: %s', dirs.fix);
end
if ~exist(dirs.mat, 'dir')
    error('Mat directory does not exist: %s', dirs.mat);
end

[resfiles, sub_ses_res, select_sess] = get_eye_data_files(dirs.fix);
pairs_FT = sub_ses_res(select_sess, :);
pairs_FT = pairs_FT(~ismember(pairs_FT(:,1), exclude_sub), :);

[fixTable, Nsubj] = build_fixTable(select_sess, exclude_sub, sub_ses_res, dirs, resfiles, learn_stage_n, last_trial, cfg.skip_corr);
[fixTable, start_FT, dur_FT, angles_FT, ~, ~, sub_FT, ses_FT] = filter_fixTable_for_analysis(fixTable, cfg.R_max);

cleanMask = ~fixTable.dropFix;
r_FT = fixTable.r(cleanMask);

fprintf('Data prepared: version=%s, subjects=%d, cleaned fixations=%d\n', cfg.version, Nsubj, nnz(cleanMask));

%% ========================= Batch Plot =========================
for w = 1:numel(cfg.window_sizes_ms)
    win_ms = cfg.window_sizes_ms(w);
    outDir = fullfile(saveRoot, sprintf('win_%04dms_step_%04dms', win_ms, cfg.step_ms));
    if ~exist(outDir, 'dir')
        mkdir(outDir);
    end

    starts = cfg.start_ms:cfg.step_ms:(cfg.end_ms - win_ms);
    nWins = numel(starts);
    summary = table('Size', [nWins, 4], ...
        'VariableTypes', {'double', 'double', 'double', 'string'}, ...
        'VariableNames', {'start_ms', 'end_ms', 'nFix', 'png_name'});

    fprintf('\nWindow size: %d ms, step: %d ms, number of windows: %d\n', win_ms, cfg.step_ms, nWins);

    for iWin = 1:nWins
        win_left = starts(iWin);
        win_right = win_left + win_ms;

        win_select = (win_left <= start_FT & start_FT <= win_right) | ...
                     (win_left <= start_FT + dur_FT & start_FT + dur_FT <= win_right);
        win_select = win_select & (cfg.R_min <= r_FT & r_FT <= cfg.R_max);

        nFix = nnz(win_select);
        [centers, ~, ~, subjCurves] = analyze_angle_curve(...
            angles_FT(win_select), sub_FT(win_select), ses_FT(win_select), ...
            pairs_FT, Nsubj, cfg.ang_bin_size, 360, cfg.norm_mode, cfg.start_angle);

        plot_angle_curve(subjCurves, centers, 360, cfg.ang_bin_size, CardColor, GapColor, ObliColor, cfg.norm_label, [0.007, 0.04]);
        title(sprintf('%s Full 360 Angle Scan (%d-%d ms), nFix=%d, nSubj=%d', ...
            upper(cfg.version), win_left, win_right, nFix, Nsubj), 'Interpreter', 'none');
        set(gcf, 'Color', 'w', 'Name', sprintf('%s_%d_%dms', cfg.version, win_left, win_right), 'NumberTitle', 'off');

        pngName = sprintf('scan_%04d_%04dms.png', win_left, win_right);
        exportgraphics(gcf, fullfile(outDir, pngName), 'Resolution', 300);
        close(gcf);

        summary.start_ms(iWin) = win_left;
        summary.end_ms(iWin) = win_right;
        summary.nFix(iWin) = nFix;
        summary.png_name(iWin) = string(pngName);

        fprintf('  [%3d/%3d] %d-%d ms, nFix=%d\n', iWin, nWins, win_left, win_right, nFix);
    end

    writetable(summary, fullfile(outDir, sprintf('window_summary_%dms.csv', win_ms)));
    fprintf('Saved %d PNGs to: %s\n', nWins, outDir);
end

fprintf('\nDone. Output root: %s\n', saveRoot);

%% ========================= Local Helpers =========================
function [dirs, exclude_sub, learn_stage_n, last_trial] = local_get_version_settings(ver, customVerDir, getThrData)
    if nargin < 2
        customVerDir = '';
    end
    if nargin < 3
        getThrData = true;
    end

    if isempty(customVerDir)
        verDir = sprintf('E:\\Desktop\\临时文件\\AttenSamp\\FreeView_%s', ver);
    else
        verDir = customVerDir;
    end

    matDir = '\\Data';
    if getThrData
        learn_stage_n = 0;
    else
        learn_stage_n = 72;
    end
    last_trial = 480;

    switch lower(ver)
        case 'v1'
            exclude_sub = 22;
        case 'v1.5'
            exclude_sub = [8, 15];
            if getThrData
                thrDir = '\\Analysis\\Threshold\\Processed_data\\fixDet';
                learn_stage_n = 0;
                last_trial = 480;
            else
                learn_stage_n = 0;
                last_trial = 480 - 72;
            end
        case 'v2'
            exclude_sub = [1, 2, 3, 10];
            matDir = '\\Data\\Formal';
        otherwise
            error('Unsupported version: %s', ver);
    end

    dirs = struct();
    dirs.mat = fullfile(verDir, matDir);
    dirs.fix = fullfile(verDir, 'Analysis', 'Processed_data', 'fixDet');

    if strcmpi(ver, 'v1.5') && getThrData
        dirs.thrDir = fullfile(verDir, thrDir);
    end
end
