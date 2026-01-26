% EyeLink analysis pipeline (v2, batch, core analyses)
% Assumes getFix_eyelink.m has already processed raw EDF files under
% FreeView_Paper_sup/Processed_data as Dat_Sub*_Ses*.mat with expT appended.
% This script runs ANA preprocessing (build_fixTable) across all subjects,
% and performs core analyses aligned with stat_Nfix_time.m (heatmap/basic
% distributions/angle scans/trial-level gaussian sliding window).

%% User config -----------------------------------------------------------
rootDir    = fileparts(mfilename('fullpath')); % FreeView_Paper_sup
projRoot   = fileparts(rootDir);               % FreeView
processedDir = fullfile(rootDir, 'Processed_data');
outDir     = fullfile(rootDir, 'Results');
maxTrials  = 480;                              % v2 total formal trials
learn_stage_n = 0;                             % v2 no learning stage
outPrefix = 'v2_full';                         % file name prefix for outputs
vers       = {'v2'};                           % label for figure titles
if ~exist(outDir,'dir'); mkdir(outDir); end
if ~exist(processedDir,'dir'); mkdir(processedDir); end

%% Paths for ANA helpers -------------------------------------------------
addpath(genpath(fullfile(projRoot, 'FreeView_Paper', 'Methods', 'function_library_cus', 'ANA')));
addpath(genpath(fullfile(projRoot, 'function_library')));
addpath(genpath(fullfile(projRoot, 'function_library_cus')));

%% Step: Preprocessing + core analyses ----------------------------------
% Global switches (align with stat_Nfix_time; limit to available helpers)
RUN_DATA_PREP = true;
PLOT_BASIC_STATS = true;            % fixation/trial distributions summary
PLOT_HEATMAP = true;                % 1-4s heatmap
PLOT_ANG_Z = false;                  % angle scan in zScore (360/45/90)
PLOT_ANG_PROP = true;               % angle scan in proportion (360/45)
PLOT_TIME_SERIES_GAU = true;        % trial-level gaussian sliding
PLOT_16BIN = true;                 % 16-bin bar plots
PLOT_AXIS_EFFECT = true;           % Axis proportion violin
PLOT_OBLI_EFFECT = true;           % Cardinal vs Oblique violin
PLOT_PIE = false;                   % Pie charts for window effects
SAVE_FIGURES = true;                % write PNGs to outDir

% Analysis parameters (as in stat_Nfix_time)
heat_binSize = 25;
n_bin_FT = 16;
CardColor = [251,  4, 255]/255;
ObliColor = [160, 95, 255]/255;
GapColor  = [57, 198, 255]/255;
AxisColor = mix_RGB_by_HSL(CardColor, ObliColor);
angbinSize = 11.25;      % deg
shift_FT = 360 / n_bin_FT / 2;
edges_FT = linspace(0,360,n_bin_FT+1) - shift_FT;
R_min = 1; R_max = 7.5;  % deg visual angle
win_left = 0; win_right = 4000; % ms

% Subject/session discovery (batch)
datFiles = dir(fullfile(processedDir, 'Dat_Sub*_Ses*.mat'));
if isempty(datFiles)
    error('No Dat_Sub*_Ses*.mat found in %s. Run getFix_eyelink.m first.', processedDir);
end

% Use ANA helper to gather resfiles/sub_ses_res/select_sess from processedDir
dirs = struct();
dirs.mat = processedDir;  % both mat and fix under Processed_data for sup
dirs.fix = processedDir;
[resfiles, sub_ses_res, select_sess] = get_eye_data_files(dirs.fix);
if isempty(resfiles)
    error('get_eye_data_files found nothing under %s.', dirs.fix);
end

if RUN_DATA_PREP
    fprintf('Building fixTable across %d subject-session pairs...\n', numel(select_sess));
    exclude_sub = [3];
    last_trial = maxTrials;
    skip_corr = false;
    [fixTable, Nsubj, img_width, img_height, ut, center, digPlace] = build_fixTable( ...
        select_sess, exclude_sub, sub_ses_res, dirs, resfiles, learn_stage_n, last_trial, skip_corr);
    outFixCsv = fullfile(outDir, sprintf('%s_ALL_fixTable_%dSubj.csv', outPrefix, Nsubj));
    try
        writetable(fixTable, outFixCsv);
        fprintf('Wrote %s\n', outFixCsv);
    catch ME
%         warning('Failed to write fixTable CSV: %s', ME.message);
    end
    % Prepare pairs_FT for angle analysis (match stat_Nfix_time)
    pairs_FT = sub_ses_res(select_sess, :);
    pairs_FT = pairs_FT(~ismember(pairs_FT(:,1), exclude_sub), :);
end

% Prepare cmap16_FT (match stat_Nfix_time)
cmap16_FT = repmat([CardColor; GapColor; ObliColor; GapColor], ceil(n_bin_FT/4), 1);
cmap16_FT = cmap16_FT(1:n_bin_FT,:);

%% Filter + shared variables --------------------------------------------
% filter and basic masks, then run selected analyses
[fixTable, start_FT, dur_FT, angles_FT, xpos_FT, ypos_FT, sub_FT, ses_FT, tri_FT, dnfix_FT, subj_stats] = ...
    filter_fixTable_for_analysis(fixTable, R_max); %#ok<ASGLU>

% Select 1–4s window and radius band
win_select_Fixs = (win_left<=start_FT & start_FT<=win_right) | (win_left<=start_FT+dur_FT & start_FT+dur_FT<=win_right);
win_select_Fixs = win_select_Fixs & (fixTable.r(~fixTable.dropFix)>=R_min & fixTable.r(~fixTable.dropFix) <= R_max);

% cleaned table (non-dropped)
CleanedTable = fixTable(~fixTable.dropFix, :);

%% Basic distributions ---------------------------------------------------
if PLOT_BASIC_STATS
    plot_fixTable_distributions(CleanedTable, img_width, img_height, heat_binSize, ut, center, vers{1});
    set(gcf, 'Name', [vers{1}, ' --fixation统计'], 'NumberTitle', 'off');
    plot_trial_distributions(CleanedTable);
    set(gcf, 'Name', [vers{1}, ' --trial统计'], 'NumberTitle', 'off');
    if SAVE_FIGURES
        saveas(gcf, fullfile(outDir, sprintf('%s_trial_stats.png', outPrefix)));
    end
end

%% Heatmap (1–4s) --------------------------------------------------------
if PLOT_HEATMAP
    plot_fixTable_heatmap(xpos_FT(win_select_Fixs), ypos_FT(win_select_Fixs), heat_binSize, img_width, img_height);
    set(gcf, 'Name', [vers{1}, ' --1-4s热图'], 'NumberTitle', 'off'); title('Fixation Position Heatmap (1–4s)');
    if SAVE_FIGURES
        saveas(gcf, fullfile(outDir, sprintf('%s_heatmap_1_4s.png', outPrefix)));
    end
end

%% ---- 角度扫描曲线 (zScore) ----
if PLOT_ANG_Z && RUN_DATA_PREP % 角度扫描需要预处理数据
    startAngle = -22.5/2;
    angles_FT_ = angles_FT(win_select_Fixs);
    sub_FT_   = sub_FT(win_select_Fixs);
    ses_FT_   = ses_FT(win_select_Fixs);
    % 角度滑动统计未做时间加权（与柱状图/时程曲线不等价）！且边界处理是宽松的，即只要注视点的开始或结束时间在时间窗内即被考虑
    [centers360_FT, m360_FT, se360_FT, subCounts360n_FT] = analyze_angle_curve(angles_FT_, sub_FT_, ses_FT_, pairs_FT, Nsubj, angbinSize, 360, 'zScore',startAngle);
    plot_angle_curve(subCounts360n_FT, centers360_FT, 360, angbinSize, CardColor, GapColor, ObliColor, 'zScore',[-1, 1]);
    set(gcf, 'Name', [vers{1}, ' --全角度扫描'], 'NumberTitle', 'off');
    if SAVE_FIGURES
        saveas(gcf, fullfile(outDir, sprintf('%s_ang_z_360.png', outPrefix)));
    end

    % 45° fold: 将360°的曲线分为8段，每段宽度为45°，对每段内的数据取平均
    % 计算45°fold：将centers360_FT mod 45后相同的点分组，行内平均
    binSize = angbinSize;
    foldPeriod = 45;
    mod_angles = mod(centers360_FT - startAngle, foldPeriod);
    [uniq_mod, ~, ic] = unique(round(mod_angles, 8)); % 防止浮点误差
    nbins45 = numel(uniq_mod);
    subCounts45n_FT = zeros(Nsubj, nbins45);
    for k = 1:nbins45
        idx = (ic == k);
        subCounts45n_FT(:,k) = mean(subCounts360n_FT(:,idx), 2, 'omitnan');
    end
    centers45 = uniq_mod + startAngle;
    plot_angle_curve(subCounts45n_FT, centers45, 45, binSize, CardColor, GapColor, ObliColor, 'zScore',[-.4, .8]);
    set(gcf, 'Name', [vers{1}, ' --45°fold扫描'], 'NumberTitle', 'off');
    if SAVE_FIGURES
        saveas(gcf, fullfile(outDir, sprintf('%s_ang_z_45.png', outPrefix)));
    end

    % 90° fold: 将360°的曲线分为4段，每段宽度为90°，对每段内的数据取平均
    foldPeriod = 90;
    binSize = angbinSize;
    mod_angles90 = mod(centers360_FT - startAngle, foldPeriod);
    [uniq_mod90, ~, ic90] = unique(round(mod_angles90, 8)); % 防止浮点误差
    nbins90 = numel(uniq_mod90);
    subCounts90n_FT = zeros(Nsubj, nbins90);
    for k = 1:nbins90
        idx = (ic90 == k);
        subCounts90n_FT(:,k) = mean(subCounts360n_FT(:,idx), 2, 'omitnan');
    end
    centers90 = uniq_mod90 + startAngle;
    plot_angle_curve(subCounts90n_FT, centers90, 90, binSize, CardColor, GapColor, ObliColor, 'zScore',[-.4, .8]);
    set(gcf, 'Name', [vers{1}, ' --90°fold扫描'], 'NumberTitle', 'off');
    if SAVE_FIGURES
        saveas(gcf, fullfile(outDir, sprintf('%s_ang_z_90.png', outPrefix)));
    end
end

%% ---- 角度占比扫描曲线 (Proportion) ----
if PLOT_ANG_PROP && RUN_DATA_PREP % 角度扫描需要预处理数据
    angles_FT_ = angles_FT(win_select_Fixs);
    sub_FT_   = sub_FT(win_select_Fixs);
    ses_FT_   = ses_FT(win_select_Fixs);
    startAngle = -22.5/2;
    % 角度滑动统计未做时间加权（与柱状图/时程曲线不等价）！且边界处理是宽松的，即只要注视点的开始或结束时间在时间窗内即被考虑
    [centers360_FT, m360_FT, se360_FT, subCounts360n_FT] = analyze_angle_curve(angles_FT_, sub_FT_, ses_FT_, pairs_FT, Nsubj, angbinSize, 360, 'sum1', startAngle);
    plot_angle_curve(subCounts360n_FT, centers360_FT, 360, angbinSize, CardColor, GapColor, ObliColor, 'Proportion',[0.02, 0.04]);
    set(gcf, 'Name', [vers{1}, ' --全角度占比扫描'], 'NumberTitle', 'off');
    if SAVE_FIGURES
        saveas(gcf, fullfile(outDir, sprintf('%s_ang_prop_360.png', outPrefix)));
    end
    % 45° fold: 将360°的曲线分为8段，每段宽度为45°，直接求比例和
    binSize = angbinSize;
    foldPeriod = 45;
    mod_angles = mod(centers360_FT - startAngle, foldPeriod);
    [uniq_mod, ~, ic] = unique(round(mod_angles, 8)); % 防止浮点误差
    nbins45 = numel(uniq_mod);
    subCounts45n_FT = zeros(Nsubj, nbins45);
    for k = 1:nbins45
        idx = (ic == k);
        subCounts45n_FT(:,k) = sum(subCounts360n_FT(:,idx), 2, 'omitnan');
    end
    centers45 = uniq_mod + startAngle;
    plot_angle_curve(subCounts45n_FT, centers45, 45, binSize, CardColor, GapColor, ObliColor, 'Proportion',[0.22, 0.30], false);
    set(gcf, 'Name', [vers{1}, ' --45°fold占比扫描'], 'NumberTitle', 'off');
    if SAVE_FIGURES
        saveas(gcf, fullfile(outDir, sprintf('%s_ang_prop_45.png', outPrefix)));
    end
    % ---------- Fit N°-period cosine with fixed constant = chance ----------
%     % v2版本使用默认统计角度和r刻度标签位置
%     statAngle = 0;
%     rtickTheta = 15;
%     rtickHAlign = 'right';
%     rtickVAlign = 'top';
%     [fit45, fit45_stats] = fit_periodic_cosine_phase(subCounts45n_FT, centers45, 45, 0.25, ...
%         'plot', true, 'histColor', [0.5 0.5 0.5], 'scatterColor', [0 0 0], 'ExpColor', AxisColor, ...
%         'statAngle', statAngle, 'rtickTheta', rtickTheta, 'rtickHAlign', rtickHAlign, 'rtickVAlign', rtickVAlign);
%     set(gcf, 'Name', [vers{1}, ' --45°fold相位vtest'], 'NumberTitle', 'off');
%     if SAVE_FIGURES
%         saveas(gcf, fullfile(outDir, sprintf('%s_ang_prop_45_fit.png', outPrefix)));
%     end
    AngelCurve45.data = subCounts45n_FT;
    AngelCurve45.centers = centers45;
    % 90° fold: 将360°的曲线分为4段，每段宽度为90°，直接求比例和
    foldPeriod = 90;
    binSize = angbinSize;
    mod_angles90 = mod(centers360_FT - startAngle, foldPeriod);
    [uniq_mod90, ~, ic90] = unique(round(mod_angles90, 8)); % 防止浮点误差
    nbins90 = numel(uniq_mod90);
    subCounts90n_FT = zeros(Nsubj, nbins90);
    for k = 1:nbins90
        idx = (ic90 == k);
        subCounts90n_FT(:,k) = sum(subCounts360n_FT(:,idx), 2, 'omitnan');
    end
    centers90 = uniq_mod90 + startAngle;
    plot_angle_curve(subCounts90n_FT, centers90, 90, binSize, CardColor, GapColor, ObliColor, 'Proportion',[0.105, 0.162], false);
    set(gcf, 'Name', [vers{1}, ' --90°fold占比扫描'], 'NumberTitle', 'off');
    if SAVE_FIGURES
        saveas(gcf, fullfile(outDir, sprintf('%s_ang_prop_90.png', outPrefix)));
    end
    AngelCurve90.data = subCounts90n_FT;
    AngelCurve90.centers = centers90;
end

%% ---- 柱状图 / 饼图统计 ----
if PLOT_16BIN || PLOT_AXIS_EFFECT || PLOT_OBLI_EFFECT || PLOT_PIE
    % 统计每个被试在选定时间窗内的注视点角度分布（16-bin计数）
    bin_counts = zeros(Nsubj, n_bin_FT);
    for si = 1:Nsubj
        bin_counts(si,:) = histcounts(mod(angles_FT(win_select_Fixs & sub_FT == pairs_FT(si,1)) + shift_FT, 360) - shift_FT, edges_FT);
    end

    % 16bin柱状图数据
    mean_sectors = normalize_by_dim(bin_counts - 0.5 * (bin_counts(:,[end,1:end-1]) + bin_counts(:,[2:end,1])), 'zScore'); % 环形邻居去趋势 then normalize
    mean_sectors_prop = normalize_by_dim(bin_counts, 'sum1'); % 转为比例
    % Axis/Gaps/Card/Obli效应
    Axis_Effect = sum(bin_counts(:,1:2:16),2) ./ sum(bin_counts,2);   % 占比数据
    Obli_Effect = [sum(bin_counts(:,1:4:16),2), sum(bin_counts(:,3:4:16),2)]./ sum(bin_counts(:,1:2:16),2); % 占比数据
    Obli_Effect_inall = [sum(bin_counts(:,1:4:16),2), sum(bin_counts(:,3:4:16),2)]./ sum(bin_counts,2);       % 占比数据（整体）
    Gap_Mean = sum(bin_counts(:,2:2:16),2) ./ sum(bin_counts,2);
end

% ---- 饼图绘制 ----
if PLOT_PIE
    % axis vs gap
    all_effects = [mean(Obli_Effect_inall(:,1)); mean(Obli_Effect_inall(:,2)); mean(Gap_Mean)];
    plot_pie(all_effects, [CardColor; ObliColor; GapColor], 'ax', [], 'show_text', true, 'legend', {'Card.','Obli.','Gap'});
    set(gcf, 'Name', [vers{1}, ' --选定时窗饼图'], 'NumberTitle', 'off');
    if SAVE_FIGURES; saveas(gcf, fullfile(outDir, sprintf('%s_pie_card_obli_gap.png', outPrefix))); end

    % axis vs gap
    all_effects = [mean(Axis_Effect); mean(Gap_Mean)];
    plot_pie(all_effects, [AxisColor; GapColor], 'ax', [], 'show_text', true,  'legend', {'Axis','Gap'});
    set(gcf, 'Name', [vers{1}, ' --选定时窗饼图（Axis）'], 'NumberTitle', 'off');
    if SAVE_FIGURES; saveas(gcf, fullfile(outDir, sprintf('%s_pie_axis_gap.png', outPrefix))); end

    % card vs obli
    all_effects = [mean(Obli_Effect(:,1)); mean(Obli_Effect(:,2))];
    plot_pie(all_effects, [CardColor; ObliColor], 'ax', [], 'show_text', true, 'legend', {'Card.','Obli.'});
    set(gcf, 'Name', [vers{1}, ' --选定时窗饼图（Card vs Obli）'], 'NumberTitle', 'off');
    if SAVE_FIGURES; saveas(gcf, fullfile(outDir, sprintf('%s_pie_card_obli.png', outPrefix))); end
end

% ---- 16bin分别统计 ----
if PLOT_16BIN
    labels_16bin = arrayfun(@(x) sprintf('%.0f°',x), (edges_FT(1:end-1)+edges_FT(2:end))/2, 'uni', 0);
    plot_bar_multi(mean_sectors, cmap16_FT, labels_16bin, ...
        'ylabel', 'Detrend Z Score', ...
        'xlabel', '', ...
        'title', sprintf('Mean Sector Density (%d–%d ms)', win_left, win_right), ...
        'xtickMode', 'deg', ...
        'showInd', true, ...
        'showIndNum', false);
    set(gcf, 'Name', [vers{1}, ' --16扇区分别统计柱状图'], 'NumberTitle', 'off');
    if SAVE_FIGURES; saveas(gcf, fullfile(outDir, sprintf('%s_bar_16bin_zscore.png', outPrefix))); end

    plot_bar_multi(mean_sectors_prop, cmap16_FT, labels_16bin, ...
        'ylabel', 'Proportion', ...
        'xlabel', '', ...
        'title', sprintf('Mean Sector Density (%d–%d ms)', win_left, win_right), ...
        'xtickMode', 'deg', ...
        'showInd', true, ...
        'showIndNum', false, ...
        'baseline',  1/16, ...
        'baselineLabel', 'Chance');
    set(gcf, 'Name', [vers{1}, ' --16扇区分别统计柱状图（比例）'], 'NumberTitle', 'off');
    if SAVE_FIGURES; saveas(gcf, fullfile(outDir, sprintf('%s_bar_16bin_prop.png', outPrefix))); end
end

% ---- 轴主效应统计 ----
if PLOT_AXIS_EFFECT
    plot_violin_prop(Axis_Effect, AxisColor, {'Axis'}, ...
        'ylabel', 'Proportion', ...
        'xlabel', '', ...
        'title', sprintf('Axis Effect (%d–%d ms)', win_left, win_right), ...
        'showInd', true, ...
        'showIndNum', false, ...
        'showIndLink', false, ...
        'chanceLevel', 0.5);
    set(gcf, 'Name', [vers{1}, ' --AG总效应'], 'NumberTitle', 'off');
    if SAVE_FIGURES; saveas(gcf, fullfile(outDir, sprintf('%s_axis_effect.png', outPrefix))); end
end

% ---- 斜主效应统计 ----
if PLOT_OBLI_EFFECT
    plot_violin_prop(Obli_Effect(:,1), cmap16_FT(1,:), {'Card'}, ...
        'ylabel', 'Proportion', ...
        'xlabel', '', ...
        'title', sprintf('Cardinal Effect (%d–%d ms)', win_left, win_right), ...
        'showInd', true, ...
        'showIndNum', false, ...
        'showIndLink', false, ...
        'chanceLevel', 0.5, ...
        'chanceLabel', 'Chance');
    set(gcf, 'Name', [vers{1}, ' --CO总效应'], 'NumberTitle', 'off');
    if SAVE_FIGURES; saveas(gcf, fullfile(outDir, sprintf('%s_obli_effect_card.png', outPrefix))); end

    plot_violin_prop(Obli_Effect_inall, cmap16_FT([1,3],:), {'Card','Obli'}, ...
        'ylabel', 'Proportion', ...
        'xlabel', '', ...
        'title', sprintf('Cardinal Effect (%d–%d ms)', win_left, win_right), ...
        'showInd', true, ...
        'showIndNum', false, ...
        'showIndLink', true, ...
        'chanceLevel', mean(Gap_Mean)/2, ...
        'chanceLabel', 'Gap/2');
    set(gcf, 'Name', [vers{1}, ' --CO总效应SI'], 'NumberTitle', 'off');
    if SAVE_FIGURES; saveas(gcf, fullfile(outDir, sprintf('%s_obli_effect_si.png', outPrefix))); end
end

%% Trial-level Gaussian sliding window ----------------------------------
if PLOT_TIME_SERIES_GAU
    total_trial = maxTrials - learn_stage_n;
    % Build tailMask using max present trial per subject-session in CleanedTable
    uniq_pairs = unique([CleanedTable.subID, CleanedTable.sesID], 'rows', 'stable');
    tailMask = false(height(uniq_pairs), total_trial);
    % Map from pair index to max trial number observed
    for i = 1:size(uniq_pairs,1)
        si = uniq_pairs(i,1); se = uniq_pairs(i,2);
        tri_of_pair = CleanedTable.tri(CleanedTable.subID==si & CleanedTable.sesID==se);
        if isempty(tri_of_pair)
            continue;
        end
        maxTri = max(tri_of_pair);
        if maxTri < total_trial
            tailMask(i, maxTri+1:total_trial) = true;
        end
    end
    cfg = struct();
    cfg.edges_FT = edges_FT;
    cfg.shift_FT = shift_FT;
    cfg.cmap16_FT = repmat([CardColor; GapColor; ObliColor; GapColor], ceil(n_bin_FT/4), 1);
    cfg.win_left = win_left; cfg.win_right = win_right;
    cfg.step_trials = 1; cfg.total_trial = total_trial;
    cfg.digPlace = digPlace; cfg.mode = 'continuous';
    cfg.doStats = true; cfg.learn_stage_n = learn_stage_n;
    cfg.normMode_FT = 'zScore'; cfg.xlabel = 'Trial number';
    cfg.ylabel = 'Proportion'; cfg.doSmooth = false;
    cfg.ver = vers{1}; cfg.AxisColor = AxisColor;
    cfg.Card_in_Axis = true; cfg.cutoff = 300;
    cfg.tailMask = tailMask;

    GauWin = 40; cfg_gau = cfg; cfg_gau.win_trials = GauWin;
    [Taxis, Tgap, Tcard, Tobli, xWin] = compute_gaussian_window_series(CleanedTable, cfg_gau, total_trial, tailMask);
    plot_sliding_window_analysis(Taxis, Tgap, Tcard, Tobli, xWin, cfg_gau);
    set(gcf, 'Name', [vers{1}, ' --trial高斯滑窗'], 'NumberTitle', 'off');
    if SAVE_FIGURES
        saveas(gcf, fullfile(outDir, sprintf('%s_trial_gaussian.png', outPrefix)));
    end
end

fprintf('Analysis complete. Inputs: %s; Outputs: %s\n', processedDir, outDir);
