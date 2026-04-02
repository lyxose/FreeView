function make_paper_composite_figures(whichSet)
% Make manuscript-ready composite figures by combining subplots directly in MATLAB.
% This script does NOT modify original analysis scripts or stitch existing PNG files.

close all; clc;

if nargin < 1 || isempty(whichSet)
    whichSet = 'all';
end
whichSet = lower(string(whichSet));
fprintf('make_paper_composite_figures start, mode=%s\n', whichSet);

rootDir = fileparts(mfilename('fullpath'));
addpath(fullfile(rootDir, 'Methods', 'function_library_cus', 'ANA'));

outDir = fullfile(rootDir, 'Results', 'Paper_Composites_PDF');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

% Enable extra figures when explicitly requested.
RUN_EXTRA_FIGURES = any(whichSet == ["all", "supp", "supplemental"]); % Figure 4/S1-S5

% User-tunable controls for Figure 2/3/5 layout and display.
ROW_HEIGHT_SCALE = [2/3, 1, 2/3];  % [row1,row2,row3], normalized internally to fill canvas
SHOW_16SECTOR_INDIVIDUAL = true;   % manual switch for individual dots in 16-sector panel
SECTOR_SCATTER_MARKER = 2.0;       % shared marker size for Figure 2 panels B and G
SECTOR_SCATTER_ALPHA = 0.58;       % shared marker alpha for Figure 2 panels B and G
B_PATCH_ALPHA = 1.00;              % Figure 2 panel B background patch alpha (1=opaque)

% -------------------- shared parameters (aligned with stat_Nfix_time.m) --------------------
CardColor = [61 124 76] / 255;
ObliColor = [69 87 191] / 255;
GapColor  = [121 137 170] / 255;
AxisColor = mix_RGB_by_HSL(CardColor, ObliColor);

n_bin_FT = 16;
shift_FT = 360 / n_bin_FT / 2;
edges_FT = linspace(0, 360, n_bin_FT + 1) - shift_FT;
heat_binSize = 25;
angbinSize = 11.25;
keep_Time = 4000;
timeRes_FT = 10;
R_min = 1;
R_max = 7.5;
win_left = 1000;
win_right = 4000;

cmap16_FT = repmat([CardColor; GapColor; ObliColor; GapColor], ceil(n_bin_FT/4), 1);
cmap16_FT = cmap16_FT(1:n_bin_FT,:);

% Version metadata aligned with stat_Nfix_time.m behavior when getThrData=true.
vers = {'v1', 'v1.5', 'v2'};
meta = struct();
meta.v1   = struct('csvRel', fullfile('Results', 'v1',   'ALL_fixTable_26Subj.csv'), 'learn_stage_n', 0, 'last_trial', 480);
meta.v1_5 = struct('csvRel', fullfile('Results', 'v1.5', 'ALL_fixTable_26Subj.csv'), 'learn_stage_n', 0, 'last_trial', 480);
meta.v2   = struct('csvRel', fullfile('Results', 'v2',   'ALL_fixTable_25Subj.csv'), 'learn_stage_n', 0, 'last_trial', 480);

% -------------------- load required versions from existing cleaned csv --------------------
needSupp = any(whichSet == ["all", "supp", "supplemental"]);
needV1 = any(whichSet == ["all", "fig2", "2", "exp1a"]) || needSupp;
needV15 = any(whichSet == ["all", "fig3", "3", "exp1b"]) || needSupp;
needV2 = any(whichSet == ["all", "fig5", "5", "exp2"]) || needSupp;

D = struct();
if needV1
    fprintf('Loading v1 csv...\n');
    csvPath = fullfile(rootDir, meta.v1.csvRel);
    assert(exist(csvPath, 'file') == 2, 'Missing csv: %s', csvPath);
    D.v1 = preprocess_from_fixTable(readtable(csvPath));
    fprintf('Loaded v1 csv.\n');
end
if needV15
    fprintf('Loading v1.5 csv...\n');
    csvPath = fullfile(rootDir, meta.v1_5.csvRel);
    assert(exist(csvPath, 'file') == 2, 'Missing csv: %s', csvPath);
    D.v1_5 = preprocess_from_fixTable(readtable(csvPath));
    fprintf('Loaded v1.5 csv.\n');
end
if needV2
    fprintf('Loading v2 csv...\n');
    csvPath = fullfile(rootDir, meta.v2.csvRel);
    assert(exist(csvPath, 'file') == 2, 'Missing csv: %s', csvPath);
    D.v2 = preprocess_from_fixTable(readtable(csvPath));
    fprintf('Loaded v2 csv.\n');
end

% -------------------- Figure 2/3/5 main composite figures --------------------
verMain = 'v1';
vfMain = strrep(verMain, '.', '_');
if isfield(D, vfMain)
    dm = D.(vfMain);
elseif isfield(D, 'v1_5')
    dm = D.v1_5;
else
    dm = D.v2;
end

if any(whichSet == ["all", "fig2", "2", "exp1a"])
    create_main_effect_figure(D.v1, 'Figure2_Fixation_Distribution_Composite.pdf', true, ...
    CardColor, ObliColor, GapColor, AxisColor, cmap16_FT, edges_FT, heat_binSize, angbinSize, R_max, outDir, 'Figure 2', 'image', ...
    ROW_HEIGHT_SCALE, SHOW_16SECTOR_INDIVIDUAL, SECTOR_SCATTER_MARKER, SECTOR_SCATTER_ALPHA, B_PATCH_ALPHA);
end
if any(whichSet == ["all", "fig3", "3", "exp1b"])
    create_main_effect_figure(D.v1_5, 'Figure3_Fixation_Distribution_Composite.pdf', false, ...
        CardColor, ObliColor, GapColor, AxisColor, cmap16_FT, edges_FT, heat_binSize, angbinSize, R_max, outDir, 'Figure 3', 'image', ...
        ROW_HEIGHT_SCALE, SHOW_16SECTOR_INDIVIDUAL, SECTOR_SCATTER_MARKER, SECTOR_SCATTER_ALPHA, B_PATCH_ALPHA);
end
if any(whichSet == ["all", "fig5", "5", "exp2"])
    create_main_effect_figure(D.v2, 'Figure5_Fixation_Distribution_Composite.pdf', false, ...
    CardColor, ObliColor, GapColor, AxisColor, cmap16_FT, edges_FT, heat_binSize, angbinSize, R_max, outDir, 'Figure 5', 'image', ...
    ROW_HEIGHT_SCALE, SHOW_16SECTOR_INDIVIDUAL, SECTOR_SCATTER_MARKER, SECTOR_SCATTER_ALPHA, B_PATCH_ALPHA);
end

if ~RUN_EXTRA_FIGURES || any(whichSet == ["fig2", "2", "fig3", "3", "fig5", "5", "exp1b", "exp2"])
    fprintf('\nGenerated requested Figure 2/3/5 only.\nOutput: %s\n', outDir);
    return;
end

% -------------------- Figure 4 (trial sequence, EXP1b / v1.5) --------------------
d4 = dm;
if isfield(D, 'v1_5')
    d4 = D.v1_5;
end
f4 = figure('Color', 'w', 'Units', 'normalized', 'Position', [0.21 0.19 0.64*(2/3) 0.68*(2/3)]);
tl4 = tiledlayout(f4, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

cfg_win = build_trial_cfg(d4, edges_FT, shift_FT, cmap16_FT, win_left, win_right, AxisColor, meta.v1_5.learn_stage_n, meta.v1_5.last_trial);
cfg_gau = cfg_win;
cfg_gau.win_trials = 40;
[TseriesAxisG, TseriesGapG, TseriesCardG, TseriesObliG, xWinG] = compute_gaussian_window_series(d4.visTable, cfg_gau, cfg_win.total_trial); %#ok<ASGLU>

ax = nexttile(tl4, 1);
cfg_axis = struct('doStats', true, 'xlabel', 'Trial number', 'ylabel', 'Proportion', ...
    'chanceLevel', 0.5, 'chanceLabel', 'Chance', 'statTail', 'both', 'axesHandle', ax);
plot_single_prop(TseriesAxisG, xWinG, AxisColor, 'Axis proportion', cfg_axis);
title(ax, 'Axis Effect Across Trial Sequence');
add_panel_label(ax, 'A');

ax = nexttile(tl4, 2);
cfg_card = struct('doStats', true, 'xlabel', 'Trial number', 'ylabel', 'Proportion', ...
    'chanceLevel', 0.5, 'chanceLabel', 'Chance', 'statTail', 'right', 'axesHandle', ax);
plot_single_prop(TseriesCardG, xWinG, CardColor, 'Card proportion', cfg_card);
title(ax, 'Cardinal Effect Across Trial Sequence');
add_panel_label(ax, 'B');

remove_legends_from_figure(f4);

set(f4, 'Renderer', 'painters');
enforce_uniform_fontsize(f4, 11);
exportgraphics(f4, fullfile(outDir, 'Figure4_Trial_Sequence_Composite.pdf'), 'ContentType', 'vector');
exportgraphics(f4, fullfile(outDir, 'Figure4_Trial_Sequence_Composite.png'), 'Resolution', 600);
savefig(f4, fullfile(outDir, 'Figure4_Trial_Sequence_Composite.fig'));
close(f4);

% -------------------- Supplemental Figure S3 (trial sequence, EXP1a / v1) --------------------
dS3 = d4;
if isfield(D, 'v1')
    dS3 = D.v1;
end
fS3 = figure('Color', 'w', 'Units', 'normalized', 'Position', [0.21 0.19 0.64*(2/3) 0.68*(2/3)]);
tlS3 = tiledlayout(fS3, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

cfg_s3 = build_trial_cfg(dS3, edges_FT, shift_FT, cmap16_FT, win_left, win_right, AxisColor, meta.v1.learn_stage_n, meta.v1.last_trial);
cfg_s3g = cfg_s3;
cfg_s3g.win_trials = 40;
[S3_axis, ~, S3_card, ~, xWinS3] = compute_gaussian_window_series(dS3.visTable, cfg_s3g, cfg_s3.total_trial);

ax = nexttile(tlS3, 1);
cfg_axisS3 = struct('doStats', true, 'xlabel', 'Trial number', 'ylabel', 'Proportion', ...
    'chanceLevel', 0.5, 'chanceLabel', 'Chance', 'statTail', 'both', 'axesHandle', ax);
plot_single_prop(S3_axis, xWinS3, AxisColor, 'Axis proportion', cfg_axisS3);
title(ax, 'Axis Effect Across Trial Sequence');
add_panel_label(ax, 'A');

ax = nexttile(tlS3, 2);
cfg_cardS3 = struct('doStats', true, 'xlabel', 'Trial number', 'ylabel', 'Proportion', ...
    'chanceLevel', 0.5, 'chanceLabel', 'Chance', 'statTail', 'right', 'axesHandle', ax);
plot_single_prop(S3_card, xWinS3, CardColor, 'Card proportion', cfg_cardS3);
title(ax, 'Cardinal Effect Across Trial Sequence');
add_panel_label(ax, 'B');

    remove_legends_from_figure(fS3);

set(fS3, 'Renderer', 'painters');
enforce_uniform_fontsize(fS3, 11);
exportgraphics(fS3, fullfile(outDir, 'Supplemental_Figure_S3_Trial_Sequence_exp1a.pdf'), 'ContentType', 'vector');
exportgraphics(fS3, fullfile(outDir, 'Supplemental_Figure_S3_Trial_Sequence_exp1a.png'), 'Resolution', 600);
savefig(fS3, fullfile(outDir, 'Supplemental_Figure_S3_Trial_Sequence_exp1a.fig'));
close(fS3);

% -------------------- Figure 6 (trial sequence, EXP2; reuse Figure 4 logic) --------------------
if isfield(D, 'v2')
    d6 = D.v2;
    f6 = figure('Color', 'w', 'Units', 'normalized', 'Position', [0.21 0.19 0.64*(2/3) 0.68*(2/3)]);
    tl6 = tiledlayout(f6, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

    cfg_win6 = build_trial_cfg(d6, edges_FT, shift_FT, cmap16_FT, win_left, win_right, AxisColor, meta.v2.learn_stage_n, meta.v2.last_trial);
    cfg_gau6 = cfg_win6;
    cfg_gau6.win_trials = 40;
    [TseriesAxis6, ~, TseriesCard6, ~, xWin6] = compute_gaussian_window_series(d6.visTable, cfg_gau6, cfg_win6.total_trial);

    ax = nexttile(tl6, 1);
    cfg_axis6 = struct('doStats', true, 'xlabel', 'Trial number', 'ylabel', 'Proportion', ...
        'chanceLevel', 0.5, 'chanceLabel', 'Chance', 'statTail', 'both', 'axesHandle', ax);
    plot_single_prop(TseriesAxis6, xWin6, AxisColor, 'Axis proportion', cfg_axis6);
    title(ax, 'Axis Effect Across Trial Sequence (EXP2)');
    add_panel_label(ax, 'A');

    ax = nexttile(tl6, 2);
    cfg_card6 = struct('doStats', true, 'xlabel', 'Trial number', 'ylabel', 'Proportion', ...
        'chanceLevel', 0.5, 'chanceLabel', 'Chance', 'statTail', 'right', 'axesHandle', ax);
    plot_single_prop(TseriesCard6, xWin6, CardColor, 'Card proportion', cfg_card6);
    title(ax, 'Cardinal Effect Across Trial Sequence (EXP2)');
    add_panel_label(ax, 'B');

    remove_legends_from_figure(f6);

    set(f6, 'Renderer', 'painters');
    enforce_uniform_fontsize(f6, 11);
    exportgraphics(f6, fullfile(outDir, 'Figure6_Trial_Sequence_Composite.pdf'), 'ContentType', 'vector');
    exportgraphics(f6, fullfile(outDir, 'Figure6_Trial_Sequence_Composite.png'), 'Resolution', 600);
    savefig(f6, fullfile(outDir, 'Figure6_Trial_Sequence_Composite.fig'));
    close(f6);
end

% -------------------- Figure 7 (cross-version comparisons) --------------------
f7 = figure('Color', 'w', 'Units', 'normalized', 'Position', [0.16 0.24 0.78*(2/3) 0.64*(2/3)]);
tl7 = tiledlayout(f7, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

compare_groups = {'v1', 'v1_5', 'v2'};
axisStruct = struct('v1', D.v1.Axis_Effect, 'v1_5', D.v1_5.Axis_Effect, 'v2', D.v2.Axis_Effect);
obliStruct = struct('v1', D.v1.Obli_Effect, 'v1_5', D.v1_5.Obli_Effect, 'v2', D.v2.Obli_Effect);

[effect_axis, group_axis, labels_axis] = CrossV_collect_data(axisStruct, compare_groups);
labels_axis = cellfun(@(v) exp_short_label(v), labels_axis, 'UniformOutput', false);

ax = nexttile(tl7, 1);
plot_by_temp_figure(ax, @() CrossV_plot_violin_anova(effect_axis, group_axis, labels_axis, 'Axis Proportion', 0.5));
title(ax, 'Axis Effect');
add_panel_label(ax, 'A');

[effect_obli, group_obli, labels_obli] = CrossV_collect_data(obliStruct, compare_groups);
labels_obli = cellfun(@(v) exp_short_label(v), labels_obli, 'UniformOutput', false);

ax = nexttile(tl7, 2);
plot_by_temp_figure(ax, @() CrossV_plot_violin_anova(effect_obli, group_obli, labels_obli, 'Card. Proportion', 0.5));
title(ax, 'Cardinal Effect');
add_panel_label(ax, 'B');

export_pdf_vector(f7, fullfile(outDir, 'Figure7_Cross_TaskSpace_Composite.pdf'));

% -------------------- Supplemental Figure S1 --------------------
fs1 = figure('Color', 'w', 'Units', 'normalized', 'Position', [0.08 0.18 0.84 0.56]);
tls1 = tiledlayout(fs1, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
for i = 1:numel(vers)
    verc = vers{i};
    vf = strrep(verc, '.', '_');
    dd = D.(vf);
    ax = nexttile(tls1, i);
    mask01 = (dd.start >= 0) & (dd.start < 1000) & (dd.r >= R_min) & (dd.r <= R_max);
    plot_heatmap_on_axes(ax, dd.xpos(mask01), dd.ypos(mask01), heat_binSize, dd.img_width, dd.img_height, ObliColor, CardColor);
    title(ax, sprintf('%s (0-1 s)', exp_short_label(verc)));
    add_panel_label(ax, char('A' + i - 1));
end
enforce_uniform_fontsize(fs1, 11);
export_pdf_vector(fs1, fullfile(outDir, 'Supplemental_Figure_S1_Heatmap_0to1s.pdf'));

% -------------------- Supplemental Figure S2 --------------------
fs2 = figure('Color', 'w', 'Units', 'normalized', 'Position', [0.15 0.20 0.84*(2/3) 0.56*(2/3)]);
tls2 = tiledlayout(fs2, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
for i = 1:numel(vers)
    verc = vers{i};
    vf = strrep(verc, '.', '_');
    dd = D.(vf);

    scan = compute_angle_scan_curves(dd.angles(dd.winSel), dd.sub(dd.winSel), dd.ses(dd.winSel), dd.pairs, dd.Nsubj, angbinSize);

    axTmp = nexttile(tls2, i);
    axPos = get(axTmp, 'Position');
    delete(axTmp);
    ax = polaraxes('Parent', fs2, 'Position', axPos);
    fit_periodic_cosine_phase_copy(scan.sub45, scan.centers45, 45, 0.25, ...
        'plot', true, 'ax', ax, 'histColor', [0.5 0.5 0.5], 'scatterColor', [0 0 0], 'ExpColor', AxisColor, ...
        'statAngle', ternary(strcmpi(vf, 'v1_5'), 22.5, 0), ...
        'rtickTheta', ternary(strcmpi(vf, 'v1_5'), 5, 15), ...
        'rtickHAlign', ternary(strcmpi(vf, 'v1_5'), 'left', 'right'), ...
        'rtickVAlign', 'top');
    title(ax, sprintf('%s: 45° fold phase', exp_short_label(verc)));
    add_panel_label(ax, char('A' + i - 1));
end
enforce_uniform_fontsize(fs2, 11);
export_pdf_vector(fs2, fullfile(outDir, 'Supplemental_Figure_S2_Phase_Distribution.pdf'));

% -------------------- Supplemental Figure S2b (Additional Fold Symmetries: 3, 5, 6, 7fold) --------------------
fs2b = figure('Color', 'w', 'Units', 'normalized', 'Position', [0.08 0.10 0.9 0.75]);
tls2b = tiledlayout(fs2b, 3, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

fold_configs = {
    struct('name', '3fold', 'field_centers', 'centers120', 'field_sub', 'sub120', 'period', 120, 'align_90', true), ...
    struct('name', '5fold', 'field_centers', 'centers72', 'field_sub', 'sub72', 'period', 72, 'align_90', true), ...
    struct('name', '6fold', 'field_centers', 'centers60', 'field_sub', 'sub60', 'period', 60, 'align_90', true), ...
    struct('name', '7fold', 'field_centers', 'centers51_43', 'field_sub', 'sub51_43', 'period', 360/7, 'align_90', true), ...
    struct('name', '4fold', 'field_centers', 'centers90', 'field_sub', 'sub90', 'period', 90, 'align_90', false), ...
    struct('name', '8fold', 'field_centers', 'centers45', 'field_sub', 'sub45', 'period', 45, 'align_90', false)
};

for i = 1:numel(fold_configs)
    cfg = fold_configs{i};
    ax = nexttile(tls2b, i);
    
    verc = 'v2';  % Use v2 data
    vf = strrep(verc, '.', '_');
    dd = D.(vf);
    
    scan = compute_angle_scan_curves(dd.angles(dd.winSel), dd.sub(dd.winSel), dd.ses(dd.winSel), dd.pairs, dd.Nsubj, angbinSize);
    
    centers_fold = scan.(cfg.field_centers);
    sub_fold = scan.(cfg.field_sub);
    
    if cfg.align_90
        plot_angle_curve_full_nosample(ax, sub_fold, centers_fold, cfg.period, CardColor, GapColor, ObliColor, [], true);
        title(ax, sprintf('%s (%.1f° fold, aligned 90°)', cfg.name, cfg.period));
    else
        plot_angle_curve_full_nosample(ax, sub_fold, centers_fold, cfg.period, CardColor, GapColor, ObliColor, [], false);
        title(ax, sprintf('%s (%.1f° fold)', cfg.name, cfg.period));
    end
    add_panel_label(ax, char('A' + i - 1));
end

enforce_uniform_fontsize(fs2b, 11);
export_pdf_vector(fs2b, fullfile(outDir, 'Supplemental_Figure_S2b_Additional_Fold_Symmetries.pdf'));

% -------------------- Supplemental Figure S4 --------------------
fs4 = figure('Color', 'w', 'Units', 'normalized', 'Position', [0.15 0.20 0.84*(2/3) 0.56*(2/3)]);
tls4 = tiledlayout(fs4, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
for i = 1:numel(vers)
    verc = vers{i};
    vf = strrep(verc, '.', '_');
    dd = D.(vf);

    ax = nexttile(tls4, i);
    valid = dd.saccLen ~= 0;
    histogram(ax, dd.saccAngle(valid), -180:5:180, 'FaceColor', [0.35 0.35 0.35], 'EdgeColor', 'none');
    xlabel(ax, 'Saccade Angle (deg)');
    ylabel(ax, 'Count');
    xlim(ax, [-180 180]);
    title(ax, sprintf('%s', exp_short_label(verc)));
    set(ax, 'FontSize', 11, 'LineWidth', 1.0);
    box(ax, 'off');
    add_panel_label(ax, char('A' + i - 1));
end
enforce_uniform_fontsize(fs4, 11);
export_pdf_vector(fs4, fullfile(outDir, 'Supplemental_Figure_S4_Saccade_Direction.pdf'));

% -------------------- Supplemental Figure S5 --------------------
fs5 = figure('Color', 'w', 'Units', 'normalized', 'Position', [0.05 0.06 0.9 0.86]);
tls5 = tiledlayout(fs5, 3, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

for i = 1:numel(vers)
    verc = vers{i};
    vf = strrep(verc, '.', '_');
    dd = D.(vf);

    [sub_time_bin_count, timeCenters_FT] = process_TimeCourse_count(dd.angles, dd.start, dd.dur, dd.sub, dd.ses, dd.pairs, dd.Nsubj, ...
        timeRes_FT, keep_Time, edges_FT, shift_FT, n_bin_FT);
    this_tseri = normalize_by_dim(sub_time_bin_count, 'sum1');

    seriesAxis = squeeze(sum(this_tseri(:,:,1:2:16),3));
    seriesCard = squeeze(sum(this_tseri(:,:,1:4:16),3) ./ max(seriesAxis, eps));

    ax1 = nexttile(tls5, (i-1)*2 + 1);
    cfg1 = struct('doStats', true, 'xlabel', 'Time (ms)', 'ylabel', 'Proportion', ...
        'chanceLevel', 0.5, 'chanceLabel', 'Chance', 'statTail', 'both', 'axesHandle', ax1);
    plot_single_prop(seriesAxis, timeCenters_FT, AxisColor, 'Axis', cfg1);
    xline(ax1, win_left, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.2, 'HandleVisibility', 'off');
    title(ax1, sprintf('%s: Axis Effect', exp_short_label(verc)));
    add_panel_label(ax1, char('A' + (i-1)*2));

    ax2 = nexttile(tls5, (i-1)*2 + 2);
    cfg2 = struct('doStats', true, 'xlabel', 'Time (ms)', 'ylabel', 'Proportion', ...
        'chanceLevel', 0.5, 'chanceLabel', 'Chance', 'statTail', 'right', 'axesHandle', ax2);
    plot_single_prop(seriesCard, timeCenters_FT, CardColor, 'Cardinal', cfg2);
    xline(ax2, win_left, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.2, 'HandleVisibility', 'off');
    title(ax2, sprintf('%s: Cardinal Effect', exp_short_label(verc)));
    add_panel_label(ax2, char('A' + (i-1)*2 + 1));

    if mod(i, 1) == 0
        ylabel(ax2, '');
    end
end

remove_legends_from_figure(fs5);

enforce_uniform_fontsize(fs5, 11);
export_pdf_vector(fs5, fullfile(outDir, 'Supplemental_Figure_S5_WithinTrial_TimeCourse.pdf'));

fprintf('\nComposite figures exported to:\n%s\n', outDir);

end

% ============================== helpers ==============================

function create_main_effect_figure(dm, outName, includeBG, CardColor, ObliColor, GapColor, AxisColor, cmap16_FT, edges_FT, heat_binSize, angbinSize, R_max, outDir, figTag, exportMode, rowHeightScale, show16SectorInd, sectorScatterMarker, sectorScatterAlpha, bPatchAlpha)
    % Build Figure 2-style composites with precise row ratio control.
    fprintf('%s build start\n', figTag);
    % Resize figure: width to half, height to 50% of previous setting.
    % Use a slightly taller canvas when B/G are removed to avoid bottom clipping.
    baseW = round(2300 * 0.5);
    figHScale = 0.50;
    if ~includeBG
        figHScale = 0.56;
    end
    figH = round((900 * 1.70) * figHScale);
    fig = figure('Color', 'w', 'Units', 'pixels', 'Position', [40 40 baseW figH], 'Visible', 'off');

    % Margins and row geometry.
    lm = 0.05; rm = 0.03;
    tm = 0.05;
    if includeBG
        bm = 0.06;
    else
        bm = 0.095;
    end
    rowGap = 0.072;
    hAvail = 1 - tm - bm - 2 * rowGap;
    rowHWeights = rowHeightScale(:)';
    rowHs = hAvail * rowHWeights / sum(rowHWeights);
    h1 = rowHs(1); h2 = rowHs(2); h3 = rowHs(3);
    y3 = bm;
    y2 = y3 + h3 + rowGap;
    y1 = y2 + h2 + rowGap;

    if includeBG
        row1IDs = {'A', 'B', 'C'}; row1R = [2, 1, 3];
        row2IDs = {'D', 'E', 'F'}; row2R = [1, 2, 1];
        row3IDs = {'G', 'H', 'I', 'J'}; row3R = [1.5, 4, 0.7, 1.3];
        labelOrder = {'A','B','C','D','E','F','G','H','I','J'};
    else
        % Remove B and G; keep relative widths and fill each row.
        row1IDs = {'A', 'C'}; row1R = [2, 3];
        row2IDs = {'D', 'E', 'F'}; row2R = [1, 2, 1];
        row3IDs = {'H', 'I', 'J'}; row3R = [4, 0.7, 1.3];
        labelOrder = {'A','C','D','E','F','H','I','J'};
    end

    labelMap = containers.Map();
    for i = 1:numel(labelOrder)
        labelMap(labelOrder{i}) = char('A' + i - 1);
    end

    row1Pos = make_row_positions(row1R, y1, h1, lm, rm, 0.038);
    row2Pos = make_row_positions(row2R, y2, h2, lm, rm, 0.034);
    row3Pos = make_row_positions(row3R, y3, h3, lm, rm, 0.036);

    % Slightly shift B/G left to avoid overlap with the right panel y-labels.
    if includeBG
        row1Pos{2}(1) = row1Pos{2}(1) - 0.022;
        row3Pos{1}(1) = row3Pos{1}(1) - 0.023;
        if strcmpi(figTag, 'Figure 2')
            row1Pos{2}(1) = row1Pos{2}(1) - 0.012; % move panel B content left by ~one letter width
        end
    end

    fprintf('%s scanData start\n', figTag);
    scanData = compute_angle_scan_curves(dm.angles(dm.winSel), dm.sub(dm.winSel), dm.ses(dm.winSel), dm.pairs, dm.Nsubj, angbinSize);
    fprintf('%s scanData done\n', figTag);

    % Draw row 1
    for i = 1:numel(row1IDs)
        id = row1IDs{i};
        fprintf('%s panel %s start\n', figTag, id);
        ax = axes('Parent', fig, 'Position', row1Pos{i});
        draw_panel(id, ax);
        labelDx = 0;
        if strcmpi(figTag, 'Figure 2') && strcmp(id, 'B')
            labelDx = 0.012; % move panel-B label right by ~one letter width
        end
        add_row_aligned_label(fig, ax, labelMap(id), y1 + h1, i == 1, labelDx);
        fprintf('%s panel %s done\n', figTag, id);
    end

    % Draw row 2
    for i = 1:numel(row2IDs)
        id = row2IDs{i};
        fprintf('%s panel %s start\n', figTag, id);
        ax = axes('Parent', fig, 'Position', row2Pos{i});
        draw_panel(id, ax);
        add_row_aligned_label(fig, ax, labelMap(id), y2 + h2, i == 1, 0);

        % Inset pies for the three violin plots in middle row.
        if includeBG
            if strcmp(id, 'D')
                add_pie_inset(ax, [mean(dm.Axis_Effect), mean(dm.Gap_Mean)], [AxisColor; GapColor]);
            elseif strcmp(id, 'E')
                add_pie_inset(ax, [mean(dm.Obli_Effect_inall(:,1)); mean(dm.Obli_Effect_inall(:,2)); mean(dm.Gap_Mean)], [CardColor; ObliColor; GapColor]);
            elseif strcmp(id, 'F')
                add_pie_inset(ax, [mean(dm.Obli_Effect(:,1)); mean(dm.Obli_Effect(:,2))], [CardColor; ObliColor]);
            end
        end
        if any(strcmp(id, {'D','E','F'}))
            shift_dprime_annotation(ax, -0.12);
        end
        if any(strcmp(id, {'D','F'}))
            % Widen single-group violin by reducing x-span while preserving left margin.
            xlim(ax, [0.35, 1.65]);
            set_dprime_right_margin(ax, 0.90);
        end
        fprintf('%s panel %s done\n', figTag, id);
    end

    % Draw row 3
    for i = 1:numel(row3IDs)
        id = row3IDs{i};
        fprintf('%s panel %s start\n', figTag, id);
        ax = axes('Parent', fig, 'Position', row3Pos{i});
        draw_panel(id, ax);
        add_row_aligned_label(fig, ax, labelMap(id), y3 + h3, i == 1, 0);
        fprintf('%s panel %s done\n', figTag, id);
    end

    fprintf('%s export start: %s\n', figTag, outName);
    enforce_uniform_fontsize(fig, 11);
    export_pdf_by_mode(fig, fullfile(outDir, outName), exportMode);
    fprintf('%s export done: %s\n', figTag, outName);

    function draw_panel(id, ax)
        switch id
            case 'A'
                plot_heatmap_on_axes(ax, dm.xpos(dm.winSel), dm.ypos(dm.winSel), heat_binSize, dm.img_width, dm.img_height, ObliColor, CardColor);
                title(ax, 'Heatmap');
            case 'B'
                draw_sector_scatter_panel(ax, dm, R_max, cmap16_FT, edges_FT, true, sectorScatterMarker * 0.72, sectorScatterAlpha * 0.55, true, bPatchAlpha);
                title(ax, 'Sectors');
            case 'C'
                plot_16sector_prop_from_chance(ax, dm.mean_sectors_prop, cmap16_FT, dm.labels_16bin, 1/16, show16SectorInd);
                title(ax, '16-Sector Proportion');
            case 'D'
                plot_by_temp_figure(ax, @() plot_violin_prop(dm.Axis_Effect, AxisColor, {'Axis'}, ...
                    'ylabel', 'Proportion', 'xlabel', '', 'title', '', 'showInd', true, ...
                    'showIndNum', false, 'showIndLink', false, 'chanceLevel', 0.5));
                title(ax, 'Axis Effect');
            case 'E'
                plot_by_temp_figure(ax, @() plot_violin_prop(dm.Obli_Effect_inall, cmap16_FT([1,3],:), {'Card', 'Obli'}, ...
                    'ylabel', 'Proportion', 'xlabel', '', 'title', '', 'showInd', true, ...
                    'showIndNum', false, 'showIndLink', true, 'chanceLevel', mean(dm.Gap_Mean)/2, 'chanceLabel', 'Gap/2'));
                if strcmpi(figTag, 'Figure 3')
                    % Figure 3 panel D: make marginal p-value annotation muted and lower it to avoid overlap.
                    style_marginal_violin_sig(ax, 0.058, 0.012, -0.010);
                end
                title(ax, 'Card vs Obli');
            case 'F'
                plot_by_temp_figure(ax, @() plot_violin_prop(dm.Obli_Effect(:,1), cmap16_FT(1,:), {'Card'}, ...
                    'ylabel', 'Proportion', 'xlabel', '', 'title', '', 'showInd', true, ...
                    'showIndNum', false, 'showIndLink', false, 'chanceLevel', 0.5, 'chanceLabel', 'Chance'));
                title(ax, 'Card in Axis');
            case 'G'
                draw_scan_schematic(ax, dm.angles(dm.winSel), dm.r(dm.winSel), R_max, CardColor, ObliColor, GapColor, sectorScatterMarker, sectorScatterAlpha);
                title(ax, 'Scan Schematic');
            case 'H'
                plot_angle_curve_full_nosample(ax, scanData.sub360, scanData.centers360, 360, CardColor, GapColor, ObliColor, [0.02, 0.04]);
                title(ax, '360° Scan');
            case 'I'
                plot_angle_curve_full_nosample(ax, scanData.sub45, scanData.centers45, 45, CardColor, GapColor, ObliColor, [0.22, 0.30]);
                ylabel(ax, '');
                title(ax, '45° Fold');
            case 'J'
                plot_angle_curve_full_nosample(ax, scanData.sub90, scanData.centers90, 90, CardColor, GapColor, ObliColor, [0.105, 0.162]);
                ylabel(ax, '');
                title(ax, '90° Fold');
        end
    end
end

function rowPos = make_row_positions(ratios, y, h, lm, rm, gap)
    n = numel(ratios);
    usableW = 1 - lm - rm - (n - 1) * gap;
    widths = usableW * ratios / sum(ratios);
    rowPos = cell(1, n);
    x = lm;
    for i = 1:n
        rowPos{i} = [x, y, widths(i), h];
        x = x + widths(i) + gap;
    end
end

function add_pie_inset(parentAx, values, colors)
    fig = ancestor(parentAx, 'figure');
    pos = get(parentAx, 'Position');

    % Keep identical inset size across D/E/F regardless of parent width.
    w = 0.085 * 0.816;
    h = 0.105 * 0.816;
    x = pos(1) + 0.008;
    d = min(w, h);
    y = pos(2) + pos(4) - h - 0.01 + 0.125 * d;

    axIn = axes('Parent', fig, 'Position', [x y w h], 'Color', 'none');
    plot_pie(values(:), colors, 'ax', axIn, 'show_text', false);
    set(axIn, 'Color', 'none');
end

function add_row_aligned_label(fig, ax, labelText, rowTopY, isRowFirst, xShift)
    if nargin < 5
        isRowFirst = false;
    end
    if nargin < 6
        xShift = 0;
    end
    pos = get(ax, 'Position');
    if isRowFirst
        x = 0.006;
    else
        x = max(0.005, pos(1) - 0.018);
    end
    x = x + xShift;
    y = min(0.995, rowTopY + 0.006);
    annotation(fig, 'textbox', [x y 0.02 0.02], 'String', labelText, ...
        'LineStyle', 'none', 'FontWeight', 'bold', 'FontSize', 11, ...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');
end

function draw_sector_scatter_panel(ax, dm, R_max, cmap16, edges_FT, allowDownsample, markerSize, markerAlpha, showAngleLabels, bPatchAlpha)
    xpos = dm.xpos(dm.winSel);
    ypos = dm.ypos(dm.winSel);
    if allowDownsample
        nMax = 65000;
        if numel(xpos) > nMax
            rp = randperm(numel(xpos), nMax);
            xpos = xpos(rp);
            ypos = ypos(rp);
        end
    end
    % Draw sector-colored scatter with lighter background, matching PLOT_SECTORS style.
    center_pix = [dm.img_width, dm.img_height] / 2;
    ut = UT(53.2, dm.img_width, 57, false);
    RpatchPix = ut.deg2pix(R_max);
    theta_patch_res = 100;
    light_cmap = min(1, 0.72 + 0.28 * cmap16);
    ringInner = 0.14 * RpatchPix;

    hold(ax, 'on');
    for si = 1:16
        a1 = edges_FT(si);
        a2 = edges_FT(si + 1);
        if a2 < a1, a2 = a2 + 360; end
        aa = linspace(a1, a2, theta_patch_res);
        xo = center_pix(1) + RpatchPix * cosd(aa);
        yo = center_pix(2) - RpatchPix * sind(aa);
        xi = center_pix(1) + ringInner * cosd(fliplr(aa));
        yi = center_pix(2) - ringInner * sind(fliplr(aa));
        patch(ax, [xo, xi], [yo, yi], light_cmap(si,:), 'EdgeColor', [1 1 1], 'LineWidth', 0.45, 'FaceAlpha', bPatchAlpha);
    end

    % Keep a bright center so annular segmentation remains visible behind points.
    th0 = linspace(0, 2*pi, 180);
    patch(ax, center_pix(1) + ringInner * cos(th0), center_pix(2) - ringInner * sin(th0), [1 1 1], ...
        'EdgeColor', 'none', 'FaceAlpha', min(0.75, bPatchAlpha));

    angA = mod(atan2d(center_pix(2) - ypos, xpos - center_pix(1)), 360);
    sectorIdx = zeros(numel(angA), 1, 'uint8');
    for si = 1:16
        l = mod(edges_FT(si), 360);
        r = mod(edges_FT(si + 1), 360);
        if l < r
            mask = angA >= l & angA < r;
        else
            mask = angA >= l | angA < r;
        end
        sectorIdx(mask) = si;
    end
    ptColors = cmap16(double(sectorIdx), :);
    scatter(ax, xpos, ypos, markerSize, ptColors, 'filled', 'MarkerFaceAlpha', markerAlpha, 'MarkerEdgeAlpha', 0.10);

    if showAngleLabels
        for a = 0:22.5:337.5
            tx = center_pix(1) + 1.24 * RpatchPix * cosd(a);
            ty = center_pix(2) - 1.24 * RpatchPix * sind(a);
            text(ax, tx, ty, sprintf('%.1f°', a), 'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'middle', 'FontSize', 8.5, 'FontWeight', 'bold', 'Color', [0.16 0.16 0.16]);
        end
    end
    hold(ax, 'off');

    axis(ax, 'equal');
    xlim(ax, center_pix(1) + [-1.30, 1.30] * RpatchPix);
    ylim(ax, center_pix(2) + [-1.30, 1.30] * RpatchPix);
    set(ax, 'YDir', 'reverse', 'XTick', [], 'YTick', [], 'Visible', 'off');
end

function plot_angle_curve_light(ax, data, centers, foldPeriod, CardColor, GapColor, ObliColor, yRange)
    % Lightweight replacement of plot_angle_curve for stable vector export.
    m = mean(data, 1, 'omitnan');
    se = std(data, 0, 1, 'omitnan') ./ sqrt(max(1, size(data,1)));

    % Decimate to keep vector PDF manageable.
    if foldPeriod == 360
        step = 20;
    elseif foldPeriod == 90
        step = 3;
    else
        step = 2;
    end
    idx = 1:step:numel(centers);
    x = centers(idx);
    mm = m(idx);
    ss = se(idx);

    cla(ax);
    hold(ax, 'on');
    fill(ax, [x, fliplr(x)], [mm + ss, fliplr(mm - ss)], [0.7 0.7 0.7], ...
        'EdgeColor', 'none', 'FaceAlpha', 0.35);
    plot(ax, x, mm, 'Color', [0.3 0.3 0.3], 'LineWidth', 1.6);

    if nargin >= 8 && ~isempty(yRange)
        ylim(ax, yRange);
    end
    yl = ylim(ax);

    if foldPeriod == 360
        for a = 0:45:315
            if mod(a, 90) == 0
                lc = CardColor;
            else
                lc = ObliColor;
            end
            xline(ax, a, '--', 'Color', lc, 'LineWidth', 0.8, 'HandleVisibility', 'off');
            xline(ax, a + 22.5, '--', 'Color', GapColor, 'LineWidth', 0.8, 'HandleVisibility', 'off');
        end
        xticks(ax, 0:45:315);
        xticklabels(ax, arrayfun(@(v) sprintf('%.1f°', v), 0:45:315, 'UniformOutput', false));
        xlabel(ax, 'Angle Mod 360°');
    elseif foldPeriod == 45
        for a = [0, 22.5]
            xline(ax, a, '--', 'Color', ternary(a==0, CardColor, GapColor), 'LineWidth', 0.9, 'HandleVisibility', 'off');
        end
        xticks(ax, [0, 22.5]);
        xticklabels(ax, {'0.0°','22.5°'});
        xlabel(ax, 'Angle Mod 45°');
    elseif foldPeriod == 90
        for a = [0, 22.5, 45, 67.5]
            if a == 0
                lc = CardColor;
            elseif a == 45
                lc = ObliColor;
            else
                lc = GapColor;
            end
            xline(ax, a, '--', 'Color', lc, 'LineWidth', 0.9, 'HandleVisibility', 'off');
        end
        xticks(ax, [0, 22.5, 45, 67.5]);
        xticklabels(ax, {'0.0°','22.5°','45.0°','67.5°'});
        xlabel(ax, 'Angle Mod 90°');
    end

    xlim(ax, [min(x), max(x)]);
    ylim(ax, yl);
    ylabel(ax, 'Count Proportion');
    set(ax, 'FontSize', 10.5, 'LineWidth', 1.0);
    box(ax, 'off');
    hold(ax, 'off');
end

function d = preprocess_from_fixTable(T)
    if islogical(T.dropFix)
        keepMask = ~T.dropFix;
    else
        keepMask = T.dropFix == 0;
    end
    vis = T(keepMask, :);
    d.visTable = vis;

    d.img_width = 1920;
    d.img_height = 1080;

    d.start = vis.startT;
    d.dur = vis.dur;
    d.angles = vis.theta;
    d.xpos = vis.xpos;
    d.ypos = vis.ypos;
    d.sub = vis.subID;
    d.ses = vis.sessID;
    d.r = vis.r;
    d.saccLen = vis.saccLen;
    d.saccAngle = vis.saccAngle;

    d.pairs = unique([d.sub, d.ses], 'rows', 'stable');
    d.Nsubj = size(d.pairs, 1);

    win_left = 1000;
    win_right = 4000;
    R_min = 1;
    R_max = 7.5;

    d.winSel = ((win_left <= d.start & d.start <= win_right) | ...
               (win_left <= d.start + d.dur & d.start + d.dur <= win_right)) & ...
               (d.r >= R_min & d.r <= R_max);

    n_bin_FT = 16;
    shift_FT = 360 / n_bin_FT / 2;
    edges_FT = linspace(0, 360, n_bin_FT + 1) - shift_FT;

    bin_counts = zeros(d.Nsubj, n_bin_FT);
    for si = 1:d.Nsubj
        m = d.winSel & d.sub == d.pairs(si,1) & d.ses == d.pairs(si,2);
        bin_counts(si,:) = histcounts(mod(d.angles(m) + shift_FT, 360) - shift_FT, edges_FT);
    end

    d.mean_sectors_prop = normalize_by_dim(bin_counts, 'sum1');
    d.Axis_Effect = sum(bin_counts(:,1:2:16), 2) ./ max(sum(bin_counts,2), eps);
    d.Obli_Effect = [sum(bin_counts(:,1:4:16),2), sum(bin_counts(:,3:4:16),2)] ./ max(sum(bin_counts(:,1:2:16),2), eps);
    d.Obli_Effect_inall = [sum(bin_counts(:,1:4:16),2), sum(bin_counts(:,3:4:16),2)] ./ max(sum(bin_counts,2), eps);
    d.Gap_Mean = sum(bin_counts(:,2:2:16),2) ./ max(sum(bin_counts,2), eps);

    d.labels_16bin = arrayfun(@(x) sprintf('%.0f°', x), (edges_FT(1:end-1)+edges_FT(2:end))/2, 'uni', 0);
end

function cfg_win = build_trial_cfg(d, edges_FT, shift_FT, cmap16_FT, win_left, win_right, AxisColor, learn_stage_n, last_trial)
    total_trial = last_trial - learn_stage_n;

    % Build tail mask from each subject-session max trial.
    pairN = size(d.pairs, 1);
    maxNtrial = nan(pairN, 1);
    for i = 1:pairN
        m = d.visTable.subID == d.pairs(i,1) & d.visTable.sessID == d.pairs(i,2);
        if any(m)
            maxNtrial(i) = max(d.visTable.Ntrial(m));
        end
    end
    tailMask = false(pairN, total_trial);
    for i = 1:pairN
        if ~isnan(maxNtrial(i)) && maxNtrial(i) < total_trial
            tailMask(i, maxNtrial(i)+1:end) = true;
        end
    end

    cfg_win = struct();
    cfg_win.edges_FT = edges_FT;
    cfg_win.shift_FT = shift_FT;
    cfg_win.cmap16_FT = cmap16_FT;
    cfg_win.win_left = win_left;
    cfg_win.win_right = win_right;
    cfg_win.step_trials = 1;
    cfg_win.total_trial = total_trial;
    cfg_win.digPlace = floor(log10(max(d.visTable.Ntrial)));
    cfg_win.mode = 'continuous';
    cfg_win.doStats = true;
    cfg_win.learn_stage_n = learn_stage_n;
    cfg_win.normMode_FT = 'zScore';
    cfg_win.xlabel = 'Trial number';
    cfg_win.ylabel = 'Proportion';
    cfg_win.doSmooth = false;
    cfg_win.ver = 'v1';
    cfg_win.AxisColor = AxisColor;
    cfg_win.Card_in_Axis = true;
    cfg_win.cutoff = 300;
    cfg_win.tailMask = tailMask;
end

function out = ternary(cond, a, b)
    if cond
        out = a;
    else
        out = b;
    end
end

function lbl = exp_short_label(verc)
    v = lower(string(verc));
    if v == "v1"
        lbl = 'EXP1a';
    elseif v == "v1.5"
        lbl = 'EXP1b';
    else
        lbl = 'EXP2';
    end
end

function scan = compute_angle_scan_curves(angles, sub, ses, pairs, nsbj, angbinSize)
    startAngle = -22.5/2;
    [centers360, ~, ~, sub360] = analyze_angle_curve(angles, sub, ses, pairs, nsbj, angbinSize, 360, 'sum1', startAngle);

    % Original fold configs aligned to 0° (4fold, 8fold)
    fold_helper_0 = @(period) fold_angle_data(centers360, sub360, nsbj, startAngle, period);
    
    % New fold configs aligned to 90° with start_angle info for plotting
    fold_helper_90 = @(period) fold_angle_data(centers360, sub360, nsbj, 90 - period/2, period);
    
    % 3fold (120°) - aligned to 90°
    [centers120, sub120] = fold_helper_90(120);
    
    % 4fold (90°) - aligned to 0°
    [centers90, sub90] = fold_helper_0(90);
    
    % 5fold (72°) - aligned to 90°
    [centers72, sub72] = fold_helper_90(72);
    
    % 6fold (60°) - aligned to 90°
    [centers60, sub60] = fold_helper_90(60);
    
    % 7fold (51.43°) - aligned to 90°
    [centers51_43, sub51_43] = fold_helper_90(360/7);
    
    % 8fold (45°) - aligned to 0°
    [centers45, sub45] = fold_helper_0(45);

    scan = struct('centers360', centers360, 'sub360', sub360, ...
                  'centers45', centers45, 'sub45', sub45, ...
                  'centers90', centers90, 'sub90', sub90, ...
                  'centers120', centers120, 'sub120', sub120, ...
                  'centers72', centers72, 'sub72', sub72, ...
                  'centers60', centers60, 'sub60', sub60, ...
                  'centers51_43', centers51_43, 'sub51_43', sub51_43);
end

function [centers_folded, sub_folded] = fold_angle_data(centers360, sub360, nsbj, startAngle, foldPeriod)
    mod_angles = mod(centers360 - startAngle, foldPeriod);
    [uniq_mod, ~, ic] = unique(round(mod_angles, 8));
    sub_folded = zeros(nsbj, numel(uniq_mod));
    for k = 1:numel(uniq_mod)
        idx = (ic == k);
        sub_folded(:,k) = sum(sub360(:,idx), 2, 'omitnan');
    end
    centers_folded = uniq_mod + startAngle;
end

function plot_heatmap_on_axes(ax, xpos, ypos, heat_binSize, img_width, img_height, ObliColor, CardColor)
    ex = [1, img_width];
    ey = [1, img_height];

    nx = max(10, round((ex(2)-ex(1))/heat_binSize));
    ny = max(10, round((ey(2)-ey(1))/heat_binSize));
    edges_x = linspace(ex(1), ex(2), nx+1);
    edges_y = linspace(ey(1), ey(2), ny+1);

    [count2D,~,~] = histcounts2(xpos, ypos, edges_x, edges_y);
    density = count2D' / max(1, sum(count2D(:)));

    imagesc(ax, edges_x, edges_y, density);
    axis(ax, 'image');
    set(ax, 'YDir', 'reverse');
    nMap = 256;
    lightCard = min(1, CardColor + [0.22, 0.28, 0.22]);
    mapBG = [linspace(ObliColor(1), lightCard(1), nMap)', ...
             linspace(ObliColor(2), lightCard(2), nMap)', ...
             linspace(ObliColor(3), lightCard(3), nMap)'];
    colormap(ax, mapBG);
    cmax = prctile(density(:), 100);
    if isempty(cmax) || cmax == 0
        clim(ax, 'auto');
    else
        clim(ax, [min(density(:)), cmax]);
    end
    colorbar(ax);
    xlabel(ax, 'X Position (pixel)');
    ylabel(ax, 'Y Position (pixel)');
end

function draw_sector_partition_scheme(ax, video_width, video_height, R_max, cmap16, sector_edges)
    center_pix = [video_width, video_height] / 2;
    ut = UT(53.2, video_width, 57, false);
    RpatchPix = ut.deg2pix(R_max);

    theta_patch_res = 120;
    hold(ax, 'on');

    for si = 1:16
        a1 = sector_edges(si);
        a2 = sector_edges(si + 1);
        if a2 < a1
            a2 = a2 + 360;
        end
        aa = linspace(a1, a2, theta_patch_res);
        xp = center_pix(1) + RpatchPix * cosd(aa);
        yp = center_pix(2) - RpatchPix * sind(aa);
        P = [center_pix; [xp(:), yp(:)]];
        patch(ax, 'XData', P(:,1), 'YData', P(:,2), 'FaceColor', cmap16(si,:), ...
            'EdgeColor', [1 1 1], 'LineWidth', 0.8, 'FaceAlpha', 0.9);
    end

    for si = 1:16
        midA = mod((sector_edges(si)+sector_edges(si+1))/2, 360);
        x = center_pix(1) + 1.08*RpatchPix*cosd(midA);
        y = center_pix(2) - 1.08*RpatchPix*sind(midA);
        text(ax, x, y, sprintf('%.1f°', midA), 'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'middle', 'FontSize', 8, 'Color', [0.1 0.1 0.1]);
    end

    axis(ax, 'equal');
    xlim(ax, center_pix(1) + [-1.35, 1.35]*RpatchPix);
    ylim(ax, center_pix(2) + [-1.35, 1.35]*RpatchPix);
    set(ax, 'YDir', 'reverse', 'XTick', [], 'YTick', [], 'Visible', 'off');

    text(ax, center_pix(1), center_pix(2) + 1.23*RpatchPix, 'Green: Cardinal; Blue: Oblique; Gray: Gap', ...
        'HorizontalAlignment', 'center', 'FontSize', 9, 'Color', [0.1 0.1 0.1]);
    hold(ax, 'off');
end

function draw_scan_schematic(ax, anglesDeg, radiiDeg, R_max, CardColor, ObliColor, GapColor, markerSize, markerAlpha)
    cla(ax);
    hold(ax, 'on');

    % Reduce point count to avoid heavy export objects while preserving distribution.
    nMaxPts = 14000;
    nAll = numel(anglesDeg);
    if nAll > nMaxPts
        rp = randperm(nAll, nMaxPts);
        anglesDeg = anglesDeg(rp);
        radiiDeg = radiiDeg(rp);
    end

    rr = max(0, min(radiiDeg ./ max(R_max, eps), 1));
    xx = rr .* cosd(anglesDeg);
    yy = rr .* sind(anglesDeg);
    scatter(ax, xx, yy, markerSize, [0.72 0.72 0.72], 'filled', ...
        'MarkerFaceAlpha', markerAlpha, 'MarkerEdgeAlpha', 0);

    % Outer circle.
    th = linspace(0, 2*pi, 720);
    plot(ax, cos(th), sin(th), '-', 'Color', [0.72 0.72 0.72], 'LineWidth', 2.0);

    % Scan window example wedge around 0°.
    wHalf = 11.25/2;
    tw = deg2rad(linspace(-wHalf, wHalf, 60));
    patch(ax, [0, cos(tw), 0], [0, sin(tw), 0], [0.08 0.08 0.08], ...
        'FaceAlpha', 0.42, 'EdgeColor', 'none');

    % Sector guides every 22.5°.
    cardGuide = min(1, 0.35 + 0.65 * CardColor);
    obliGuide = min(1, 0.35 + 0.65 * ObliColor);
    gapGuide  = min(1, 0.35 + 0.65 * GapColor);
    for a = 0:22.5:337.5
        if any(abs(a - [0, 90, 180, 270]) < 1e-9)
            lc = cardGuide;
        elseif any(abs(a - [45, 135, 225, 315]) < 1e-9)
            lc = obliGuide;
        else
            lc = gapGuide;
        end
        plot(ax, [0 cosd(a)], [0 sind(a)], '--', 'Color', lc, 'LineWidth', 1.8);
    end

    % Angle labels around circle.
    for a = 0:22.5:337.5
        tx = 1.15 * cosd(a);
        ty = 1.15 * sind(a);
        text(ax, tx, ty, sprintf('%.1f°', a), 'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'middle', 'FontSize', 8.5, 'FontWeight', 'bold', 'Color', [0.1 0.1 0.1]);
    end

    % Curved arrow outside right side.
    tArrow = linspace(-25, 25, 40);
    xa = 1.23 * cosd(tArrow);
    ya = 1.23 * sind(tArrow);
    plot(ax, xa, ya, '-', 'Color', [0.35 0.35 0.35], 'LineWidth', 3);
    quiver(ax, xa(end-1), ya(end-1), xa(end)-xa(end-1), ya(end)-ya(end-1), 0, ...
        'Color', [0.35 0.35 0.35], 'LineWidth', 2.2, 'MaxHeadSize', 6);

    axis(ax, 'equal');
    xlim(ax, [-1.35 1.35]);
    ylim(ax, [-1.25 1.25]);
    axis(ax, 'off');
    hold(ax, 'off');
end

function plot_16sector_prop_from_chance(ax, data, colors, labels, baseline, showIndividual)
    [nSub, nB] = size(data);
    m = mean(data, 1, 'omitnan');
    se = std(data, 0, 1, 'omitnan') ./ sqrt(max(1, nSub));

    cla(ax);
    hold(ax, 'on');
    hb = bar(ax, 1:nB, m, 'FaceColor', 'flat', 'BaseValue', baseline, 'LineWidth', 1.1);
    for ii = 1:nB
        hb.CData(ii,:) = colors(min(ii,size(colors,1)),:);
    end
    hb.FaceAlpha = 0.35;
    hb.EdgeColor = 'flat';
    errorbar(ax, 1:nB, m, se, 'LineStyle', 'none', 'Color', [0.2 0.2 0.2], 'LineWidth', 1.0);

    % Individual subject dots (toggle).
    if showIndividual
        for ib = 1:nB
            xj = ib + 0.06 * randn(nSub,1);
            scatter(ax, xj, data(:,ib), 16, repmat(colors(ib,:), nSub, 1), ...
                'filled', 'MarkerFaceAlpha', 0.6, 'MarkerEdgeAlpha', 0.35);
        end
    end

    % Baseline and significance marks.
    yline(ax, baseline, '--', 'Color', [0.45 0.45 0.45], 'LineWidth', 1.0, 'Label', 'Chance', ...
        'LabelHorizontalAlignment', 'left');
    pvals = nan(1, nB);
    for i = 1:nB
        try
            [~, pvals(i)] = ttest(data(:,i), baseline);
        catch
            pvals(i) = NaN;
        end
    end
    [~, ~, ~, adj_p] = fdr_bh(pvals, 0.05, 'pdep', 'yes');
    bump = 0.04 * max(range(ylim(ax)), eps);
    for i = 1:nB
        s = sig_symbol(adj_p(i));
        if isempty(s)
            s = 'n.s.';
        end
        text(ax, i, m(i) + sign(m(i)-baseline) * (se(i) + bump), s, ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
            'FontSize', 10.5, 'FontWeight', 'bold');
    end

    xticks(ax, 1:2:nB);
    xticklabels(ax, labels(1:2:nB));
    ylabel(ax, 'Proportion');
    box(ax, 'off');
    set(ax, 'FontSize', 11, 'LineWidth', 1.0);
    hold(ax, 'off');
end

function add_panel_label(ax, labelText)
    text(ax, -0.12, 1.08, labelText, 'Units', 'normalized', ...
        'FontWeight', 'bold', 'FontSize', 13, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');
end

function export_pdf_vector(figHandle, outPdf)
    set(figHandle, 'Renderer', 'painters');
    exportgraphics(figHandle, outPdf, 'ContentType', 'vector');
    exportgraphics(figHandle, replace(outPdf, '.pdf', '.png'), 'Resolution', 600);
    savefig(figHandle, replace(outPdf, '.pdf', '.fig'));
    close(figHandle);
end

function export_pdf_by_mode(figHandle, outPdf, mode)
    if nargin < 3 || isempty(mode)
        mode = 'vector';
    end
    mode = lower(string(mode));
    if mode == "image"
        set(figHandle, 'Renderer', 'opengl');
        exportgraphics(figHandle, outPdf, 'ContentType', 'image', 'Resolution', 600);
        exportgraphics(figHandle, replace(outPdf, '.pdf', '.png'), 'Resolution', 600);
        savefig(figHandle, replace(outPdf, '.pdf', '.fig'));
        close(figHandle);
    else
        export_pdf_vector(figHandle, outPdf);
    end
end

function plot_by_temp_figure(targetAx, drawFcn)
    figsBefore = findall(0, 'Type', 'figure');
    drawFcn();
    figsAfter = findall(0, 'Type', 'figure');

    newFigs = setdiff(figsAfter, figsBefore);
    if isempty(newFigs)
        return;
    end

    srcFig = newFigs(1);
    srcAxes = findall(srcFig, 'Type', 'axes');
    srcAxes = srcAxes(arrayfun(@(a) ~strcmpi(get(a, 'Tag'), 'legend'), srcAxes));
    if isempty(srcAxes)
        close(srcFig);
        return;
    end

    % Use the first non-legend axis as source.
    srcAx = srcAxes(end);

    cla(targetAx);
    copyobj(allchild(srcAx), targetAx);

    set(targetAx, 'XLim', get(srcAx, 'XLim'));
    set(targetAx, 'YLim', get(srcAx, 'YLim'));
    set(targetAx, 'XScale', get(srcAx, 'XScale'));
    set(targetAx, 'YScale', get(srcAx, 'YScale'));
    set(targetAx, 'XTick', get(srcAx, 'XTick'));
    set(targetAx, 'YTick', get(srcAx, 'YTick'));
    set(targetAx, 'XTickLabel', get(srcAx, 'XTickLabel'));
    set(targetAx, 'YTickLabel', get(srcAx, 'YTickLabel'));
    set(targetAx, 'FontSize', get(srcAx, 'FontSize'));
    set(targetAx, 'LineWidth', get(srcAx, 'LineWidth'));
    set(targetAx, 'Box', get(srcAx, 'Box'));

    xlabel(targetAx, get(get(srcAx, 'XLabel'), 'String'));
    ylabel(targetAx, get(get(srcAx, 'YLabel'), 'String'));
    title(targetAx, get(get(srcAx, 'Title'), 'String'));

    close(srcFig);
end

function [fitResult, stats] = fit_periodic_cosine_phase_copy(subCountsFold, centersFold, period, chance, varargin)
    p = inputParser;
    addParameter(p, 'plot', true, @(x)islogical(x)||isnumeric(x));
    addParameter(p, 'histColor', [0.5 0.5 0.5]);
    addParameter(p, 'scatterColor', [1 0 0]);
    addParameter(p, 'ExpColor', [1 0 0]);
    addParameter(p, 'ax', []);
    addParameter(p, 'statAngle', 0);
    addParameter(p, 'rtickTheta', 120);
    addParameter(p, 'rtickHAlign', 'right');
    addParameter(p, 'rtickVAlign', 'top');
    parse(p, varargin{:});

    doPlot = logical(p.Results.plot);
    histColor = p.Results.histColor;
    scatterColor = p.Results.scatterColor;
    ExpColor = p.Results.ExpColor;
    ax = p.Results.ax;
    statAngle = p.Results.statAngle;
    rtickTheta = p.Results.rtickTheta;
    rtickHAlign = p.Results.rtickHAlign;
    rtickVAlign = p.Results.rtickVAlign;

    [Nsubj, ~] = size(subCountsFold);
    omega = 2*pi / period;
    theta_deg = centersFold(:)';
    cos_term = cos(omega * theta_deg);
    sin_term = sin(omega * theta_deg);

    phases = nan(Nsubj,1);
    amps = nan(Nsubj,1);
    beta12 = nan(Nsubj,2);

    for si = 1:Nsubj
        y = subCountsFold(si,:);
        valid = ~isnan(y);
        X = [cos_term(valid)', sin_term(valid)'];
        b = X \ (y(valid)' - chance);
        beta1 = b(1);
        beta2 = b(2);
        R = hypot(beta1, beta2);
        phi = atan2(beta2, beta1);
        phi = mod(phi + pi, 2*pi) - pi;
        phases(si) = phi;
        amps(si) = R;
        beta12(si,:) = [beta1, beta2];
    end

    theta0_deg = phases / omega;
    theta0_deg_mod = mod(theta0_deg, period);

    valid_idx = ~isnan(phases);
    angles_test = phases(valid_idx);
    n_angles = numel(angles_test);
    phi0 = omega * mod(statAngle, period);

    if n_angles >= 2
        if exist('circ_vtest','file') == 2
            [vtest_p, vtest_stat] = circ_vtest(angles_test, phi0);
            vtest_method = 'circ_vtest';
        elseif exist('circ_rtest','file') == 2
            [vtest_p, vtest_stat] = circ_rtest(angles_test - phi0);
            vtest_method = 'circ_rtest (fallback)';
        else
            R = abs(sum(exp(1i*(angles_test - phi0))));
            Rbar_tmp = R / n_angles;
            Z = n_angles * Rbar_tmp^2;
            p_approx = exp(-Z) * (1 + (2*Z - Z^2)/(4*n_angles) - (24*Z - 132*Z^2 + 76*Z^3 - 9*Z^4)/(288*n_angles^2));
            vtest_p = min(max(p_approx,0),1);
            vtest_stat = Z;
            vtest_method = 'rayleigh_approx (no CircStat)';
        end
        mean_dir = angle(mean(exp(1i*angles_test)));
        Rbar = abs(mean(exp(1i*angles_test)));
        mean_theta0_deg = mod((mean_dir / omega), period);
    else
        vtest_p = NaN;
        vtest_stat = NaN;
        vtest_method = 'insufficient';
        mean_dir = NaN;
        Rbar = NaN;
        mean_theta0_deg = NaN;
    end

    fitResult.phase_rad = phases;
    fitResult.phase_deg = rad2deg(phases);
    fitResult.amp = amps;
    fitResult.beta = beta12;
    fitResult.theta0_deg = theta0_deg;
    fitResult.theta0_deg_mod = theta0_deg_mod;
    fitResult.vtest_method = vtest_method;
    fitResult.vtest_p = vtest_p;
    fitResult.vtest_stat = vtest_stat;
    fitResult.mean_phi = mean_dir;
    fitResult.mean_theta0_deg = mean_theta0_deg;
    fitResult.Rbar = Rbar;
    fitResult.valid_idx = valid_idx;

    stats.n_valid = sum(valid_idx);
    stats.mean_amp = mean(amps, 'omitnan');
    stats.std_amp = std(amps, 0, 'omitnan');
    stats.mean_phi = mean_dir;
    stats.mean_theta0_deg = mean_theta0_deg;
    stats.Rbar = Rbar;
    stats.vtest_p = vtest_p;
    stats.vtest_stat = vtest_stat;
    stats.vtest_method = vtest_method;

    if doPlot
        if isempty(ax) || ~isvalid(ax)
            f = figure('Color', 'none', 'Position', [120,80,270,240]); %#ok<NASGU>
            ax = polaraxes;
        end

        theta_plot = theta0_deg_mod / period * 2*pi;
        nBins = 16;
        binWidth = 2*pi / nBins;
        binEdges = linspace(0, 2*pi, nBins+1) - binWidth/2;

        polarhistogram(ax, theta_plot(~isnan(theta_plot)), binEdges, ...
            'Normalization', 'probability', 'FaceAlpha', 0.5, 'FaceColor', histColor, 'EdgeColor', histColor);
        hold(ax, 'on');

        r_plot = amps;
        rmin = 1e-16;
        r_plot_log = log10(max(r_plot, rmin));
        r_plot_log_scaled = (r_plot_log + 3) / 10;
        polarplot(ax, theta_plot, r_plot_log_scaled, '.', 'MarkerSize', 12, 'Color', scatterColor);

        theta_stat = mod(statAngle, period) / period * 2*pi;
        rmax = 0.32;
        polarplot(ax, [theta_stat theta_stat], [0 rmax], '--', 'LineWidth', 2, 'Color', ExpColor);

        ax.RLim = [0 0.3];
        log_ticks = log10([0.001 0.01 0.1 1]);
        log_ticks_scaled = (log_ticks + 3) / 10;
        log_vals = {'0.001(0)', '0.01(0.1)', '0.1(0.2)', '1(0.3)'};
        ax.RTick = log_ticks_scaled;
        ax.RTickLabel = {};

        theta = deg2rad(rtickTheta * (360/period));
        for ii = 1:length(log_vals)
            text(ax, theta, log_ticks_scaled(ii), log_vals{ii}, 'HorizontalAlignment', rtickHAlign, ...
                'VerticalAlignment', rtickVAlign, 'FontSize', 9);
        end

        title(ax, sprintf('Phases folded to 0-%d°, n=%d', period, sum(valid_idx)));
        tt = 0:40:320;
        set(ax, 'ThetaTick', tt);
        set(ax, 'ThetaTickLabel', arrayfun(@(x) sprintf('%d°', x/(360/period)), tt, 'UniformOutput', false));

        r_p = ax.RLim(2) * 1.03;
        p_str = sprintf('V test\np = %.3f', vtest_p);
        p_col = [0 0 0];
        if vtest_p < 0.05
            p_col = [1 0 0];
        end
        text(ax, 59/360*2*pi, r_p, p_str, 'Color', p_col, 'FontWeight', 'bold', 'FontSize', 8, ...
            'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');
        hold(ax, 'off');
    end
end

function plot_angle_curve_full_nosample(ax, data, centers, foldPeriod, CardColor, GapColor, ObliColor, yRange, alignTo90)
    if nargin < 9
        alignTo90 = false;  % Default: align to 0°
    end
    
    m = mean(data, 1, 'omitnan');
    se = std(data, 0, 1, 'omitnan') ./ sqrt(max(1, size(data,1)));

    cla(ax);
    hold(ax, 'on');
    fill(ax, [centers, fliplr(centers)], [m + se, fliplr(m - se)], [0.75 0.75 0.75], ...
        'EdgeColor', 'none', 'FaceAlpha', 0.35);
    plot(ax, centers, m, '-', 'Color', [0.3 0.3 0.3], 'LineWidth', 1.4);

    if nargin >= 8 && ~isempty(yRange)
        ylim(ax, yRange);
    end

    if foldPeriod == 360
        for a = 0:45:315
            if mod(a, 90) == 0
                xline(ax, a, '--', 'Color', CardColor, 'LineWidth', 0.8, 'HandleVisibility', 'off');
            else
                xline(ax, a, '--', 'Color', ObliColor, 'LineWidth', 0.8, 'HandleVisibility', 'off');
            end
            xline(ax, a + 22.5, '--', 'Color', GapColor, 'LineWidth', 0.7, 'HandleVisibility', 'off');
        end
        xticks(ax, 0:45:315);
        xticklabels(ax, arrayfun(@(v) sprintf('%.1f°', v), 0:45:315, 'UniformOutput', false));
        xlabel(ax, 'Angle Mod 360°');
    elseif alignTo90  % New configs: aligned to 90°, starting from 75°
        if foldPeriod == 120
            % 3fold with 90° alignment
            xline(ax, 90, '--', 'Color', CardColor, 'LineWidth', 0.9, 'HandleVisibility', 'off');
            xline(ax, 150, '--', 'Color', GapColor, 'LineWidth', 0.9, 'HandleVisibility', 'off');
            xticks(ax, [75, 90, 150]);
            xticklabels(ax, {'75.0°','90.0°','150.0°'});
            xlabel(ax, 'Angle (3fold, aligned to 90°)');
            xlim(ax, [75, 195]);
        elseif foldPeriod == 72
            % 5fold with 90° alignment
            xline(ax, 90, '--', 'Color', CardColor, 'LineWidth', 0.9, 'HandleVisibility', 'off');
            xline(ax, 126, '--', 'Color', GapColor, 'LineWidth', 0.9, 'HandleVisibility', 'off');
            xticks(ax, [75, 90, 126]);
            xticklabels(ax, {'75.0°','90.0°','126.0°'});
            xlabel(ax, 'Angle (5fold, aligned to 90°)');
            xlim(ax, [75, 147]);
        elseif foldPeriod == 60
            % 6fold with 90° alignment
            xline(ax, 90, '--', 'Color', CardColor, 'LineWidth', 0.9, 'HandleVisibility', 'off');
            xline(ax, 120, '--', 'Color', GapColor, 'LineWidth', 0.9, 'HandleVisibility', 'off');
            xticks(ax, [75, 90, 120]);
            xticklabels(ax, {'75.0°','90.0°','120.0°'});
            xlabel(ax, 'Angle (6fold, aligned to 90°)');
            xlim(ax, [75, 135]);
        elseif abs(foldPeriod - 360/7) < 0.1
            % 7fold with 90° alignment
            xline(ax, 90, '--', 'Color', CardColor, 'LineWidth', 0.9, 'HandleVisibility', 'off');
            half_fold = foldPeriod / 2;
            xline(ax, 90 + half_fold, '--', 'Color', GapColor, 'LineWidth', 0.9, 'HandleVisibility', 'off');
            xticks(ax, [75, 90, 90 + half_fold]);
            xticklabels(ax, {sprintf('%.1f°', 75), sprintf('%.1f°', 90), sprintf('%.1f°', 90 + half_fold)});
            xlabel(ax, sprintf('Angle (7fold, aligned to 90°)'));
            xlim(ax, [75, 90 + foldPeriod]);
        end
    elseif foldPeriod == 120
        % 3fold (old style, aligned to 0°)
        xline(ax, 0, '--', 'Color', CardColor, 'LineWidth', 0.9, 'HandleVisibility', 'off');
        xline(ax, 60, '--', 'Color', GapColor, 'LineWidth', 0.9, 'HandleVisibility', 'off');
        xticks(ax, [0, 60]);
        xticklabels(ax, {'0.0°','60.0°'});
        xlabel(ax, 'Angle Mod 120°');
    elseif foldPeriod == 90
        xline(ax, 0, '--', 'Color', CardColor, 'LineWidth', 0.9, 'HandleVisibility', 'off');
        xline(ax, 22.5, '--', 'Color', GapColor, 'LineWidth', 0.9, 'HandleVisibility', 'off');
        xline(ax, 45, '--', 'Color', ObliColor, 'LineWidth', 0.9, 'HandleVisibility', 'off');
        xline(ax, 67.5, '--', 'Color', GapColor, 'LineWidth', 0.9, 'HandleVisibility', 'off');
        xticks(ax, [0, 22.5, 45, 67.5]);
        xticklabels(ax, {'0.0°','22.5°','45.0°','67.5°'});
        xlabel(ax, 'Angle Mod 90°');

        % Restore pairwise significance bracket: 0° vs 45°.
        [~, idx0] = min(abs(centers - 0));
        [~, idx45] = min(abs(centers - 45));
        vals0 = data(:, idx0);
        vals45 = data(:, idx45);
        try
            [~, p_pair] = ttest(vals0, vals45);
        catch
            p_pair = NaN;
        end
        bump = 0.06 * range(ylim(ax));
        y1 = m(idx0) + se(idx0) + bump;
        y2 = m(idx45) + se(idx45) + bump;
        yb = max([y1, y2]) + bump * 0.45;
        yb = max(yb, 0.14);
        sigColor = [1, 0, 0];
        if ~isnan(p_pair) && p_pair >= 0.05 && p_pair < 0.10
            sigColor = [0.72, 0.43, 0.43];
        end
        plot(ax, [centers(idx0), centers(idx0), centers(idx45), centers(idx45)], [y1, yb, yb, y2], '-', 'Color', sigColor, 'LineWidth', 1.4);
        text(ax, mean([centers(idx0), centers(idx45)]), yb + 0.85*bump, sprintf('p=%.3f', p_pair), ...
            'HorizontalAlignment', 'center', 'Color', sigColor, 'FontWeight', 'bold', 'FontSize', 9);
    elseif foldPeriod == 72
        % 5fold (old style, aligned to 0°)
        xline(ax, 0, '--', 'Color', CardColor, 'LineWidth', 0.9, 'HandleVisibility', 'off');
        xline(ax, 36, '--', 'Color', GapColor, 'LineWidth', 0.9, 'HandleVisibility', 'off');
        xticks(ax, [0, 36]);
        xticklabels(ax, {'0.0°','36.0°'});
        xlabel(ax, 'Angle Mod 72°');
    elseif foldPeriod == 60
        % 6fold (old style, aligned to 0°)
        xline(ax, 0, '--', 'Color', CardColor, 'LineWidth', 0.9, 'HandleVisibility', 'off');
        xline(ax, 30, '--', 'Color', GapColor, 'LineWidth', 0.9, 'HandleVisibility', 'off');
        xticks(ax, [0, 30]);
        xticklabels(ax, {'0.0°','30.0°'});
        xlabel(ax, 'Angle Mod 60°');
    elseif abs(foldPeriod - 360/7) < 0.1
        % 7fold (old style, aligned to 0°)
        xline(ax, 0, '--', 'Color', CardColor, 'LineWidth', 0.9, 'HandleVisibility', 'off');
        half_fold = foldPeriod / 2;
        xline(ax, half_fold, '--', 'Color', GapColor, 'LineWidth', 0.9, 'HandleVisibility', 'off');
        xticks(ax, [0, half_fold]);
        xticklabels(ax, {sprintf('%.1f°', 0), sprintf('%.1f°', half_fold)});
        xlabel(ax, sprintf('Angle Mod %.1f°', foldPeriod));
    elseif foldPeriod == 45
        xline(ax, 0, '--', 'Color', CardColor, 'LineWidth', 0.9, 'HandleVisibility', 'off');
        xline(ax, 22.5, '--', 'Color', GapColor, 'LineWidth', 0.9, 'HandleVisibility', 'off');
        xticks(ax, [0, 22.5]);
        xticklabels(ax, {'0.0°','22.5°'});
        xlabel(ax, 'Angle Mod 45°');

        % Restore pairwise significance bracket: 0° vs 22.5°.
        [~, idx0] = min(abs(centers - 0));
        [~, idx225] = min(abs(centers - 22.5));
        vals0 = data(:, idx0);
        vals225 = data(:, idx225);
        try
            [~, p_pair] = ttest(vals0, vals225);
        catch
            p_pair = NaN;
        end
        bump = 0.06 * range(ylim(ax));
        y1 = m(idx0) + se(idx0) + bump;
        y2 = m(idx225) + se(idx225) + bump;
        yb = max([y1, y2]) + bump * 0.45;
        plot(ax, [centers(idx0), centers(idx0), centers(idx225), centers(idx225)], [y1, yb, yb, y2], 'r-', 'LineWidth', 1.4);
        text(ax, mean([centers(idx0), centers(idx225)]), yb + 0.85*bump, sprintf('p=%.3f', p_pair), ...
            'HorizontalAlignment', 'center', 'Color', 'r', 'FontWeight', 'bold', 'FontSize', 9);
    end

    if ~alignTo90 || (foldPeriod ~= 120 && foldPeriod ~= 72 && foldPeriod ~= 60 && abs(foldPeriod - 360/7) >= 0.1)
        xlim(ax, [min(centers), max(centers)]);
    end
    ylabel(ax, 'Count Proportion');
    set(ax, 'FontSize', 10, 'LineWidth', 1.0);
    box(ax, 'off');
    hold(ax, 'off');
end

function enforce_uniform_fontsize(figHandle, maxFs)
    if nargin < 2
        maxFs = 11;
    end

    axList = findall(figHandle, 'Type', 'axes');
    for i = 1:numel(axList)
        ax = axList(i);
        try
            ax.FontSize = min(ax.FontSize, maxFs);
        catch
        end
        try
            tt = get(ax, 'Title');
            if ~isempty(tt)
                tt.FontSize = maxFs;
            end
        catch
        end
        try
            xl = get(ax, 'XLabel');
            yl = get(ax, 'YLabel');
            zl = get(ax, 'ZLabel');
            if ~isempty(xl), xl.FontSize = min(xl.FontSize, maxFs); end
            if ~isempty(yl), yl.FontSize = min(yl.FontSize, maxFs); end
            if ~isempty(zl), zl.FontSize = min(zl.FontSize, maxFs); end
        catch
        end
    end

    cbList = findall(figHandle, 'Type', 'ColorBar');
    for i = 1:numel(cbList)
        try
            cbList(i).FontSize = min(cbList(i).FontSize, maxFs);
        catch
        end
    end

    txtList = findall(figHandle, 'Type', 'text');
    for i = 1:numel(txtList)
        try
            txtList(i).FontSize = min(txtList(i).FontSize, maxFs);
        catch
        end
    end
end

function shift_dprime_annotation(ax, dxNorm)
    if nargin < 2
        dxNorm = -0.12;
    end
    txt = findall(ax, 'Type', 'text');
    for i = 1:numel(txt)
        s = string(txt(i).String);
        if contains(lower(s), "d'")
            oldUnits = txt(i).Units;
            txt(i).Units = 'normalized';
            p = txt(i).Position;
            p(1) = p(1) + dxNorm;
            txt(i).Position = p;
            txt(i).Units = oldUnits;
        end
    end
end

function style_marginal_violin_sig(ax, pTarget, tol, dy)
    % Re-style and reposition marginal significance (e.g., p~0.058) on violin panels.
    if nargin < 2, pTarget = 0.058; end
    if nargin < 3, tol = 0.012; end
    if nargin < 4, dy = -0.010; end
    sigColor = [0.72, 0.43, 0.43];

    txt = findall(ax, 'Type', 'text');
    for i = 1:numel(txt)
        s = char(string(txt(i).String));
        tk = regexp(s, 'p\s*=\s*([0-9]*\.?[0-9]+)', 'tokens', 'once');
        if isempty(tk)
            continue;
        end
        pval = str2double(tk{1});
        if isnan(pval) || abs(pval - pTarget) > tol
            continue;
        end

        oldPos = txt(i).Position;
        txt(i).Position = [oldPos(1), oldPos(2) + dy, oldPos(3)];
        txt(i).Color = sigColor;
        txt(i).FontWeight = 'bold';

        % Shift the nearest upper bracket line together with the p-text.
        ln = findall(ax, 'Type', 'line');
        for j = 1:numel(ln)
            xd = ln(j).XData;
            yd = ln(j).YData;
            if isempty(xd) || isempty(yd)
                continue;
            end
            if max(yd) > oldPos(2) - 0.03 && max(yd) < oldPos(2) + 0.02
                ln(j).YData = yd + dy;
                ln(j).Color = sigColor;
                ln(j).LineWidth = max(1.2, ln(j).LineWidth);
            end
        end
    end
end

function set_dprime_right_margin(ax, xNorm)
    if nargin < 2
        xNorm = 0.90;
    end
    txt = findall(ax, 'Type', 'text');
    for i = 1:numel(txt)
        s = string(txt(i).String);
        if contains(lower(s), "d'")
            oldUnits = txt(i).Units;
            txt(i).Units = 'normalized';
            p = txt(i).Position;
            p(1) = xNorm;
            txt(i).Position = p;
            txt(i).HorizontalAlignment = 'right';
            txt(i).Units = oldUnits;
        end
    end
end

function remove_legends_from_figure(figHandle)
    lg = findall(figHandle, 'Type', 'Legend');
    if ~isempty(lg)
        delete(lg);
    end
end