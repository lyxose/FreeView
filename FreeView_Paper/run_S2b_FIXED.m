%% FIXED S2b Script - Comprehensive version with all corrections
% This is a complete replacement script that incorporates all fixes:
% 1. Correct startAngle = 75 for aligned_to_90 folds
% 2. Dynamic binning precision for small periods
% 3. Peak/valley value annotations

function run_S2b_fold_symmetries_fixed(whichVersion)
if nargin < 1 || isempty(whichVersion)
    whichVersion = 'v2';
end

fprintf('[INFO] Starting S2b generation with ALL FIXES (version: %s)\n', whichVersion);

close all; clc;

%% Parameters
rootDir = fileparts(mfilename('fullpath'));
addpath(fullfile(rootDir, 'Methods', 'function_library_cus', 'ANA'));

outDir = fullfile(rootDir, 'Results', 'Paper_Composites_PDF');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

CardColor = [61 124 76] / 255;
ObliColor = [69 87 191] / 255;
GapColor  = [121 137 170] / 255;

n_bin_FT = 16;
shift_FT = 360 / n_bin_FT / 2;
edges_FT = linspace(0, 360, n_bin_FT + 1) - shift_FT;
angbinSize = 11.25;
keep_Time = 4000;
timeRes_FT = 10;
R_min = 1;
R_max = 7.5;
win_left = 1000;
win_right = 4000;

%% Load data
versMetadata = struct();
versMetadata.v1 = struct('path', 'v1', 'file', 'ALL_fixTable_26Subj.csv');
versMetadata.v1_5 = struct('path', 'v1.5', 'file', 'ALL_fixTable_26Subj.csv');
versMetadata.v2 = struct('path', 'v2', 'file', 'ALL_fixTable_25Subj.csv');

if whichVersion == "v1" || whichVersion == "v1.5"
    metaKey = strrep(whichVersion, '.', '_');
else
    metaKey = 'v2';
end

versMeta = versMetadata.(metaKey);
csvPath = fullfile(rootDir, 'Results', versMeta.path, versMeta.file);

fprintf('[INFO] Loading %s data from: %s\n', whichVersion, csvPath);
D = struct();
D.(metaKey) = preprocess_from_fixTable(readtable(csvPath));

fprintf('[OK] Data loaded (%d fixations)\n', height(D.(metaKey).visTable));

%% Generate S2b figure
fprintf('[INFO] Generating S2b figure...\n');

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

dd = D.(metaKey);

for i = 1:numel(fold_configs)
    cfg = fold_configs{i};
    ax = nexttile(tls2b, i);
    
    fprintf('[INFO]   Generating %s panel...\n', cfg.name);
    
    scan = compute_angle_scan_curves_fixed(dd.angles(dd.winSel), dd.sub(dd.winSel), dd.ses(dd.winSel), dd.pairs, dd.Nsubj, angbinSize);
    
    centers_fold = scan.(cfg.field_centers);
    sub_fold = scan.(cfg.field_sub);
    
    if cfg.align_90
        plot_angle_curve_fixed(ax, sub_fold, centers_fold, cfg.period, CardColor, GapColor, ObliColor, [], true);
        title(ax, sprintf('%s (%.1f fold, aligned 90)', cfg.name, cfg.period));
    else
        plot_angle_curve_fixed(ax, sub_fold, centers_fold, cfg.period, CardColor, GapColor, ObliColor, [], false);
        title(ax, sprintf('%s (%.1f fold)', cfg.name, cfg.period));
    end
    add_panel_label(ax, char('A' + i - 1));
end

enforce_uniform_fontsize(fs2b, 11);

fprintf('[INFO] Exporting S2b figure (%s)...\n', whichVersion);
outFileName = sprintf('Supplemental_Figure_S2b_Additional_Fold_Symmetries_%s_FIXED.pdf', strrep(string(whichVersion), '.', '_'));
export_pdf_vector(fs2b, fullfile(outDir, outFileName));

fprintf('[OK] S2b figure exported to: %s\n', outDir);
fprintf('[INFO] Files generated:\n');
fprintf('     - %s\n', outFileName);
fprintf('[OK] COMPLETE!\n');

% ========================= Helper functions =========================

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
        bin_counts(si, :) = histcounts(d.angles(m), edges_FT);
    end
    d.count360 = bin_counts;
end

function scan = compute_angle_scan_curves_fixed(angles, sub, ses, pairs, nsbj, angbinSize)
    startAngle = -22.5/2;
    [centers360, ~, ~, sub360] = analyze_angle_curve(angles, sub, ses, pairs, nsbj, angbinSize, 360, 'sum1', startAngle);

    % FIXED: fold_helper_90 now uses startAngle = 75 (not 90 - period/2)
    fold_helper_0 = @(period) fold_angle_data_fixed(centers360, sub360, nsbj, startAngle, period);
    fold_helper_90 = @(period) fold_angle_data_fixed(centers360, sub360, nsbj, 75, period);  % FIXED: was 90 - period/2
    
    [centers120, sub120] = fold_helper_90(120);
    [centers90, sub90] = fold_helper_0(90);
    [centers72, sub72] = fold_helper_90(72);
    [centers60, sub60] = fold_helper_90(60);
    [centers51_43, sub51_43] = fold_helper_90(360/7);
    [centers45, sub45] = fold_helper_0(45);

    scan = struct('centers360', centers360, 'sub360', sub360, ...
                  'centers45', centers45, 'sub45', sub45, ...
                  'centers90', centers90, 'sub90', sub90, ...
                  'centers120', centers120, 'sub120', sub120, ...
                  'centers72', centers72, 'sub72', sub72, ...
                  'centers60', centers60, 'sub60', sub60, ...
                  'centers51_43', centers51_43, 'sub51_43', sub51_43);
end

function [centers_folded, sub_folded] = fold_angle_data_fixed(centers360, sub360, nsbj, startAngle, foldPeriod)
    mod_angles = mod(centers360 - startAngle, foldPeriod);
    
    % FIXED: Adjust binning precision for small periods (like 7fold ~51.43 degrees)
    if foldPeriod < 55
        roundPrec = 6;  % Coarser binning to reduce noise for ultra-small periods
    else
        roundPrec = 8;
    end
    
    [uniq_mod, ~, ic] = unique(round(mod_angles, roundPrec));
    sub_folded = zeros(nsbj, numel(uniq_mod));
    for k = 1:numel(uniq_mod)
        idx = (ic == k);
        sub_folded(:,k) = sum(sub360(:,idx), 2, 'omitnan');
    end
    centers_folded = uniq_mod + startAngle;
end

function plot_angle_curve_fixed(ax, data, centers, foldPeriod, CardColor, GapColor, ObliColor, yRange, alignTo90)
    if nargin < 9
        alignTo90 = false;
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

    if alignTo90
        if foldPeriod == 120
            % 3fold
            xline(ax, 90, '--', 'Color', CardColor, 'LineWidth', 0.9, 'HandleVisibility', 'off');
            xline(ax, 150, '--', 'Color', GapColor, 'LineWidth', 0.9, 'HandleVisibility', 'off');
            [~, idx90] = min(abs(centers - 90));
            [~, idx150] = min(abs(centers - 150));
            val90 = m(idx90);
            val150 = m(idx150);
            text(ax, 90, val90 + 0.002, sprintf('%.4f', val90), 'HorizontalAlignment', 'center', 'FontSize', 9, 'Color', CardColor);
            text(ax, 150, val150 - 0.002, sprintf('%.4f', val150), 'HorizontalAlignment', 'center', 'FontSize', 9, 'Color', GapColor);
            xticks(ax, [75, 90, 150]);
            xticklabels(ax, {'75.0 deg','90.0 deg','150.0 deg'});
            xlabel(ax, 'Angle (3fold, aligned to 90 deg)');
            xlim(ax, [75, 195]);
            
        elseif foldPeriod == 72
            % 5fold
            xline(ax, 90, '--', 'Color', CardColor, 'LineWidth', 0.9, 'HandleVisibility', 'off');
            xline(ax, 126, '--', 'Color', GapColor, 'LineWidth', 0.9, 'HandleVisibility', 'off');
            [~, idx90] = min(abs(centers - 90));
            [~, idx126] = min(abs(centers - 126));
            val90 = m(idx90);
            val126 = m(idx126);
            text(ax, 90, val90 + 0.0002, sprintf('%.4f', val90), 'HorizontalAlignment', 'center', 'FontSize', 9, 'Color', CardColor);
            text(ax, 126, val126 - 0.0002, sprintf('%.4f', val126), 'HorizontalAlignment', 'center', 'FontSize', 9, 'Color', GapColor);
            xticks(ax, [75, 90, 126]);
            xticklabels(ax, {'75.0 deg','90.0 deg','126.0 deg'});
            xlabel(ax, 'Angle (5fold, aligned to 90 deg)');
            xlim(ax, [75, 147]);
            
        elseif foldPeriod == 60
            % 6fold
            xline(ax, 90, '--', 'Color', CardColor, 'LineWidth', 0.9, 'HandleVisibility', 'off');
            xline(ax, 120, '--', 'Color', GapColor, 'LineWidth', 0.9, 'HandleVisibility', 'off');
            [~, idx90] = min(abs(centers - 90));
            [~, idx120] = min(abs(centers - 120));
            val90 = m(idx90);
            val120 = m(idx120);
            text(ax, 90, val90 + 0.001, sprintf('%.4f', val90), 'HorizontalAlignment', 'center', 'FontSize', 9, 'Color', CardColor);
            text(ax, 120, val120 - 0.001, sprintf('%.4f', val120), 'HorizontalAlignment', 'center', 'FontSize', 9, 'Color', GapColor);
            xticks(ax, [75, 90, 120]);
            xticklabels(ax, {'75.0 deg','90.0 deg','120.0 deg'});
            xlabel(ax, 'Angle (6fold, aligned to 90 deg)');
            xlim(ax, [75, 135]);
            
        elseif abs(foldPeriod - 360/7) < 0.1
            % 7fold - FIXED: with better binning precision
            half_fold = foldPeriod / 2;
            xline(ax, 90, '--', 'Color', CardColor, 'LineWidth', 0.9, 'HandleVisibility', 'off');
            xline(ax, 90 + half_fold, '--', 'Color', GapColor, 'LineWidth', 0.9, 'HandleVisibility', 'off');
            [~, idx90] = min(abs(centers - 90));
            [~, idxGap] = min(abs(centers - (90 + half_fold)));
            val90 = m(idx90);
            valGap = m(idxGap);
            text(ax, 90, val90 + 0.0003, sprintf('%.4f', val90), 'HorizontalAlignment', 'center', 'FontSize', 8, 'Color', CardColor);
            text(ax, 90 + half_fold, valGap - 0.0003, sprintf('%.4f', valGap), 'HorizontalAlignment', 'center', 'FontSize', 8, 'Color', GapColor);
            xticks(ax, [75, 90, 90 + half_fold]);
            xticklabels(ax, {sprintf('%.1f deg', 75), sprintf('%.1f deg', 90), sprintf('%.1f deg', 90 + half_fold)});
            xlabel(ax, sprintf('Angle (7fold, aligned to 90 deg)'));
            xlim(ax, [75, 90 + foldPeriod]);
        end
    else
        if foldPeriod == 90
            xline(ax, 0, '--', 'Color', CardColor, 'LineWidth', 0.9, 'HandleVisibility', 'off');
            xline(ax, 22.5, '--', 'Color', GapColor, 'LineWidth', 0.9, 'HandleVisibility', 'off');
            xline(ax, 45, '--', 'Color', ObliColor, 'LineWidth', 0.9, 'HandleVisibility', 'off');
            xline(ax, 67.5, '--', 'Color', GapColor, 'LineWidth', 0.9, 'HandleVisibility', 'off');
            xticks(ax, [0, 22.5, 45, 67.5]);
            xticklabels(ax, {'0.0 deg','22.5 deg','45.0 deg','67.5 deg'});
            xlabel(ax, 'Angle Mod 90 deg');
        elseif foldPeriod == 45
            xline(ax, 0, '--', 'Color', CardColor, 'LineWidth', 0.9, 'HandleVisibility', 'off');
            xline(ax, 22.5, '--', 'Color', GapColor, 'LineWidth', 0.9, 'HandleVisibility', 'off');
            xticks(ax, [0, 22.5]);
            xticklabels(ax, {'0.0 deg','22.5 deg'});
            xlabel(ax, 'Angle Mod 45 deg');
        end
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
            set(get(ax, 'XLabel'), 'FontSize', maxFs);
            set(get(ax, 'YLabel'), 'FontSize', maxFs);
            set(get(ax, 'Title'), 'FontSize', maxFs);
        catch
        end
    end
    txList = findall(figHandle, 'Type', 'text');
    for i = 1:numel(txList)
        tx = txList(i);
        try
            set(tx, 'FontSize', min(get(tx, 'FontSize'), maxFs));
        catch
        end
    end
end

function export_pdf_vector(figHandle, outPdf)
    set(figHandle, 'Renderer', 'painters');
    exportgraphics(figHandle, outPdf, 'ContentType', 'vector');
    exportgraphics(figHandle, replace(outPdf, '.pdf', '.png'), 'Resolution', 600);
    savefig(figHandle, replace(outPdf, '.pdf', '.fig'));
    close(figHandle);
end

function add_panel_label(ax, labelText)
    text(ax, -0.12, 1.08, labelText, 'Units', 'normalized', ...
        'FontWeight', 'bold', 'FontSize', 13, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');
end

end
