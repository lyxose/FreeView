function plot_sliding_window_analysis_corr(TseriesAxis, TseriesGap, TseriesCard, TseriesObli, xWin, cfg)
% 上图：Axis-Gap（v1_5时为Gap-Axis，需翻转），下图：Card-Obli
% 并计算两组数据的相关性（Pearson）
% 两个子图都做显著性检验（是否显著大于0）

    if ~isfield(cfg, 'doStats'), cfg.doStats = true; end
    if ~isfield(cfg, 'cmap16_FT'), cfg.cmap16_FT = lines(16); end
    if ~isfield(cfg, 'ver'), cfg.ver = ''; end

    % Axis-Gap或Gap-Axis（根据ver翻转）
    if strcmpi(cfg.ver, 'v1_5')
        ag_series = TseriesGap - TseriesAxis;
        ag_label = 'Gap - Axis';
    else
        ag_series = TseriesAxis - TseriesGap;
        ag_label = 'Axis - Gap';
    end

    % Card-Obli
    co_series = TseriesCard - TseriesObli;
    co_label = 'Card - Obli';

    figure;
    tl = tiledlayout(2,1, 'Padding', 'compact', 'TileSpacing', 'compact');

    % 上方：Axis-Gap或Gap-Axis
    ax1 = nexttile(1);
    cfg1 = cfg;
    cfg1.axesHandle = ax1;
    cfg1.statTail = 'right'; % 显著大于0
    plot_single(ag_series, xWin, cfg.cmap16_FT(1,:), ag_label, cfg1);
    title(ax1, ag_label);
    xlabel(ax1, 'Trial number (window center)');
    ylabel(ax1, '');

    xl = xlim(ax1); yl = ylim(ax1);
    x_text = xl(2);
    y_text = yl(2) + 0.03 * (yl(2) - yl(1));
    txt_str = sprintf('Win. size: %d', cfg.win_trials);
    text(ax1, x_text, y_text, txt_str, ...
        'HorizontalAlignment', 'right', ...
        'VerticalAlignment', 'bottom', ...
        'FontSize', max(8, get(ax1, 'FontSize')-2), ...
        'FontWeight', 'normal', ...
        'Color', [0.2 0.2 0.2]);

    % 下方：Card-Obli
    ax2 = nexttile(2);
    cfg2 = cfg;
    cfg2.axesHandle = ax2;
    cfg2.statTail = 'right'; % 显著大于0
    plot_single(co_series, xWin, cfg.cmap16_FT(3,:), co_label, cfg2);
    title(ax2, co_label);
    xlabel(ax2, 'Trial number (window center)');
    ylabel(ax2, '');

    % 共用y轴标签
    xlFontSize = get(ax1, 'FontSize');
    ylabel(tl, cfg.ylabel, 'FontSize', xlFontSize, 'FontWeight', 'normal');

    % 计算相关性（组均值）
    mean_ag = mean(ag_series, 1, 'omitnan');
    mean_co = mean(co_series, 1, 'omitnan');
    valid_mask = isfinite(mean_ag) & isfinite(mean_co);
    if nnz(valid_mask) > 2
        [r, p] = corr(mean_ag(valid_mask)', mean_co(valid_mask)', 'type', 'Pearson');
        fprintf('Correlation between %s and %s: r = %.3f, p = %.4g\n', ag_label, co_label, r, p);
        % 在下方子图右上方（外部）标注相关性
        ax2_pos = get(ax2, 'Position');
        ann_x = ax2_pos(1) + ax2_pos(3) + 0.04;
        ann_y = ax2_pos(2) + ax2_pos(4) + 0.2;
        annotation(gcf, 'textbox', [ann_x, ann_y, 0.22, 0.08], ...
            'String', sprintf('Corr: r = %.3f\np = %.4g', r, p), ...
            'FontSize', 12, 'FontWeight', 'bold', ...
            'EdgeColor', 'none', 'Color', [0.2 0.2 0.2], ...
            'HorizontalAlignment', 'right', 'VerticalAlignment', 'top');
    end
end
