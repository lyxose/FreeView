function plot_sliding_window_analysis_obli(TseriesObli, xWin, cfg)
% 绘制Oblique比例滑动窗口分析结果，并与0.25的chance比较
% TseriesObli: nsbj × nWin
% xWin: 1 × nWin
% cfg: 配置结构体，需包含 .cmap16_FT, .win_trials, .ylabel, .xlabel, .cutoff

    if ~isfield(cfg, 'doStats'), cfg.doStats = true; end
    if ~isfield(cfg, 'cmap16_FT'), cfg.cmap16_FT = lines(16); end
    if ~isfield(cfg, 'ylabel'), cfg.ylabel = 'Proportion'; end
    if ~isfield(cfg, 'xlabel'), cfg.xlabel = 'Trial number (window center)'; end
    if ~isfield(cfg, 'win_trials'), cfg.win_trials = []; end
    if ~isfield(cfg, 'cutoff'), cfg.cutoff = size(TseriesObli, 2); end

    % 设置图窗大小，保证坐标轴区域为11:5比例
    fig = figure;
    ax = gca;
    % 计算合适的figure和axes位置
    fig_w = 570; fig_h = 210;

    % 配置plot_single_prop参数
    cfg1 = cfg;
    cfg1.axesHandle = ax;
    cfg1.chanceLevel = 0.25;
    cfg1.chanceLabel = 'Chance';
    cfg1.statTail = 'right'; % 检验是否显著大于0.25
    
    plot_single_prop(TseriesObli, xWin, cfg.cmap16_FT(3,:), 'Oblique proportion', cfg1);
    set(fig, 'Units', 'pixels', 'Position', [100 100 fig_w fig_h]);

    title(ax, 'Oblique Proportion');
    xlabel(ax, cfg.xlabel);
    ylabel(ax, cfg.ylabel);

    % 标记窗口大小
    if ~isempty(cfg.win_trials)
        xl = xlim(ax); yl = ylim(ax);
        x_text = xl(2);
        y_text = yl(2) + 0.03 * (yl(2) - yl(1));
        txt_str = sprintf('Win. size: %d', cfg.win_trials);
        text(ax, x_text, y_text, txt_str, ...
            'HorizontalAlignment', 'left', ...
            'VerticalAlignment', 'bottom', ...
            'FontSize', max(8, get(ax, 'FontSize')-2), ...
            'FontWeight', 'normal', ...
            'Color', [0.2 0.2 0.2]);
    end

    % Mann-Kendall趋势检验
    disp('Mann-Kendall trend test for TseriesObli:');
    meanSeries = mean(TseriesObli, 1, 'omitnan');
    [h, p, S, Z, sen, n_eff, acf1] = mktrend_mod(meanSeries(1:cfg.cutoff));
    fprintf('Sen''s slope: %.6f\n', sen);
    if h == 1
        if Z > 0
            trendStr = 'monotonically increasing';
        else
            trendStr = 'monotonically decreasing';
        end
    else
        trendStr = 'no significant trend';
    end
    fprintf('TseriesObli: %s (p = %.4f, Z = %.3f)\n', trendStr, p, Z);
end
