function plot_sliding_window_analysis(TseriesAxis, TseriesGap, TseriesCard, TseriesObli, xWin, cfg)
% PLOT_SLIDING_WINDOW_ANALYSIS 绘制滑动窗口分析结果
%
% 输入：
%   TseriesAxis, TseriesGap, TseriesCard, TseriesObli : nsbj × nWin 矩阵
%   xWin : 1 × nWin 窗口中心位置向量
%   cfg  - 配置结构体，包含：
%       .cmap16_FT : 16×3 颜色映射
%       .mode      : 'continuous' 或 'block'
%       .doStats   : 是否执行统计检验（默认 true）
%       .xlabel     : x轴标签
%       .ylabel     : y轴标签
%       .win_trials : 窗口大小（单位：trial数）

    % ---- 参数验证 ----
    if ~isfield(cfg, 'doStats'), cfg.doStats = true; end
    if ~isfield(cfg, 'cmap16_FT')
        cfg.cmap16_FT = lines(16);
    end

    % ---- 合并为一张图的上下两子图 ----
    figure;
    tl = tiledlayout(2,1, 'Padding', 'compact', 'TileSpacing', 'compact');

    % 上方：Axis proportion
    ax1 = nexttile(1);
    cfg1 = cfg;
    cfg1.axesHandle = ax1;
    cfg1.statTail = 'both'; % 检验是否显著大于0.5
    plot_single_prop(TseriesAxis, xWin, cfg.cmap16_FT(1,:), 'Axis proportion', cfg1);
    title(ax1, 'Axis Proportion');
    xlabel(ax1, 'Trial number (window center)');
    ylabel(ax1, '');
    % 在右上角（标题右侧上方）标记窗口大小
    xl = xlim(ax1); yl = ylim(ax1);
    x_text = xl(2);
    y_text = yl(2) + 0.03 * (yl(2) - yl(1));
    txt_str = sprintf('Win. size: %d', cfg.win_trials);
    text(ax1, x_text, y_text, txt_str, ...
        'HorizontalAlignment', 'left', ...
        'VerticalAlignment', 'bottom', ...
        'FontSize', max(8, get(ax1, 'FontSize')-2), ...
        'FontWeight', 'normal', ...
        'Color', [0.2 0.2 0.2]);

    % 下方：Card proportion
    ax2 = nexttile(2);
    cfg2 = cfg;
    cfg2.axesHandle = ax2;
    cfg2.statTail = 'right'; % 只检验Card > 0.5是否显著
    plot_single_prop(TseriesCard, xWin, cfg.cmap16_FT(3,:), 'Card proportion', cfg2);
    title(ax2, 'Card Proportion');
    xlabel(ax2, 'Trial number (window center)');
    ylabel(ax2, '');

    % 共用y轴标签，横跨上下两个子图
    xlFontSize = get(ax1, 'FontSize');
    ylabel(tl, cfg.ylabel, 'FontSize', xlFontSize, 'FontWeight', 'normal');

end
