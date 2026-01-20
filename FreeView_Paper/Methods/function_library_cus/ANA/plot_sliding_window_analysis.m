function plot_sliding_window_analysis(TseriesAxis, TseriesGap, TseriesCard, TseriesObli, xWin, cfg)
% PLOT_SLIDING_WINDOW_ANALYSIS 绘制滑动窗口分析结果
%
% 输入：
%   TseriesAxis, TseriesGap, TseriesCard, TseriesObli : nsbj × nWin 矩阵
%   xWin : 1 × nWin 窗口中心位置向量
%   cfg  - 配置结构体，包含：
%       .cmap16_FT : 16×3 颜色映射
%       .AxisColor : Axis 颜色（可选，默认混合 Card 和 Oblique 颜色）
%       .mode      : 'continuous' 或 'block'
%       .doStats   : 是否执行统计检验（默认 true）
%       .xlabel     : x轴标签
%       .ylabel     : y轴标签
%       .win_trials : 窗口大小（单位：trial数）
%       .Card_in_Axis : 下方子图绘制 Card 在 Axis 中的比例（true），
%                        还是 Card vs Oblique 对比（false，默认）
%       .cutoff    : 用于趋势检验的最大trial数（默认使用全部）

    % ---- 参数验证 ----
    if ~isfield(cfg, 'doStats'), cfg.doStats = true; end
    if ~isfield(cfg, 'cmap16_FT')
        cfg.cmap16_FT = lines(16);
    end
    if ~isfield(cfg, 'Card_in_Axis'), cfg.Card_in_Axis = false; end
    if ~isfield(cfg, 'AxisColor')
        cfg.AxisColor = mix_RGB_by_HSL(cfg.cmap16_FT(1,:), cfg.cmap16_FT(3,:));
    end
    if ~isfield(cfg, 'xlabel'), cfg.xlabel = 'Trial number (window center)'; end
    if ~isfield(cfg, 'ylabel'), cfg.ylabel = ''; end
    if ~isfield(cfg, 'mode'), cfg.mode = 'continuous'; end
    if ~isfield(cfg, 'win_trials'), cfg.win_trials = []; end
    if ~isfield(cfg, 'cutoff'), cfg.cutoff = size(TseriesAxis, 2); end
    % ---- 合并为一张图的上下两子图 ----
    figure;
    tl = tiledlayout(2,1, 'Padding', 'compact', 'TileSpacing', 'compact');

    % 上方：Axis proportion
    ax1 = nexttile(1);
    cfg1 = cfg;
    cfg1.axesHandle = ax1;
    cfg1.statTail = 'both'; % 检验是否显著大于0.5
    plot_single_prop(TseriesAxis, xWin, cfg.AxisColor, 'Axis proportion', cfg1);
    title(ax1, 'Axis Effect');
    xlabel(ax1, 'Trial number (window center)');
    ylabel(ax1, '');
    % 在右上角（标题右侧上方）标记窗口大小

    if cfg.Card_in_Axis
        % 下方：Card proportion
        ax2 = nexttile(2);
        cfg2 = cfg;
        cfg2.axesHandle = ax2;
        cfg2.statTail = 'right'; % 只检验Card > 0.5是否显著
        plot_single_prop(TseriesCard, xWin, cfg.cmap16_FT(1,:), 'Card proportion', cfg2);
        title(ax2, 'Card. Effect');
        xlabel(ax2, 'Trial number (window center)');
        ylabel(ax2, '');
    else
        % 下方：Card vs Oblique对比
        ax2 = nexttile(2);
        cfg2 = cfg;
        cfg2.axesHandle = ax2;
        cfg2.chanceLevel = TseriesGap./2; % 以Gap/2作为基线
        cfg2.chanceLabel = 'Gap/2';
        cfg2.statTail = 'both'; % 双尾检验Card vs Oblique
        plot_comparison(TseriesCard, TseriesObli, xWin, cfg.cmap16_FT([1,3],:), {'Card','Oblique'}, cfg2);
        title(ax2, '');
    end
    % 共用y轴标签，横跨上下两个子图
    xlFontSize = get(ax1, 'FontSize');
    ylabel(tl, cfg.ylabel, 'FontSize', xlFontSize, 'FontWeight', 'normal');
    
    % ---- Mann-Kendall 趋势检验 ----
    disp('Mann-Kendall trend test results:');
    varNames = {'TseriesAxis'};%, 'TseriesGap', 'TseriesCard', 'TseriesObli'};
    dataVars = {TseriesAxis};%, TseriesGap, TseriesCard, TseriesObli};
    for i = 1:numel(dataVars)
        % 对每个变量，先对被试取均值（忽略NaN），再做趋势检验
        meanSeries = mean(dataVars{i}, 1, 'omitnan');
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
        fprintf('%s: %s (p = %.4f, Z = %.3f)\n', varNames{i}, trendStr, p, Z);
    end
end

