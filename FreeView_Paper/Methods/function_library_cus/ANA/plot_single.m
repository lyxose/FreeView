function plot_single(series, xWin, color, label, cfg)
    % 绘制单一时程数据，进行与0的显著性比较
    % 输入：
    % series: nSubj x nTimepoints
    % xWin: 1 x nTimepoints
    % color: 绘图颜色
    % label: 图例标签
    % cfg: 配置结构体，包含以下字段：
    %   doStats: 是否进行统计检验（布尔值）
    %   statTail: 统计检验的尾部选项（'both', 'left', 'right'）
    %   axesHandle: （可选）指定绘图的坐标轴句柄

    nsbj = size(series, 1);
    m = mean(series, 1, 'omitnan');
    se = std(series, 0, 1, 'omitnan') / sqrt(nsbj);

    stat = struct;
    if cfg.doStats
        try
            % 支持cfg.statTail指定统计方向
            if isfield(cfg, 'statTail') && ~isempty(cfg.statTail)
                tailOpt = cfg.statTail;
            else
                tailOpt = 'both';
            end
            if exist('cluster_based_permu_independent', 'file') == 2
                [stat, ~] = cluster_based_permu_independent(series, zeros(size(series))); % 默认双尾
            else
                [stat, ~] = cluster_based_permutation_test(series, zeros(size(series)), 'tail', tailOpt);
            end
        catch ME
            warning(ME.identifier, '%s', ME.message);
            [stat, ~] = cluster_based_permutation_test(series, zeros(size(series)), 'tail', tailOpt);
        end
    end

    if isfield(cfg, 'axesHandle') && ~isempty(cfg.axesHandle) && ishandle(cfg.axesHandle)
        ax = cfg.axesHandle;
        axes(ax); cla(ax); hold(ax, 'on');
    else
        figure; ax = gca; hold(ax, 'on');
    end

    fill(ax, [xWin, fliplr(xWin)], [m+se, fliplr(m-se)], color, ...
         'FaceAlpha', 0.2, 'EdgeColor', 'none');
    plot(ax, xWin, m, '-', 'Color', color, 'LineWidth', 2);

    if isfield(stat, 'prob')
        sig = stat.prob < 0.05;
        if any(sig)
            d = diff([0, sig, 0]);
            sIdx = find(d == 1); eIdx = find(d == -1) - 1;
            yl = ylim(ax); y0 = yl(1) + 0.03 * (yl(2) - yl(1));
            for kk = 1:numel(sIdx)
                plot(ax, [xWin(sIdx(kk)), xWin(eIdx(kk))], [y0, y0], 'r-', 'LineWidth', 4);
            end
        end
    end

    xlabel(ax, cfg.xlabel);
    ylabel(ax, cfg.ylabel);
    legend(ax, {[label,'±SE'], label}, 'Location', 'northeastoutside');
    set(ax, 'FontSize', 12, 'LineWidth', 1.2); box off; hold(ax, 'off');
end
