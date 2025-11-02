function plot_comparison(seriesA, seriesB, xWin, colors, labels, cfg)
    % 绘制两组时程数据对比，进行组间显著性比较
    % 输入：
    % seriesA, seriesB: nSubj x nTimepoints
    % xWin: 1 x nTimepoints
    % colors: 2 x 3 矩阵，分别为两组数据的绘图颜色
    % labels: cell array，包含两组数据的标签
    % cfg: 配置结构体，包含以下字段：
    %   doStats: 是否进行统计检验（布尔值）
    %   statTail: 统计检验的尾部选项（'both', 'left', 'right'）
    %   axesHandle: （可选）指定绘图的坐标轴句柄
    nsbj = size(seriesA, 1);
    mA = mean(seriesA, 1, 'omitnan'); seA = std(seriesA, 0, 1, 'omitnan') / sqrt(nsbj);
    mB = mean(seriesB, 1, 'omitnan'); seB = std(seriesB, 0, 1, 'omitnan') / sqrt(nsbj);

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
                [stat, ~] = cluster_based_permu_independent(seriesA, seriesB); % 默认双尾
            else
                [stat, ~] = cluster_based_permutation_test(seriesA, seriesB, 'tail', tailOpt);
            end
        catch ME
            warning(ME.identifier, '%s', ME.message);
            [stat, ~] = cluster_based_permutation_test(seriesA, seriesB, 'tail', tailOpt);
        end
    end

    if isfield(cfg, 'axesHandle') && ~isempty(cfg.axesHandle) && ishandle(cfg.axesHandle)
        ax = cfg.axesHandle;
        axes(ax); cla(ax); hold(ax, 'on');
    else
        figure; ax = gca; hold(ax, 'on');
    end

    fill(ax, [xWin, fliplr(xWin)], [mA+seA, fliplr(mA-seA)], colors(1,:), ...
         'FaceAlpha', 0.2, 'EdgeColor', 'none');
    fill(ax, [xWin, fliplr(xWin)], [mB+seB, fliplr(mB-seB)], colors(2,:), ...
         'FaceAlpha', 0.2, 'EdgeColor', 'none');
    plot(ax, xWin, mA, '-', 'Color', colors(1,:), 'LineWidth', 2);
    plot(ax, xWin, mB, '-', 'Color', colors(2,:), 'LineWidth', 2);

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

    ylabel(ax, cfg.ylabel);
    xlabel(ax, cfg.xlabel);

    legend(ax, {[labels{1},'±SE'], [labels{2},'±SE'], labels{1}, labels{2}}, ...
            'Location', 'northeastoutside');
    set(ax, 'FontSize', 12, 'LineWidth', 1.2); box off; hold(ax, 'off');
end
