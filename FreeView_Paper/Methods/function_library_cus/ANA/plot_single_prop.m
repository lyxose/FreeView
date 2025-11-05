function plot_single_prop(series, xWin, color, label, cfg)
    % 绘制占比数据（0~1），与chance level（默认0.5）比较显著性
    % 其他绘图形式与 plot_single 保持一致

    if ~isfield(cfg, 'chanceLevel') || isempty(cfg.chanceLevel)
        chanceLevel = 0.5;
    else
        chanceLevel = cfg.chanceLevel;
    end
    
    if ~isfield(cfg, 'chanceLabel') || isempty(cfg.chanceLabel)
        cfg.chanceLabel = 'Chance';
    end

    nsbj = size(series, 1);
    m = mean(series, 1, 'omitnan');
    se = std(series, 0, 1, 'omitnan') / sqrt(nsbj);

    stat = struct;
    if cfg.doStats
        try
            if isfield(cfg, 'statTail') && ~isempty(cfg.statTail)
                tailOpt = cfg.statTail;
            else
                tailOpt = 'both';
            end
            % 与chanceLevel比较
            [stat, ~] = cluster_based_permutation_test(series, chanceLevel*ones(size(series)), 'tail', tailOpt);
        catch ME
            warning(ME.identifier, '%s', ME.message);
            [stat, ~] = cluster_based_permutation_test(series, chanceLevel*ones(size(series)), 'tail', tailOpt);
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

    % 绘制chance level基线
    yline(ax, chanceLevel, '--', 'Color', [0.4 0.4 0.4], 'LineWidth', 1.2, 'Label', cfg.chanceLabel, 'LabelHorizontalAlignment','left');
    y0=chanceLevel;
    if isfield(stat, 'prob')
        sig = stat.prob < 0.05;
        if any(sig)
            d = diff([0, sig, 0]);
            sIdx = find(d == 1); eIdx = find(d == -1) - 1;
            yl = ylim(ax); y0 = yl(1) + 0.03 * (yl(2) - yl(1));
            for kk = 1:numel(sIdx)
                plot(ax, [xWin(sIdx(kk)), xWin(eIdx(kk))], [y0, y0], 'r-', 'LineWidth', 3);
            end
        end
    end

    xlabel(ax, cfg.xlabel);
    ylabel(ax, cfg.ylabel);
    legend(ax, {[label,'±SE'], label, cfg.chanceLabel}, 'Location', 'northeastoutside');
    set(ax, 'FontSize', 12, 'LineWidth', 1.2);

    % 自动调整y轴范围
    y_min = min([chanceLevel, min(m-se), y0]);
    y_max = max([chanceLevel, max(m+se), y0]);
    y_pad = 0.08 * (y_max - y_min + eps);
    ylim(ax, [max(0, y_min - y_pad), min(1, y_max + y_pad)]);

    box off; hold(ax, 'off');
end
