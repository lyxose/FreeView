function plot_comparison(seriesA, seriesB, xWin, colors, labels, cfg)
    % 绘制两组时程数据对比，进行组间显著性比较，并可选与基线比较
    % 支持cfg.chanceLevel和cfg.chanceLabel指定基线
    % 新增：cfg.chanceLevel可为数组(与seriesA/B同形状)，此时绘制灰色均值±SE并与A/B做配对cluster-based permutation test
    nsbj = size(seriesA, 1);
    mA = mean(seriesA, 1, 'omitnan'); seA = std(seriesA, 0, 1, 'omitnan') / sqrt(nsbj);
    mB = mean(seriesB, 1, 'omitnan'); seB = std(seriesB, 0, 1, 'omitnan') / sqrt(nsbj);

    % 处理基线参数
    drawBaseline = false;
    isChanceArray = false;
    chanceArr = [];
    if isfield(cfg, 'chanceLevel') && ~isempty(cfg.chanceLevel)
        chanceLevel = cfg.chanceLevel;
        drawBaseline = true;
        if ~isscalar(chanceLevel)
            % 数组基线路径：要求与A/B同形状
            if ~isequal(size(chanceLevel), size(seriesA)) || ~isequal(size(chanceLevel), size(seriesB))
                error('cfg.chanceLevel为数组时，尺寸必须与seriesA和seriesB一致。');
            end
            isChanceArray = true;
            chanceArr = chanceLevel;
            mC = mean(chanceArr, 1, 'omitnan');
            seC = std(chanceArr, 0, 1, 'omitnan') / sqrt(nsbj);
        end
    else
        chanceLevel = 0;
    end
    if isfield(cfg, 'chanceLabel') && ~isempty(cfg.chanceLabel)
        chanceLabel = cfg.chanceLabel;
        drawBaseline = true;
    else
        chanceLabel = 'Chance';
    end

    % 统计检验
    stat = struct; statA = struct; statB = struct;
    if cfg.doStats
        if isfield(cfg, 'statTail') && ~isempty(cfg.statTail)
            tailOpt = cfg.statTail;
        else
            tailOpt = 'both';
        end
        % 组间独立
        [stat, ~] = cluster_based_permutation_test(seriesA, seriesB, 'tail', tailOpt);
        % 与基线比较
        if drawBaseline
            if isChanceArray
                % 配对检验
                [statA, ~] = cluster_based_permutation_test(seriesA, chanceArr, 'tail', tailOpt);
                [statB, ~] = cluster_based_permutation_test(seriesB, chanceArr, 'tail', tailOpt);
            else
                % 标量基线，保持原有方式
                baseA = chanceLevel * ones(size(seriesA));
                baseB = chanceLevel * ones(size(seriesB));
                [statA, ~] = cluster_based_permutation_test(seriesA, baseA, 'tail', tailOpt);
                [statB, ~] = cluster_based_permutation_test(seriesB, baseB, 'tail', tailOpt);
            end
        end
    end


    % 绘图
    if isfield(cfg, 'axesHandle') && ~isempty(cfg.axesHandle) && ishandle(cfg.axesHandle)
        ax = cfg.axesHandle;
        axes(ax); cla(ax); hold(ax, 'on');
    else
        figure; ax = gca; hold(ax, 'on');
    end

    % 先画A/B的SE带
    hFA = fill(ax, [xWin, fliplr(xWin)], [mA+seA, fliplr(mA-seA)], colors(1,:), ...
         'FaceAlpha', 0.2, 'EdgeColor', 'none');
    hFB = fill(ax, [xWin, fliplr(xWin)], [mB+seB, fliplr(mB-seB)], colors(2,:), ...
         'FaceAlpha', 0.2, 'EdgeColor', 'none');

    % 如基线为数组，画基线(SE带+均值线)为灰色
    hFC = []; hLC = [];
    if drawBaseline && isChanceArray
        grayC = [0.5 0.5 0.5];
        hFC = fill(ax, [xWin, fliplr(xWin)], [mC+seC, fliplr(mC-seC)], grayC, ...
            'FaceAlpha', 0.15, 'EdgeColor', 'none');
    end

    % 均值曲线
    hLA = plot(ax, xWin, mA, '-', 'Color', colors(1,:), 'LineWidth', 2);
    hLB = plot(ax, xWin, mB, '-', 'Color', colors(2,:), 'LineWidth', 2);
    if drawBaseline && isChanceArray
        hLC = plot(ax, xWin, mC, '-', 'Color', [0.5 0.5 0.5], 'LineWidth', 2);
    end

    % 标量基线画水平虚线
    if drawBaseline && ~isChanceArray
        yline(ax, chanceLevel, '--', 'Color', [0.4 0.4 0.4], 'LineWidth', 1.2, ...
            'Label', chanceLabel, 'LabelHorizontalAlignment','left');
    end

    % 显著性标记
    % 计算显著性标记的y坐标，考虑多子图时以figure高度折算
    yl = ylim(ax);
    axPos = get(ax, 'Position'); % [left bottom width height] in normalized units
    fig = ancestor(ax, 'figure');
    figPos = get(fig, 'Position'); % [left bottom width height] in pixels
    axHeightPix = axPos(4) * figPos(4); % axes高度（像素）
    yRange = yl(2) - yl(1);
    % 单位像素对应的y轴长度
    yPerPix = yRange / axHeightPix;

    % 线宽（以点为单位），1点=1/72英寸，1英寸=96像素（通常），但matlab线宽单位为点
    lineW = 3; % 线宽
    % 1点 = 1.333像素
    lineW_pix = lineW * 1.333; % 线宽对应像素
    lineW_y = lineW_pix * yPerPix; % 线宽对应y轴长度

    % 以线宽为单位下移，两条线之间下移2个线宽
    y0 = yl(1) - 1 * lineW_y;
    yA = y0 - 2 * lineW_y;
    yB = yA - 2 * lineW_y;

    % 组间
    if isfield(stat, 'prob')
        sig = stat.prob < 0.05;
        if any(sig)
            d = diff([0, sig, 0]);
            sIdx = find(d == 1); eIdx = find(d == -1) - 1;
            for kk = 1:numel(sIdx)
                plot(ax, [xWin(sIdx(kk)), xWin(eIdx(kk))], [y0, y0], 'r-', 'LineWidth', lineW);
            end
        end
    end
    % A vs 基线
    if isfield(statA, 'prob')
        sigA = statA.prob < 0.05;
        if any(sigA)
            dA = diff([0, sigA, 0]);
            sA = find(dA == 1); eA = find(dA == -1) - 1;
            for kk = 1:numel(sA)
                plot(ax, [xWin(sA(kk)), xWin(eA(kk))], [yA, yA], '-', 'Color', colors(1,:), 'LineWidth', lineW);
            end
        end
    end
    % B vs 基线
    if isfield(statB, 'prob')
        sigB = statB.prob < 0.05;
        if any(sigB)
            dB = diff([0, sigB, 0]);
            sB = find(dB == 1); eB = find(dB == -1) - 1;
            for kk = 1:numel(sB)
                plot(ax, [xWin(sB(kk)), xWin(eB(kk))], [yB, yB], '-', 'Color', colors(2,:), 'LineWidth', lineW);
            end
        end
    end

    ylabel(ax, cfg.ylabel);
    xlabel(ax, cfg.xlabel);

    % 图例
    if drawBaseline && isChanceArray
        legendCell = { ...
            [labels{1},'±SE'], [labels{2},'±SE'], [chanceLabel,'±SE'], ...
            labels{1}, labels{2}, chanceLabel};
    else
        legendCell = {[labels{1},'±SE'], [labels{2},'±SE'], labels{1}, labels{2}};
        if drawBaseline
            legendCell{end+1} = chanceLabel;
        end
    end
    legend(ax, legendCell, 'Location', 'northeastoutside');
    set(ax, 'FontSize', 12, 'LineWidth', 1.2); box off; hold(ax, 'off');
    % 自动调整y轴范围，考虑显著性标记的yB
    if drawBaseline && isChanceArray
        y_min = min([min(mA-seA), min(mB-seB), min(mC-seC), yB]);
        y_max = max([max(mA+seA), max(mB+seB), max(mC+seC), y0]);
    else
        y_min = min([chanceLevel, min(mA-seA), min(mB-seB), yB]);
        y_max = max([chanceLevel, max(mA+seA), max(mB+seB), y0]);
    end
    y_pad = 0.08 * (y_max - y_min + eps);
    ylim(ax, [y_min - y_pad, y_max + y_pad]);
end
