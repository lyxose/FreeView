function plot_angle_curve(data, centers, foldPeriod, binSize, CardColor, GapColor, ObliColor, normMode, varargin)
% ---- 扇区滑动绘图函数 ----
% data: Nsubj × nbins 原始数据
% centers: bin中心
% foldPeriod: 360/90/45
% binSize: bin宽度
% CardColor, GapColor, ObliColor: 颜色
% normMode: 标准化模式
% varargin: 可选参数， y轴范围 [ymin ymax]、doPermu (是否进行置换检验，默认false)

    groupMean = mean(data, 1);
    groupSE = std(data, 0, 1) / sqrt(size(data, 1));
    nsbj_FT = size(data, 1);

    % 根据foldPeriod设置图形宽度比例
    if foldPeriod == 360
        figWidth = 9;  % 基准宽度
    elseif foldPeriod == 90
        figWidth = 4;  % 360的一半
    elseif foldPeriod == 45
        figWidth = 2;  % 360的四分之一
    else
        figWidth = 6;  % 默认值
    end
    figHeight = 4;  % 固定高度
    
    figure('Position', [100, 100, figWidth*100, figHeight*100]); 
    hold on;
    fill([centers, fliplr(centers)], [groupMean+groupSE, fliplr(groupMean-groupSE)], ...
         [0.7 0.7 0.7], 'EdgeColor','none','FaceAlpha',0.35);
    plot(centers, groupMean, 'Color',[0.3 0.3 0.3], 'LineWidth',1.8);

    % 处理可选参数
    doPermu = false;
    if ~isempty(varargin)
        if numel(varargin) >= 1
            yRange = varargin{1};
            if isnumeric(yRange) && numel(yRange)==2
                ylim(yRange);
            end
        end
        if numel(varargin) >= 2
            doPermu = varargin{2};
        end
    end
    yl = ylim;

    xpad = max(1, 0.01 * (max(centers) - min(centers)));
    xlim([min(centers)-xpad, max(centers)+xpad]);

    if foldPeriod == 360
        axCols = [CardColor; ObliColor];
        xtickVals = 0:45:315;
        xticks(xtickVals);
        xticklabels(arrayfun(@(x) sprintf('%.1f°',mod(x,360)), xtickVals, 'uni', 0));
        for o = 0:7
            ang = o*45;
            if ang <= max(centers) && ang >= min(centers)
                line([ang,ang], yl+[0.01,-0.01].*diff(yl), 'Color', axCols(mod(o,2)+1,:), 'LineStyle','--','LineWidth',0.8);
            end
        end
        for o = 0:7
            ang = o*45 + 22.5;
            if ang <= max(centers) && ang >= min(centers)
                line([ang,ang], yl+[0.01,-0.01].*diff(yl), 'Color', GapColor, 'LineStyle','--','LineWidth',0.8);
            end
        end
        title(sprintf('Fixation Density Scan, n=%d', nsbj_FT));
    elseif foldPeriod == 45
        xtickVals = [0, 22.5];
        xticks(xtickVals);
        xticklabels(arrayfun(@(x) sprintf('%.1f°',x), xtickVals, 'uni', 0));
        if min(centers)-0.5 <= -22.5 
            line([-22.5 -22.5], yl+[0.01,-0.01].*diff(yl), 'Color', GapColor, 'LineStyle','--','LineWidth',1.4);
        end
        line([0 0], yl+[0.01,-0.01].*diff(yl), 'Color', CardColor, 'LineStyle','--','LineWidth',1.4);
        line([22.5 22.5], yl+[0.01,-0.01].*diff(yl), 'Color', GapColor, 'LineStyle','--','LineWidth',1.4);
        % 仅针对proportion数据！
        if strcmpi(normMode, 'Proportion') && doPermu
            chance45 = 8 * binSize / 360;
            yline(chance45, 'Color', [0.5 0.5 0.5], 'LineStyle','--','LineWidth',1.4);
            % cluster-based permutation test vs. chance baseline
            plot_cluster_perm(data, centers, chance45, yl);
        end
        [~, idx0] = min(abs(centers - 0));
        [~, idx225] = min(abs(centers - 22.5));
        vals0 = data(:, idx0);
        vals225 = data(:, idx225);
        try
            [~, p_0_225] = ttest(vals0, vals225);
        catch
            p_0_225 = NaN;
        end
        diffs_0_225 = vals0 - vals225;
        d_0_225 = mean(diffs_0_225) / std(diffs_0_225, 'omitnan');
        bump = 0.05 * range(ylim);
        y0 = groupMean(idx0) + groupSE(idx0) + bump;
        y225 = groupMean(idx225) + groupSE(idx225) + bump;
        y_bar_0_225 = max([y0, y225, max(groupMean(idx0:idx225)+groupSE(idx0:idx225))]) + bump*0.5;
        color_sig_0_225 = ternary(p_0_225 < 0.05, [1 0 0], [0 0 0]);
        plot([centers(idx0), centers(idx0), centers(idx225), centers(idx225)], ...
            [y0, y_bar_0_225, y_bar_0_225, y225], '-', 'Color', color_sig_0_225, 'LineWidth', 1.8);
        text(mean([centers(idx0), centers(idx225)]), y_bar_0_225 + bump*0.6, ...
            sprintf('p=%.3f', p_0_225), ...
            'HorizontalAlignment','center', 'VerticalAlignment','bottom', ...
            'FontSize',12, 'FontWeight','bold', 'Color', color_sig_0_225);
        title(sprintf('Fixation Density Scan, n=%d', nsbj_FT));
    elseif foldPeriod == 90
        xtickVals = [0, 22.5, 45, 67.5];
        xticks(xtickVals);
        xticklabels(arrayfun(@(x) sprintf('%.1f°', x), xtickVals, 'uni', 0));
        for a = [-22.5, 0, 22.5, 45, 67.5, 90]
            if a <= max(centers)+0.5 && a >= min(centers)-0.5
            if a == 45
                lc = ObliColor;
            elseif ismember(a, [-22.5, 22.5, 67.5])
                lc = GapColor;
            else
                lc = CardColor;
            end
            line([a a], yl + [0.01, -0.01].*diff(yl), 'Color', lc, 'LineStyle', '--', 'LineWidth', 0.8);
            end
        end
        title(sprintf('Fixation Density Scan, n=%d', nsbj_FT));
        % 仅在90度条件下检验0度和45度的峰值差异显著性，并绘制显著性标记
        [~, idx0] = min(abs(centers - 0));
        [~, idx45] = min(abs(centers - 45));
        vals1 = data(:, idx0);
        vals2 = data(:, idx45);
        try
            [~, p_pair] = ttest(vals1, vals2);
        catch
            p_pair = NaN;
        end
        diffs = vals1 - vals2;
        d = mean(diffs) / std(diffs, 'omitnan');
        bump = 0.05 * range(ylim);
        y1 = groupMean(idx0) + groupSE(idx0) + bump;
        y2 = groupMean(idx45) + groupSE(idx45) + bump;
        y_bar = max([y1, y2, max(groupMean(idx0:idx45)+groupSE(idx0:idx45))]) + bump*0.5;
        color_sig = ternary(p_pair < 0.05, [1 0 0], [0 0 0]);
        plot([centers(idx0), centers(idx0), centers(idx45), centers(idx45)], ...
            [y1, y_bar, y_bar, y2], '-', 'Color', color_sig, 'LineWidth', 1.8);
        text(mean([centers(idx0), centers(idx45)]), y_bar + bump*0.6, ...
            sprintf('p=%.3f', p_pair), ...
            'HorizontalAlignment','center', 'VerticalAlignment','bottom', ...
            'FontSize',12, 'FontWeight','bold', 'Color', color_sig);
        if strcmpi(normMode, 'Proportion') && doPermu
            chance90 = 4 * binSize / 360;
            yline(chance90, 'Color', [0.5 0.5 0.5], 'LineStyle','--','LineWidth',1.4);
            % cluster-based permutation test vs. chance baseline for 90°
            plot_cluster_perm(data, centers, chance90, yl);
        end
    end
    ylim(yl);
    ylabel(['Count ', normMode]);
    xlabel(sprintf('Angle Mod %.0f°', foldPeriod));
    set(gca,'LineWidth',1.2); box off; hold off;
end

function plot_cluster_perm(data, centers, chance, yl)
    data_vs_chance = data - chance;
    [stat, ~] = cluster_based_permutation_test(data, chance * ones(size(data)), ...
    'n_perm', 5000, 'alpha', 0.05, 'tail', 'both');
    for ci = 1:numel(stat.cluster_inds)
    inds = stat.cluster_inds{ci};
    pval = stat.cluster_p(ci);
    if pval < 0.05
        idxs = inds(1):inds(2);
        mean_sign = mean(mean(data_vs_chance(:, idxs), 1), 'omitnan');
        if mean_sign > 0
        y_pos = chance - (yl(2)-yl(1))/50;
        else
        y_pos = chance + (yl(2)-yl(1))/50;
        end
        plot(centers(idxs), repmat(y_pos, size(idxs)), 'r-', 'LineWidth', 3);
    end
    end
end
