function CrossV_plot_violin_anova(data, group, labels, ylabel_str, chanceLevel)
    % 多组小提琴图 + 单样本t检验 + 组间ANOVA/多重比较
    % data: N x 1
    % group: N x 1, 组别(1,2,...)
    % labels: cell array, 每组标签
    % ylabel_str: y轴标签
    % chanceLevel: 参考线, 默认0.5

    if nargin < 5 || isempty(chanceLevel)
        chanceLevel = 0.5;
    end

    nGroups = numel(labels);
    group = group(:);
    data = data(:);

    % --- 3sd异常值筛选 ---
    outlier_mask = false(size(data));
    for i = 1:nGroups
        vals = data(group==i);
        mu = mean(vals,'omitnan');
        sd = std(vals,'omitnan');
        mask = abs(vals-mu) > 3*sd;
        idx = find(group==i);
        outlier_mask(idx(mask)) = true;
    end

    m = zeros(1, nGroups);
    se = zeros(1, nGroups);
    for i = 1:nGroups
        mask_in = group==i & ~outlier_mask;
        vals = data(mask_in);
        if isempty(vals)
            m(i) = NaN; se(i) = NaN;
        else
            m(i) = mean(vals, 'omitnan');
            se(i) = std(vals, 0, 'omitnan') / sqrt(max(1, numel(vals)));
        end
    end

    % 组间统计
    if nGroups > 2
        [p_interact, tbl, stats] = anova1(data(~outlier_mask), group(~outlier_mask), 'off');
        mc = multcompare(stats, 'Display', 'off');
        fprintf('%s: One-way ANOVA across groups: p = %.4g\n', ylabel_str, p_interact);
        disp('Multiple comparisons between groups:');
        disp(array2table(mc, 'VariableNames', {'Group1','Group2','LowerCI','MeanDiff','UpperCI','pValue'}));
    else
        mask1 = group==1 & ~outlier_mask;
        mask2 = group==2 & ~outlier_mask;
        mc = [1 2 NaN mean(m(1)-m(2)) NaN ranksum(data(mask1), data(mask2))];
    end

    figure;
    hold on;
    x = 1:nGroups;
    cmap = lines(nGroups);

    % 绘制chance level基线
    yline(chanceLevel, '--', 'Color', [0.4 0.4 0.4], 'LineWidth', 1.2, 'Label', 'Chance', 'LabelHorizontalAlignment','left');

    violin_width = 0.3;
    xi_all = [];
    bump = [];
    indYmax = nan(1, nGroups);

    for i = 1:nGroups
        vals = data(group==i & ~outlier_mask);
        vals = vals(~isnan(vals));
        if isempty(vals), continue; end

        % 核密度估计
        [f, xi] = ksdensity(vals, 'BoundaryCorrection','reflection');
        f = f / max(f) * violin_width;
        xi_all = [xi_all, xi(:)];

        % 小提琴
        patch([x(i) - f, x(i) + fliplr(f)], [xi, fliplr(xi)], cmap(i,:), ...
            'FaceAlpha', 0.3, 'EdgeColor', cmap(i,:), 'LineWidth', 1);

        % 箱线图元素
        q = prctile(vals, [25 50 75]);
        iqr = q(3) - q(1);
        whisker_low = max(min(vals), q(1) - 1.5*iqr);
        whisker_high = min(max(vals), q(3) + 1.5*iqr);

        % 须
        plot([x(i) x(i)], [whisker_low q(1)], 'Color', cmap(i,:)*0.6, 'LineWidth', 1.2);
        plot([x(i) x(i)], [q(3) whisker_high], 'Color', cmap(i,:)*0.6, 'LineWidth', 1.2);

        % 箱体
        box_w = 0.08;
        rectangle('Position', [x(i)-box_w/2, q(1), box_w, q(3)-q(1)], ...
            'FaceColor', [1 1 1 0.7], 'EdgeColor', cmap(i,:)*0.6, 'LineWidth', 1.5);

        % 中位数线
        plot([x(i)-box_w/2 x(i)+box_w/2], [q(2) q(2)], ...
            'Color', cmap(i,:)*0.4, 'LineWidth', 2);

        indYmax(i) = max(vals);
    end

    % individual 数据点
    jitter_amp = 0.15;
    for i = 1:nGroups
        idx = find(group==i);
        vals = data(idx);
        mask_in = ~outlier_mask(idx) & ~isnan(vals);
        mask_out = outlier_mask(idx) & ~isnan(vals);
        xj = x(i) + (rand(numel(vals),1)-0.5)*jitter_amp;
        scatter(xj(mask_in), vals(mask_in), 24, cmap(i,:), 'filled', 'MarkerFaceAlpha',0.7, 'MarkerEdgeAlpha',0.5);
        scatter(xj(mask_out), vals(mask_out), 40, [0.5 0.5 0.5], 'x', 'LineWidth',1.8,'MarkerEdgeAlpha',0.7);
    end

    % errorbar
    % errorbar(x, m, se, 'LineStyle','none','Color',[0.2 0.2 0.2],'LineWidth',1.2);

    % 单样本 t 检验（对 chanceLevel 做检验），只用非离群点
    pvals = nan(1, nGroups);
    for i = 1:nGroups
        mask_in = group==i & ~outlier_mask;
        vals = data(mask_in);
        try
            [~, pvals(i)] = ttest(vals, chanceLevel);
        catch
            pvals(i) = NaN;
        end
    end

    [~, ~, ~, pvals_fdr] = fdr_bh(pvals);

    % bump for annotation
    y_span = ylim;
    bump = (y_span(2) - y_span(1)) * 0.05;
    if bump <= 0
        bump = 0.05 * max(abs(data))+1;
    end

    % 组间显著性连线
    % 框线起始位置只考虑非离群点
    for i = 1:nGroups
        s = sig_symbol(pvals_fdr(i));
        % 框线起始位置（非离群点最大值）和须的位置
        vals = data(group==i & ~outlier_mask);
        if isempty(vals)
            continue;
        end
        q = prctile(vals, [25 50 75]);
        iqr = q(3) - q(1);
        whisker_high = min(max(vals), q(3) + 1.5*iqr);

        % 星号y坐标：须的上方，框线起始位置（indYmax(i)）的下方
        y_star = whisker_high + 0.5 * (indYmax(i) - whisker_high);

        text(x(i), y_star, s, 'HorizontalAlign','center','VerticalAlign','bottom', ...
            'FontSize',12,'FontWeight','bold');
    end

    % topUsed只考虑非离群点
    indYmax(isnan(indYmax)) = m(isnan(indYmax)) + se(isnan(indYmax));
    indYmax = indYmax + bump;
    topUsed = indYmax;

    % 自动生成所有两两组合并按组间距升序排序，避免连线交叉
    pairOrder = nchoosek(1:nGroups,2);
    pairDist = abs(pairOrder(:,2) - pairOrder(:,1));
    [~, sortIdx] = sort(pairDist); % 先画相邻组
    pairOrder = pairOrder(sortIdx,:);

    for idx = 1:size(pairOrder,1)
        gPair = pairOrder(idx,:);
        if any(gPair > nGroups)
            continue;
        end
        if nGroups == 2
            pval = mc(6);
        else
            rowIdx = find((mc(:,1)==gPair(1) & mc(:,2)==gPair(2)) | (mc(:,1)==gPair(2) & mc(:,2)==gPair(1)), 1);
            if isempty(rowIdx)
                continue;
            end
            pval = mc(rowIdx,6);
        end
        x1 = x(gPair(1)); x2 = x(gPair(2));
        y1 = topUsed(gPair(1));
        y2 = topUsed(gPair(2));
        y_bar = max([y1, y2]) + bump*0.5;
        lineColor = ternary(pval < 0.05, [1 0 0], [0 0 0]);
        plot([x1, x1, x2, x2], [y1, y_bar, y_bar, y2], '-', 'Color', lineColor, 'LineWidth', 1.6);
        text(mean([x1,x2]), y_bar + bump*0.2, sprintf('p=%.4f', pval), ...
            'HorizontalAlignment','center','VerticalAlign','bottom', ...
            'FontSize',11,'FontWeight','bold','Color', lineColor);
        topUsed(gPair(1)) = max(topUsed(gPair(1)), y_bar + bump*1.4);
        topUsed(gPair(2)) = max(topUsed(gPair(2)), y_bar + bump*1.4);
        if nGroups == 2
            break; % 只画一次
        end
    end

    % ylim整体范围应包含所有数据（包括离群点）
    yl_curr = ylim;
    data_max = max(data(~isnan(data)));
    data_min = min(data(~isnan(data)));
    newYmax = max([yl_curr(2), max(topUsed)+bump, data_max+bump]);
    newYmin = min([yl_curr(1), data_min-bump]);
    ylim([newYmin, min(1, newYmax)]);

    xticks(1:nGroups);
    xticklabels(labels);
    ylabel(ylabel_str);
    title(sprintf('%s across Groups', ylabel_str));
    set(gca,'FontSize',12,'LineWidth',1.2); box off;

    if nGroups > 2
        yl_final = ylim;
        xl = xlim;
        anova_y = max([topUsed, yl_final(2) - 0.01*(yl_final(2)-yl_final(1))]);
        if anova_y >= yl_final(2)
            ylim([yl_final(1), anova_y + bump]);
            yl_final = ylim;
            anova_y = yl_final(2) - 0.01*(yl_final(2)-yl_final(1));
        end
        color_anova = ternary(p_interact < 0.05, [1 0 0], 'k');
        text(xl(1)+0.02*(xl(2)-xl(1)), anova_y, sprintf('ANOVA p=%.4f', p_interact), ...
            'FontSize',12,'FontWeight','bold', 'HorizontalAlignment','left', 'VerticalAlignment','top','Color',color_anova);
    end
    hold off;
end