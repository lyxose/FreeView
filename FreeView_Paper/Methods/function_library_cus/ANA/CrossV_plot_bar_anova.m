function CrossV_plot_bar_anova(data, group, labels, ylabel_str, chanceLevel)
    % 绘制分组柱状图（占比数据），以chance level为基线，显著性与chance level比较
    % data: nSubj x 1
    % group: nSubj x 1, 分组编号（1,2,...）
    % labels: cell array, 分组标签
    % ylabel_str: y轴标题
    % chanceLevel: (可选) 基线，默认0.5
    % 多重检验校正：
    % 组内：fdr_bh
    % 组间：tukey-kramer

    if nargin < 5 || isempty(chanceLevel)
        chanceLevel = 0.5;
    end

    nGroups = numel(labels);
    group = group(:);
    data = data(:);

    m = zeros(1, nGroups);
    se = zeros(1, nGroups);
    indYmax = -inf(1, nGroups);
    for i = 1:nGroups
        vals = data(group==i);
        if isempty(vals)
            m(i) = NaN; se(i) = NaN;
        else
            m(i) = mean(vals, 'omitnan');
            se(i) = std(vals, 0, 'omitnan') / sqrt(max(1, numel(vals)));
            indYmax(i) = max(vals);
        end
    end
    indYmax(~isfinite(indYmax)) = m(~isfinite(indYmax)) + se(~isfinite(indYmax));

    if nGroups > 2
        [p_interact, tbl, stats] = anova1(data, group, 'off');
        mc = multcompare(stats, 'Display', 'off'); % 默认tukey-kramer校正

        fprintf('%s: One-way ANOVA across groups: p = %.4g\n', ylabel_str, p_interact);
        disp('Multiple comparisons between groups:');
        disp(array2table(mc, 'VariableNames', {'Group1','Group2','LowerCI','MeanDiff','UpperCI','pValue'}));
    else
        % 两组时只做t检验
        mc = [1 2 NaN mean(m(1)-m(2)) NaN ranksum(data(group==1), data(group==2))];
    end

    figure;
    hold on;
    x = 1:nGroups;
    cmap = lines(nGroups);

    % 绘制chance level基线
    yline(chanceLevel, '--', 'Color', [0.4 0.4 0.4], 'LineWidth', 1.2, 'Label', 'Gap', 'LabelHorizontalAlignment','left');

    % 柱子从chance level起始
    for i = 1:nGroups
        b = bar(x(i), m(i), 0.7, 'FaceColor', cmap(i,:), 'BaseValue', chanceLevel, 'EdgeColor','none');
        b.FaceAlpha = 0.5;
    end

    % errorbar
    errorbar(x, m, se, 'LineStyle','none','Color',[0.2 0.2 0.2],'LineWidth',1.2);

    % individual 数据点
    bar_alpha_indv = 0.5;
    if nGroups > 1
        jitter_amp = 0.18;
    else
        jitter_amp = 0;
    end
    indYmax = nan(1,nGroups);
    for i = 1:nGroups
        vals = data(group==i);
        if isempty(vals), continue; end
        xj = x(i) + (rand(numel(vals),1)-0.5)*jitter_amp;
        scatter(xj, vals, 32, cmap(i,:), 'filled', 'MarkerFaceAlpha',0.7, 'MarkerEdgeAlpha',0.5);
        indYmax(i) = max(vals);
    end
    y_span = ylim;
    bump = (y_span(2) - y_span(1)) * 0.05;
    if bump <= 0
        bump = 0.05 * max(abs(data))+1;
    end
    % 单样本 t 检验（对 chanceLevel 做检验），并标注显著性符号
    pvals = nan(1, nGroups);
    for i = 1:nGroups
        vals = data(group==i);
        try
            [~, pvals(i), ~, stats] = ttest(vals, chanceLevel);
        catch
            pvals(i) = NaN;
            stats.tstat = NaN;
        end
    end

    % FDR 校正
    [~, ~, ~, pvals_fdr] = fdr_bh(pvals);

    for i = 1:nGroups
        s = sig_symbol(pvals_fdr(i));
        if m(i) >= chanceLevel
            y = m(i) + se(i) + bump;
            vAlign = 'bottom';
        else
            y = m(i) - se(i) - bump;
            vAlign = 'top';
        end
        text(x(i), y, s, 'HorizontalAlign','center','VerticalAlign',vAlign,...
            'FontSize',12,'FontWeight','bold');
    end
    indYmax(isnan(indYmax)) = m(isnan(indYmax)) + se(isnan(indYmax));
    indYmax = indYmax + bump;
    topUsed = indYmax;

    pairOrder = [1 2; 2 3; 1 3];
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
        text(mean([x1,x2]), y_bar + bump*0.6, sprintf('p=%.3g', pval), ...
            'HorizontalAlignment','center','VerticalAlign','bottom', ...
            'FontSize',11,'FontWeight','bold','Color', lineColor);
        topUsed(gPair(1)) = max(topUsed(gPair(1)), y_bar + bump*1.4);
        topUsed(gPair(2)) = max(topUsed(gPair(2)), y_bar + bump*1.4);
        if nGroups == 2
            break; % 只画一次
        end
    end

    yl_curr = ylim;
    newYmax = max([yl_curr(2), max(topUsed)+bump]);
    ylim([max(0, yl_curr(1)), min(1, newYmax)]);

    xticks(1:nGroups);
    xticklabels(labels);
    ylabel(ylabel_str);
    title(sprintf('%s across Groups', ylabel_str));
    set(gca,'FontSize',12,'LineWidth',1.2); box off;

    if nGroups > 2
        yl_final = ylim;
        xl = xlim;
        anova_y = max([topUsed, yl_final(2) - 0.04*(yl_final(2)-yl_final(1))]);
        if anova_y >= yl_final(2)
            ylim([yl_final(1), anova_y + bump]);
            yl_final = ylim;
            anova_y = yl_final(2) - 0.04*(yl_final(2)-yl_final(1));
        end
        color_anova = ternary(p_interact < 0.05, [1 0 0], 'k');
        text(xl(1)+0.02*(xl(2)-xl(1)), anova_y, sprintf('ANOVA p=%.3g', p_interact), ...
            'FontSize',12,'FontWeight','bold', 'HorizontalAlignment','left', 'VerticalAlignment','top','Color',color_anova);
    end
    hold off;
end
