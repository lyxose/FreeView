function plot_sector_fisherz_clusters(stat_axis, stat_card, time_ms, fz_axis_obs, fz_card_obs, perm_z8, perm_z4, AxisColor, GapColor, ObliColor)
% PLOT_SECTOR_FISHERZ_CLUSTERS 绘制扇区相关性 Fisher z 时程及显著性簇
% stat_axis, stat_card: 结构体，包含簇置换检验结果
% time_ms: 时间轴（ms）
% fz_axis_obs, fz_card_obs: 观察到的 Fisher z 序列
% perm_z8, perm_z4: 置换得到的 Fisher z 序列
    axisBase = GapColor;
    cardBase = ObliColor;

    figure;
    % Axis subplot
    subplot(2,1,1); hold on;
    q = prctile(perm_z8, [0 95], 1);
    q_lo = q(1,:); q_hi = q(2,:);
    hFillA = fill([time_ms, fliplr(time_ms)], [q_hi, fliplr(q_lo)], axisBase, ...
        'EdgeColor','none', 'FaceAlpha',0.18);
    hMeanA = plot(time_ms, stat_axis.mu_perm, '-', 'Color', axisBase*0.9, 'LineWidth',1.4);
    hObsA = plot(time_ms, fz_axis_obs, '-', 'Color', AxisColor, 'LineWidth',1.8);
    yl = ylim; y0 = yl(1) + 0.02*(yl(2)-yl(1));
    hCluA = []; hasCluA = false;
    for ii = 1:numel(stat_axis.clu_start)
        if stat_axis.clu_p(ii) < 0.05
            xs = time_ms(stat_axis.clu_start(ii));
            xe = time_ms(stat_axis.clu_end(ii));
            h = plot([xs, xe], [y0, y0], '-', 'Color', AxisColor, 'LineWidth', 6);
            if ~hasCluA, hCluA = h; hasCluA = true; end
        end
    end
    xlim([0, time_ms(end)]);
    xlabel('Time (ms)'); ylabel('Fisher z (mean r)');
    title('Axis sector correlation & clusters');
    if hasCluA
        legend([hFillA, hMeanA, hObsA, hCluA], {'0–95% perm range','Perm mean','Observed Fisher z (Axis)','Signif. cluster (p<0.05)'}, 'Location','northeastoutside');
    else
        legend([hFillA, hMeanA, hObsA], {'0–95% perm range','Perm mean','Observed Fisher z (Axis)'}, 'Location','northeastoutside');
    end
    box off; hold off;

    % Card subplot
    subplot(2,1,2); hold on;
    q = prctile(perm_z4, [0 95], 1);
    q_lo = q(1,:); q_hi = q(2,:);
    hFillC = fill([time_ms, fliplr(time_ms)], [q_hi, fliplr(q_lo)], cardBase, ...
        'EdgeColor','none', 'FaceAlpha',0.18);
    hMeanC = plot(time_ms, stat_card.mu_perm, '-', 'Color', cardBase*0.9, 'LineWidth',1.4);
    hObsC  = plot(time_ms, fz_card_obs, '-', 'Color', AxisColor, 'LineWidth',1.8);
    yl = ylim; y0 = yl(1) + 0.02*(yl(2)-yl(1));
    hCluC = []; hasCluC = false;
    for ii = 1:numel(stat_card.clu_start)
        if stat_card.clu_p(ii) < 0.05
            xs = time_ms(stat_card.clu_start(ii));
            xe = time_ms(stat_card.clu_end(ii));
            h = plot([xs, xe], [y0, y0], '-', 'Color', AxisColor, 'LineWidth', 6);
            if ~hasCluC, hCluC = h; hasCluC = true; end
        end
    end
    xlim([0, time_ms(end)]);
    xlabel('Time (ms)'); ylabel('Fisher z (mean r)');
    title('Card sector correlation & clusters');
    if hasCluC
        legend([hFillC, hMeanC, hObsC, hCluC], {'0–95% perm range','Perm mean','Observed Fisher z (Card)','Signif. cluster (p<0.05)'}, 'Location','northeastoutside');
    else
        legend([hFillC, hMeanC, hObsC], {'0–95% perm range','Perm mean','Observed Fisher z (Card)'}, 'Location','northeastoutside');
    end
    box off; hold off;

    % Print significant clusters to console
    disp('Axis clusters (p<0.05): [start(ms) end(ms) p]');
    for ii = 1:numel(stat_axis.clu_start)
        if stat_axis.clu_p(ii) < 0.05
            fprintf('[%g  %g]  p=%.4f\n', time_ms(stat_axis.clu_start(ii)), time_ms(stat_axis.clu_end(ii)), stat_axis.clu_p(ii));
        end
    end
    disp('Card clusters (p<0.05): [start(ms) end(ms) p]');
    for ii = 1:numel(stat_card.clu_start)
        if stat_card.clu_p(ii) < 0.05
            fprintf('[%g  %g]  p=%.4f\n', time_ms(stat_card.clu_start(ii)), time_ms(stat_card.clu_end(ii)), stat_card.clu_p(ii));
        end
    end
end