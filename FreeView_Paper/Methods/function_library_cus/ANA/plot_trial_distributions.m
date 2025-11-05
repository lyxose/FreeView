function plot_trial_distributions(visTable)
% plot_trial_distributions 绘制每trial的fixation数量分布（左）与trial时间分布（右）
% 输入:
%   visTable - fixation table，需包含 TriID, Nfix, startT, dur 字段

    % 统计每trial的fixation数量分布
    fixNumPerTrial = groupsummary(visTable, 'TriID', 'max', 'Nfix');
    triFixNp95 = prctile(fixNumPerTrial.max_Nfix, 95);

    % 统计每trial的时间分布
    G = findgroups(visTable.TriID);
    trialTimePerTrial = splitapply(@(t) max(t), visTable.startT + visTable.dur, G);
    triTime95 = prctile(trialTimePerTrial, 95);

    % 合并为一张大图的左右子图
    figure;
    % 上子图：fixation number per trial
    subplot(2,1,1);
    histogram(fixNumPerTrial.max_Nfix, 0.5:1:50.5);
    hold on;
    yl = ylim;
    plot([triFixNp95, triFixNp95], yl, 'r--', 'LineWidth', 1.5, 'DisplayName', '95th Percentile');
    hold off;
    xlabel('Number of Fixations per Trial'); ylabel('Count');
    xlim([0 51]);
    title(sprintf('Fixations per Trial (N=%d trials)', height(fixNumPerTrial)));
    legend('Location','northeast');

    % 下子图：trial time distribution (ms)
    subplot(2,1,2);
    maxT = ceil(max(trialTimePerTrial(~isnan(trialTimePerTrial))));
    binEdges = 0:50:max(triTime95+3000,50);
    histogram(trialTimePerTrial, 'BinEdges', binEdges);
    hold on;
    yl2 = ylim;
    plot([triTime95, triTime95], yl2, 'r--', 'LineWidth', 1.5, 'DisplayName', '95th Percentile');
    hold off;
    xlabel('Trial Time (ms)'); ylabel('Count');
    title(sprintf('Trial Time Distribution (N=%d trials)', numel(trialTimePerTrial)));
    legend('Location','northeast');

    % 总标题
    if exist('sgtitle','file') == 2
        sgtitle('Trial-level Distributions: Fixation Count and Trial Time');
    end
end
