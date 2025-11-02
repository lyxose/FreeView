function plot_sector_stats_time_series(sub_time_bin_FT, timeCenters_FT, nsbj)
% 绘制扇区组统计（range/std）随时间的均值±SE曲线
% 输入：
%  - sub_time_bin_FT: [nsbj x nTime x nSector]
%  - timeCenters_FT: 时间轴（ms），长度nTime
%  - nsbj: 被试数
%  - nTime: 时间点数
%  - nSector: 扇区数（16）

% 计算每被试每时间点的range和std（扇区维度）
range_series = squeeze(max(sub_time_bin_FT,[],3) - min(sub_time_bin_FT,[],3)); % [nsbj x nTime]
std_series   = squeeze(std(sub_time_bin_FT,0,3)); % [nsbj x nTime]

% 组均值±SE（被试维度）
mRange = mean(range_series, 1);
seRange = std(range_series, 0, 1) / sqrt(max(1, nsbj));
mStd = mean(std_series, 1);
seStd = std(std_series, 0, 1) / sqrt(max(1, nsbj));

figure;
hold on;
% % Range 阴影与曲线
fill([timeCenters_FT, fliplr(timeCenters_FT)], [mRange+seRange, fliplr(mRange-seRange)], [0.8,0.8,1], ...
     'FaceAlpha',0.5,'EdgeColor','none');
plot(timeCenters_FT, mRange, 'b-', 'LineWidth',2);
% Std 阴影与曲线
fill([timeCenters_FT, fliplr(timeCenters_FT)], [mStd+seStd, fliplr(mStd-seStd)], [1,0.8,0.8], ...
     'FaceAlpha',0.5,'EdgeColor','none');
plot(timeCenters_FT, mStd, 'r-', 'LineWidth',2);

xlabel('Time (ms)');
ylabel('Value');
% legend({'Std ± SE','Std Mean'}, 'Location', 'best');
legend({'Range ± SE','Range Mean','Std ± SE','Std Mean'}, 'Location', 'best');
title('Range and Std of Detrended Sector Counts over Time');
set(gca,'FontSize',12,'LineWidth',1.2); box off; hold off;
end

