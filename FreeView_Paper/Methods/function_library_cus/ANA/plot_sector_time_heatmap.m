function plot_sector_time_heatmap(timeCenters_FT, nSector, sub_time_bin_FT, edges_FT, ylab)
% 绘制按时间的16扇区组均值热图
% 输入：
%  - timeCenters_FT: 时间轴（ms），长度nTime
%  - nsbj: 被试数
%  - nTime: 时间点数
%  - nSector: 扇区数（16）
%  - sub_time_bin_FT: [nsbj x nTime x nSector]，去趋势计数
%  - edges_FT: 16-bin边界（用于扇区标签）
%  - ylab: y轴标签（如 'Normalized Count'）

% 计算组均值（被试平均），结果为[nTime x nSector]
mean_over_sub = squeeze(mean(sub_time_bin_FT, 1)); % [nTime x nSector]

figure;
imagesc(timeCenters_FT, 1:nSector, mean_over_sub'); % 转置后: 16 x nTime
set(gca,'YDir','normal');
xlabel('Time (ms)'); ylabel('Sector (1-16)');
title('Mean Detrended Sector Count over Time');

% 扇区中心角度标签
if ~isempty(edges_FT) && numel(edges_FT) >= 17
    bin_centers = mod((edges_FT(1:end-1) + edges_FT(2:end))/2, 360);
else
    bin_centers = mod(linspace(0, 360, nSector+1), 360);
    bin_centers = bin_centers(1:end-1);
end
yticks(1:nSector);
yticklabels(arrayfun(@(x) sprintf('%.1f°', x), bin_centers, 'UniformOutput', false));
set(gca,'FontSize',12,'LineWidth',1.2);
colormap("cool");
cb = colorbar;
ylabel(cb, ylab);
end
