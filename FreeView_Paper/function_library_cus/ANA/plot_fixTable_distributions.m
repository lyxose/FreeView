function plot_fixTable_distributions(visTable, img_width, img_height, heat_binSize, ut, center, ver)
% plot_fixTable_distributions 可视化注视点表格的主要分布特征
% 输入:
%   visTable   - 注视点表格（已筛选/可用），需包含 dropFix, saccLen, saccAngle, r, theta, dur, xpos, ypos, tgr, tgtheta 等字段
%   img_width  - 图像宽度（像素）
%   img_height - 图像高度（像素）
%   heat_binSize    - 热图网格尺寸（像素）
%   ut         - 工具类，需支持 Pol2Rect
%   center     - 图像中心坐标 [x, y]

    % 计算任务空间的方形范围（以rmax为半径，中心为center）
    rmax = 7.5;
    if nargin < 6 || isempty(center)
        center = [img_width, img_height] / 2;
    end
    if nargin < 5 || isempty(ut)
        error('ut (工具类) 必须提供');
    end
    rmax_pix = ut.deg2pix(rmax);
    xlim_box = round(center(1) + [-1, 1] * rmax_pix);
    ylim_box = round(center(2) + [-1, 1] * rmax_pix);

    figure('Units','normalized','Position',[0 0 1 1]);
    subplot(2,3,1);
    histogram(visTable.saccLen(~visTable.saccLen==0), 0:0.1:20);
    xlabel('Saccade Length (deg)'); ylabel('Count');
    xlim([0 20]);
    title(sprintf('Saccade Length Distribution (N=%d)', height(visTable)));

    subplot(2,3,2);
    histogram(visTable.saccAngle(~visTable.saccLen==0), -180:5:180);
    xlabel('Saccade Angle (deg)'); ylabel('Count');
    xlim([-180 180]);
    title(sprintf('Saccade Angle Distribution (N=%d)', height(visTable)));

    subplot(2,3,3);
    histogram(visTable.r, 0:0.1:20);
    xlabel('Eccentricity (deg)'); ylabel('Count');
    xlim([0 20]);
    title(sprintf('Eccentricity Distribution (N=%d)', height(visTable)));

    subplot(2,3,4);
    histogram(mod(visTable.theta,360), 0:5:360);
    xlabel('Fixation Angle (deg)'); ylabel('Count');
    xlim([0 360]);
    title(sprintf('Fixation Angle Distribution (N=%d)', height(visTable)));

    subplot(2,3,5);
    histogram(visTable.dur, 0:10:1000);
    xlabel('Fixation Duration (ms)'); ylabel('Count');
    xlim([0 1000]);
    title(sprintf('Fixation Duration Distribution (N=%d)', height(visTable)));

    subplot(2,3,6);
    % 只绘制任务空间区域的热图
    nx = max(10, round((xlim_box(2)-xlim_box(1)) / heat_binSize));
    ny = max(10, round((ylim_box(2)-ylim_box(1)) / heat_binSize));
    edges_x = linspace(xlim_box(1), xlim_box(2), nx+1);
    edges_y = linspace(ylim_box(1), ylim_box(2), ny+1);
    mask_in_box = visTable.xpos >= xlim_box(1) & visTable.xpos <= xlim_box(2) & ...
                  visTable.ypos >= ylim_box(1) & visTable.ypos <= ylim_box(2);
    counts = histcounts2(visTable.xpos(mask_in_box), visTable.ypos(mask_in_box), edges_x, edges_y);
    density = counts';
    total = sum(density(:));
    if total > 0
        density = density / total;
    end
    centers_x = edges_x(1:end-1) + diff(edges_x)/2;
    centers_y = edges_y(1:end-1) + diff(edges_y)/2;
    imagesc(centers_x, centers_y, density);
    axis image;
    xlabel('X Position (pixel)'); ylabel('Y Position (pixel)');
    title(sprintf('Fixation Position Heatmap (N=%d)', sum(mask_in_box)));
    colormap("cool"); % parula
    cmax = prctile(density(:), 100);
    if cmax <= 0
        cmax = max(density(:));
    end
    if isempty(cmax) || cmax == 0
        clim auto;
    else
        clim([0, cmax]);
    end
    colorbar;

    if ~strcmpi(ver, 'v2')
        hold on;
        uTg = center + ut.Pol2Rect(unique(round([visTable.tgr, visTable.tgtheta]), 'rows', 'stable'));
        % 只显示在任务空间内的目标点
        mask_tg = uTg(:,1) >= xlim_box(1) & uTg(:,1) <= xlim_box(2) & ...
                  uTg(:,2) >= ylim_box(1) & uTg(:,2) <= ylim_box(2);
        scatter(uTg(mask_tg,1), uTg(mask_tg,2), 10, 'r', ...
            'MarkerEdgeColor','r','MarkerFaceAlpha',0.8);
        hold off;
    end
end
