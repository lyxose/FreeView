function plot_fixation_sectors(xpos, ypos, video_width, video_height, R_max, cmap16, sector_edges, ut, ax)
% 在指定的axes上绘制所有注视点在扇区背景上的分布
% 以上方为90°，右侧为0°，即极角0°在右侧，90°在上方
% 如果未传入ax,则自动新建figure和axes

if nargin < 9 || isempty(ax) || ~isvalid(ax)
    fig = figure('Name', 'Fixation Sectors', 'Color', 'w');
    ax = axes('Parent', fig);
end

% 计算中心
center_pix = [video_width, video_height] / 2;

% 计算扇区半径: 使用最外侧注视点的距离
distFromCenter = sqrt((xpos - center_pix(1)).^2 + (ypos - center_pix(2)).^2);
RpatchPix = max(distFromCenter) * 1.01; % 略大于最远注视点
theta_patch_res = 80;

% 生成浅色背景
light_cmap = min(1, 0.15 + 0.85 * cmap16);

% 绘制扇区背景
sectorPolys = cell(16, 1);
labelPos = zeros(16, 2);

for si = 1:16
    a1 = sector_edges(si);
    a2 = sector_edges(si + 1);
    if a2 < a1
        a2 = a2 + 360;
    end
    aa = linspace(a1, a2, theta_patch_res);
    % 以上方为90°，右侧为0°
    xp = center_pix(1) + RpatchPix * cosd(aa);
    yp = center_pix(2) - RpatchPix * sind(aa);
    sectorPolys{si} = [center_pix; [xp(:), yp(:)]];
    
    % 标签位置和文本(角度值)
    midA = mod((a1 + (a2 - a1) / 2), 360);
    labR = RpatchPix * 1.12;
    labelPos(si, :) = center_pix + [labR * cosd(midA), -labR * sind(midA)];
end

% 设置绘图范围为正方形,包含圆形区域和外围标签
margin = RpatchPix * 0.25;
plot_center = center_pix;
plot_radius = RpatchPix * 1.2;
xlim(ax, plot_center(1) + [-plot_radius, plot_radius]);
ylim(ax, plot_center(2) + [-plot_radius, plot_radius]);

% 设置axes属性: 保持纵横比一致,隐藏坐标轴
set(ax, 'YDir', 'reverse', 'DataAspectRatio', [1 1 1], ...
    'XTick', [], 'YTick', [], 'Box', 'off', 'Visible', 'off');
hold(ax, 'on');

% 绘制扇区
for si = 1:16
    P = sectorPolys{si};
    patch('XData', P(:, 1), 'YData', P(:, 2), ...
          'FaceColor', light_cmap(si, :), ...
          'EdgeColor', 'none', 'FaceAlpha', 0.28, 'Parent', ax);
end

% 绘制扇区标签(角度, 整数不带小数, 非整数保留一位小数)
for si = 1:16
    midA = mod((sector_edges(si) + sector_edges(si + 1)) / 2, 360);
    if abs(midA - round(midA)) < 1e-6
        lbl = sprintf('%d°', round(midA));
    else
        lbl = sprintf('%.1f°', midA);
    end
    text(labelPos(si, 1), labelPos(si, 2), lbl, ...
         'Parent', ax, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
         'FontWeight', 'bold', 'FontSize', 9, 'Color', [0.06 0.06 0.06]);
end

% 计算每个注视点的扇区归属和颜色
% 以上方为90°，右侧为0°
angA = mod(atan2d(center_pix(2) - ypos, xpos - center_pix(1)), 360);
sectorIdx = zeros(numel(angA), 1, 'uint8');

for si = 1:16
    l = mod(sector_edges(si), 360);
    r = mod(sector_edges(si + 1), 360);
    if l < r
        mask = angA >= l & angA < r;
    else
        mask = angA >= l | angA < r;
    end
    sectorIdx(mask) = si;
end

% 为注视点分配颜色
ptColors = cmap16(double(sectorIdx), :);

% 绘制注视点(缩小尺寸)
scatter(ax, xpos, ypos, 2, ptColors, 'filled', 'MarkerFaceAlpha', 0.85, 'MarkerEdgeAlpha', 0.4);

hold(ax, 'off');

end
