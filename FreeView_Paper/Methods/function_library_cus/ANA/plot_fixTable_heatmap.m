function plot_fixTable_heatmap(xpos_FT, ypos_FT, heat_binSize, img_width, img_height)
% plot_fixTable_heatmap 绘制基于注视点坐标的二维热图（密度图），y轴翻转使0在上
%
% 输入参数:
%   xpos_FT    - 注视点的 X 坐标向量（单位: pixel）
%   ypos_FT    - 注视点的 Y 坐标向量（单位: pixel）
%   heat_binSize    - 热图网格的像素大小（如 25）
%   img_width  - （可选）画布宽度（像素），用于设定边界
%   img_height - （可选）画布高度（像素），用于设定边界
%
% 功能说明:
%   - 自动根据数据范围或指定画布尺寸生成网格边界
%   - 统计每个网格内的注视点数量，归一化为密度
%   - 使用 imagesc 绘制热图，坐标比例一致
%   - y轴翻转，使0在上方
%
% 示例用法:
%   plot_fixTable_heatmap(xpos_FT, ypos_FT, 25, 1920, 1080);

    % 1. 设定画布边界
    if nargin >= 4 && ~isempty(img_width) && ~isempty(img_height)
        ex = [1, img_width];
        ey = [1, img_height];
    else
        ex = [floor(min(xpos_FT)), ceil(max(xpos_FT))];
        ey = [floor(min(ypos_FT)), ceil(max(ypos_FT))];
    end

    % 2. 生成网格边界
    nx = max(10, round((ex(2)-ex(1))/heat_binSize));
    ny = max(10, round((ey(2)-ey(1))/heat_binSize));
    edges_x = linspace(ex(1), ex(2), nx+1);
    edges_y = linspace(ey(1), ey(2), ny+1);

    % 3. 统计二维注视点分布
    [count2D,~,~] = histcounts2(xpos_FT, ypos_FT, edges_x, edges_y);

    % 4. 归一化为密度（总和为1）
    density = count2D' / max(1, sum(count2D(:)));

    % 5. 绘制热图
    figure;
    imagesc(edges_x, edges_y, density);
    axis image; colorbar;
    set(gca,'YDir','reverse'); % 翻转y轴，使0在上方
    colormap("cool");

    % 6. 设置颜色范围（自动或按最大密度）
    cmax = prctile(density(:), 100);
    if isempty(cmax) || cmax == 0
        clim auto;
    else
        caxis([min(density(:)), cmax]);
    end

    xlabel('X Position (pixel)');
    ylabel('Y Position (pixel)');
end
