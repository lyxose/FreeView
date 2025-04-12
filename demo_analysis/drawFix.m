% draw fixations
% 假设dat.fix.xpos和dat.fix.ypos是1x51的double类型数据
% 设定图像的尺寸和分辨率
img_width = 1920;  % 二维空间的宽度
img_height = 1080; % 二维空间的高度

% 假设dat.fix.xpos和dat.fix.ypos是1x51的double类型数据
xpos = dat.fix.xpos;
ypos = dat.fix.ypos;

% 1. 对二维平面中的点进行核密度估计（KDE），以估算每个位置的点出现概率
% 创建一个网格来表示二维空间
edges_x = linspace(1, img_width, 50);  % 网格的X边界
edges_y = linspace(1, img_height, 50); % 网格的Y边界

% 2. 使用二维直方图，计算数据在网格上的分布
counts = hist3([xpos', ypos'], 'Edges', {edges_x, edges_y});

% 3. 将结果转换为二维密度
density_map = counts' / sum(counts(:));  % 将计数转换为密度

% 4. 绘制热图
figure;
imagesc(edges_x, edges_y, density_map);  % 使用imagesc绘制热图
axis xy; % 使y轴方向与图像坐标系一致
colorbar; % 显示颜色条
title('Density Heatmap');
xlabel('X Position');
ylabel('Y Position');
