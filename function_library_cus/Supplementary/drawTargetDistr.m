% 设置图形窗口大小
figure(1);

hold on;
% axis equal;
% 设置圆的参数
circle_radius = 634.65 / 2;
center_x = 1920 / 2;
center_y = 1080 / 2;
img_width=1920;
img_height = 1080;

% 设置点的数量和直径
trialNum = 600;
diameter = 13.78;

R_max = 296.3384;
R_min = 84.2807;
% 随机生成极坐标下的点
theta = rand(trialNum, 1) * 2 * pi;  % 角度范围从0到2*pi
% r = sqrt(rand(num_points, 1)) * (296.3384-84.2807) + 84.2807;  % 半径按均匀分布生成(2~7)，sqrt保证点均匀分布
r = sqrt(R_min^2 + (R_max^2 - R_min^2) * rand(trialNum, 1)); % 调整半径分布
theta = results.Orient;
r=42.14*results.ECC;
% 将极坐标转换为笛卡尔坐标
tg_x = r .* cosd(theta) + center_x;
tg_y = r .* sind(theta) + center_y;

% x = x(1:100);
% y = y(1:100);
% 绘制这些点，'o'表示圆形，MarkerSize控制点的大小
% red_color = [1, 0, 0, 0.5];  % 红色半透明
scatter(tg_x, tg_y, diameter, 'filled',  'MarkerFaceAlpha', 0.5);
% 标记中心位置
rectangle('Position', [round(img_width/2)-5, round(img_height/2)-5, 10, 10], ...
          'EdgeColor', 'k', 'LineWidth', 0.7, 'Curvature', [1, 1]);  % 红色圆形
% 标记目标位置
% for ecc = [2,4,6]
% for ori = 0:45:315
%     tgWidth = ut.deg2pix(sess.expt.GaborWidth);
%     target_loc = ut.Pol2Rect([ecc,ori]);
%     target_loc = target_loc.*[1,-1]+[img_width, img_height]/2;
% 
%     rectangle('Position', [target_loc(1)-tgWidth/2, target_loc(2)-tgWidth/2, tgWidth, tgWidth], ...
%               'EdgeColor', 'k', 'LineWidth', 0.7, 'Curvature', [1, 1]);
% end
% end
% 标记背景范围
bgWidth = 634.65 ;
rectangle('Position', [img_width/2-bgWidth/2, img_height/2-bgWidth/2, bgWidth, bgWidth], ...
          'EdgeColor', 'k', 'LineWidth', 0.7, 'Curvature', [1, 1]);

% 设置坐标轴和图形属性
axis image;
axis([0 1920 0 1080]);
title('600个随机抽样的点');
xlabel('X轴');
ylabel('Y轴');
hold off;
%%
figure(2);
scatter(tg_x(1:100), tg_y(1:100), diameter, 'filled',  'MarkerFaceAlpha', 0.5);
% 标记中心位置
rectangle('Position', [round(img_width/2)-5, round(img_height/2)-5, 10, 10], ...
          'EdgeColor', 'k', 'LineWidth', 0.7, 'Curvature', [1, 1]);  % 红色圆形
% 标记目标位置
% for ecc = [2,4,6]
% for ori = 0:45:315
%     tgWidth = ut.deg2pix(sess.expt.GaborWidth);
%     target_loc = ut.Pol2Rect([ecc,ori]);
%     target_loc = target_loc.*[1,-1]+[img_width, img_height]/2;
% 
%     rectangle('Position', [target_loc(1)-tgWidth/2, target_loc(2)-tgWidth/2, tgWidth, tgWidth], ...
%               'EdgeColor', 'k', 'LineWidth', 0.7, 'Curvature', [1, 1]);
% end
% end
% 标记背景范围
bgWidth = ut.deg2pix(sess.expt.bgWidth);
rectangle('Position', [img_width/2-bgWidth/2, img_height/2-bgWidth/2, bgWidth, bgWidth], ...
          'EdgeColor', 'k', 'LineWidth', 0.7, 'Curvature', [1, 1]);

% 设置坐标轴和图形属性
axis image;
axis([0 1920 0 1080]);
title('前100个随机抽样的点');
xlabel('X轴');
ylabel('Y轴');
hold off;

figure(3);
scatter(tg_x(1:200), tg_y(1:200), diameter, 'filled',  'MarkerFaceAlpha', 0.5);
% 标记中心位置
rectangle('Position', [round(img_width/2)-5, round(img_height/2)-5, 10, 10], ...
          'EdgeColor', 'k', 'LineWidth', 0.7, 'Curvature', [1, 1]);  % 红色圆形
% 标记目标位置
% for ecc = [2,4,6]
% for ori = 0:45:315
%     tgWidth = ut.deg2pix(sess.expt.GaborWidth);
%     target_loc = ut.Pol2Rect([ecc,ori]);
%     target_loc = target_loc.*[1,-1]+[img_width, img_height]/2;
% 
%     rectangle('Position', [target_loc(1)-tgWidth/2, target_loc(2)-tgWidth/2, tgWidth, tgWidth], ...
%               'EdgeColor', 'k', 'LineWidth', 0.7, 'Curvature', [1, 1]);
% end
% end
% 标记背景范围
bgWidth = ut.deg2pix(sess.expt.bgWidth);
rectangle('Position', [img_width/2-bgWidth/2, img_height/2-bgWidth/2, bgWidth, bgWidth], ...
          'EdgeColor', 'k', 'LineWidth', 0.7, 'Curvature', [1, 1]);

% 设置坐标轴和图形属性
axis image;
axis([0 1920 0 1080]);
title('前200个随机抽样的点');
xlabel('X轴');
ylabel('Y轴');
hold off;

figure(6);
scatter(tg_x(1:30), tg_y(1:30), diameter, 'filled',  'MarkerFaceAlpha', 0.5);
% 标记中心位置
rectangle('Position', [round(img_width/2)-5, round(img_height/2)-5, 10, 10], ...
          'EdgeColor', 'k', 'LineWidth', 0.7, 'Curvature', [1, 1]);  % 红色圆形
% 标记目标位置
% for ecc = [2,4,6]
% for ori = 0:45:315
%     tgWidth = ut.deg2pix(sess.expt.GaborWidth);
%     target_loc = ut.Pol2Rect([ecc,ori]);
%     target_loc = target_loc.*[1,-1]+[img_width, img_height]/2;
% 
%     rectangle('Position', [target_loc(1)-tgWidth/2, target_loc(2)-tgWidth/2, tgWidth, tgWidth], ...
%               'EdgeColor', 'k', 'LineWidth', 0.7, 'Curvature', [1, 1]);
% end
% end
% 标记背景范围
bgWidth = ut.deg2pix(sess.expt.bgWidth);
rectangle('Position', [img_width/2-bgWidth/2, img_height/2-bgWidth/2, bgWidth, bgWidth], ...
          'EdgeColor', 'k', 'LineWidth', 0.7, 'Curvature', [1, 1]);

% 设置坐标轴和图形属性
axis image;
axis([0 1920 0 1080]);
title('前60个随机抽样的点');
xlabel('X轴');
ylabel('Y轴');
hold off;