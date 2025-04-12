% 步骤1：创建网格坐标
function map = gaumap(img_width, img_height, sigma)
[X, Y] = meshgrid(1:img_width, 1:img_height);
X = X - img_width/2; % 将原点移至图像中心
Y = Y - img_height/2;

% 步骤2：定义一维正态分布函数
% sigma = 60; % 标准差
norm1d = @(x, sigma) exp(-(x - 0).^2 / (2 * sigma^2));

% 步骤3：生成四个一维正态分布
sigma_x =  sigma; % x方向均值
sigma_y =  sigma; % y方向均值
sigma_xy = sigma;
sigma_yx = sigma;               

% 计算x和y方向的正态分布
Z_x = norm1d(X, sigma_x);
Z_y = norm1d(Y, sigma_y);
Z_xy= norm1d((X*1+Y*1)/sqrt(2), sigma_xy);
Z_yx= norm1d((X*1-Y*1)/sqrt(2), sigma_yx);
% 
% % 步骤4：旋转并叠加正态分布
% theta = pi / 4; % 旋转角度45度
% rotation_matrix = [cos(theta), -sin(theta); sin(theta), cos(theta)];
% 
% % 旋转后的坐标
% rotated_X = rotation_matrix(1,1) * X + rotation_matrix(1,2) * Y;
% rotated_Y = rotation_matrix(2,1) * X + rotation_matrix(2,2) * Y;
% 
% % 计算旋转后的正态分布
% Z_rotated = norm1d(rotated_X, mu_x) .* norm1d(rotated_Y, mu_y);
% 
% % 步骤5：归一化处理
% Z_rotated = Z_rotated / sum(Z_rotated(:));

% 显示结果
map = Z_x+Z_y+Z_yx+Z_xy;
return
% Z(Z>1)=1;
% Z = Z-norm1d(sqrt(X.^2+Y.^2),80).*4;
% map((X.^2+Y.^2)<5000)=0;
% imagesc(map);
% colormap('cool');
% colorbar;
% axis equal;
% title('旋转叠加的正态分布图');
