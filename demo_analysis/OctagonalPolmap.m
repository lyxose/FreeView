function map = OctagonalPolmap(img_width, img_height, sigma_r, sigma_theta, u_r)
% 8向极坐标高斯分布 叠加 中央高斯分布   
% 参数设置
    center = [img_width/2, img_height/2];  % 视野中心坐标
%     max_radius = norm(center);             % 最大可能半径
    angles = deg2rad(0:45:315);             % 八方向角度(弧度)
    
    % 高斯分布参数 (基于人类眼动特性)
%     sigma_r = max_radius * 0.15;    % 径向标准差 (与中央凹衰减相关)
%     sigma_theta = deg2rad(12);      % 角度标准差 (与方位选择性相关)
if nargin < 5   
    u_r = 0; % 峰值距离中心的典型眼跳幅度
end

    % 生成极坐标系网格
    [X, Y] = meshgrid(1:img_width, 1:img_height);
    [theta, r] = cart2pol(X - center(1), Y - center(2));
    
    % 构建混合高斯模型
    mixmap = zeros(size(X));
    for k = 1:8
        % 角度差异计算 (处理周期性)
        angle_diff = mod(theta - angles(k) + pi, 2*pi) - pi;  
        
        % 径向衰减因子 (Gamma函数模拟中央凹效应)
%         radial_term = (r.^2) .* exp(-r/(peak_distance/3)) ./ gamma(3);
        
        % 高斯分量
        gaussian = exp(-(angle_diff.^2)/(2*sigma_theta^2)...
                     - (r - u_r).^2/(2*sigma_r^2));
        
        mixmap = mixmap + gaussian ;%.* radial_term;
    end
    
    % 正则化并生成热图
    map = mixmap / max(mixmap(:));
    
    % 可视化
    figure
    imagesc(map)
%     colormap(flipud(hot))  % 暖色系模拟典型眼动热图
    colorbar
    axis equal tight
    title('Octagonal Gaussian Mixture Heatmap')
end