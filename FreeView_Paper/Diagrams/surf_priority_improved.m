%% ============ 3D 地形图可视化系统 ============
% 功能: 在真实图片上叠加三维地形效果，支持视觉显著性计算、
%       手动峰指定、高度调节、防穿透等功能
% 
% 使用说明:
% 1. 修改"配置参数"段落中的参数
% 2. 运行脚本
% 3. 在命令行中修改参数后重新运行可快速预览效果

%========== 配置参数 ==========
%% 1. 图片与显著性地图
img_path = "E:\AllDownloads\DOWNLOAD\ChatGPT Image 2026年2月3日 14_53_47.png";  % 替换为你的图片路径
img = imread(img_path);
img = im2double(img);                % 转为 double，方便显示

%% 2. 显著性计算与地形参数
use_saliency = true;                 % true:根据图像视觉显著性计算; false:使用随机矩阵
saliency_method = 'contrast';        % 显著性方法: 
                                     % 'contrast'   - 基于梯度对比度（推荐）
                                     % 'edges'      - 基于边缘检测
                                     % 'intensity'  - 基于局部方差
                                     % 'frequency'  - 基于频域分析
gaussian_sigma = 60;                  % 高斯平滑标准差，越大越光滑

%% 3. 手动设置峰位置（可选）
% 格式：[pixel_x, pixel_y, height_scale, local_sigma]
% 不设置则置为空矩阵 []
% 示例（请根据实际图片坐标修改）:
manual_peaks = [];
% manual_peaks = [238,497,0.9,37;
%                 238,556,0.9,37;
%                 238,627,0.9,37;
%                 842,552,1.1,65;
%                 330,900,0.6,55;
%                 460,900,0.6,55;
%                 590,900,0.6,55;
%                 720,900,0.6,55;
%                 330,1000,0.6,55;
%                 460,1000,0.6,55;
%                 590,1000,0.6,55;
%                 720,1000,0.6,55];

% manual_peaks = [100, 150, 1.5, 10; 200, 250, 2.0, 15; 300, 180, 1.2, 8];

%% 4. 三维地形高度缩放与z偏移
z_scale = 1;                       % 高度缩放因子
                                     % < 1.0: 地形较平缓
                                     % = 1.0: 正常高度
                                     % > 1.0: 地形更陡峭 (推荐 0.5~2.0)

z_offset_method = 'custom';            % z偏移方法:
                                     % 'safe'   - 自动计算，确保不穿透底图
                                     % 'custom' - 使用自定义偏移值

z_offset_custom = 0.5;               % 当选择'custom'时使用的偏移值

% Z轴显示范围（重要！控制视觉缩放）
z_axis_max = 8;                    % 固定Z轴最大值
                                     % 改变此值查看不同的高度感：
                                     % 更小值(如0.2): 山峰显得更陡峭
                                     % 更大值(如0.6): 山峰显得更平缓

%% 5. 绘图模式
plot_mode = 'color';                 % 绘图模式:
                                     % 'color' - 表面渐变色填充（推荐）
                                     % 'grid'  - 网格线框架

%% 6. 表面属性
face_alpha = 0.4;                    % 表面透明度（0~1，0=完全透明，1=不透明）
color_map = 'parula';                % 颜色映射: 'parula', 'jet', 'cool', 'hot', 'viridis' 等
edge_alpha = 0.3;                    % 网格线透明度（仅grid模式，推荐 0.2~0.5）

%========== 处理流程 ==========

%% 步骤1: 获取图像尺寸
[H, W, channels] = size(img);
fprintf('\n========== 处理开始 ==========\n');
fprintf('图像尺寸: %d × %d 像素\n', H, W);

%% 步骤2: 计算显著性地图作为地形基础
if use_saliency
    fprintf('计算显著性地图 (%s方法)...', saliency_method);
    M = compute_saliency(img, saliency_method);
    fprintf(' 完成\n');
else
    fprintf('使用随机矩阵作为地形...\n');
    M = rand(H, W);
end

%% 步骤3: 二维高斯平滑
fprintf('应用高斯平滑 (Sigma=%.1f)...', gaussian_sigma);
M_smooth = imgaussfilt(M, gaussian_sigma);
M_smooth = M_smooth(1:H, 1:W);
fprintf(' 完成\n');

%% 步骤4: 归一化到 [0, 1]
M_smooth = (M_smooth - min(M_smooth(:))) / (max(M_smooth(:)) - min(M_smooth(:)) + eps);

%% 步骤5: 添加手动指定的峰
if ~isempty(manual_peaks)
    fprintf('添加手动峰点 (%d 个)...\n', size(manual_peaks, 1));
    M_smooth = add_manual_peaks(M_smooth, manual_peaks, H, W);
    % 再次归一化
    M_smooth = (M_smooth - min(M_smooth(:))) / (max(M_smooth(:)) - min(M_smooth(:)) + eps);
end

%% 步骤6: 应用高度缩放
fprintf('应用高度缩放 (z_scale=%.2f)...', z_scale);
M_scaled = M_smooth * z_scale;
fprintf(' 完成\n');

%% 步骤7: 计算安全的 z_offset（防穿透）
if strcmp(z_offset_method, 'safe')
    % 确保地形最小值不会穿透底图
    z_min = min(M_scaled(:));
    z_max = max(M_scaled(:));
    
    % 方法：将地形重新映射到 [z_base, z_base + z_range] 范围
    % 这样可以避免任何穿透
    z_base = 0.05;              % 底部离底图的距离
    z_range = 0.3;              % 地形在z方向的厚度
    
    if z_max > z_min
        M_scaled = (M_scaled - z_min) / (z_max - z_min) * z_range + z_base;
    else
        M_scaled = ones(size(M_scaled)) * z_base;
    end
    z_offset = 0;
    fprintf('使用安全z偏移（范围：%.3f~%.3f）\n', z_base, z_base + z_range);
else
    fprintf('使用自定义z偏移 (%.4f)\n', z_offset_custom);
    z_offset = z_offset_custom;
end

%% 步骤8: 构建坐标网格
[x, y] = meshgrid(1:W, 1:H);

%% 步骤9: 绘图
fprintf('绘制地形图...\n');
figure('Color','w', 'NumberTitle','off', 'Name','3D Terrain Overlay on Image');
axes('Position',[0 0 1 1]);
hold on;

% 背景：真实图片在 z=0 平面
surf(x, y, zeros(H, W), ...
     img, ...
     'EdgeColor', 'none', ...
     'FaceColor', 'texturemap');

% 前景：三维地形
z_terrain = M_scaled + z_offset;

if strcmp(plot_mode, 'color')
    % 表面渐变色模式
    fprintf('绘图模式: 表面渐变色\n');
    surf(x, y, z_terrain, ...
         'EdgeColor', 'none', ...
         'FaceAlpha', face_alpha, ...
         'FaceColor', 'interp', ...
         'CData', M_scaled);
        colormap(color_map);
        caxis([min(M_scaled(:)), max(M_scaled(:))]); % 自动调整颜色映射范围
        colorbar('Location', 'eastoutside');
    
elseif strcmp(plot_mode, 'grid')
    % 网格线模式
    fprintf('绘图模式: 网格线\n');
    surf(x, y, z_terrain, ...
         'EdgeColor', 'black', ...
         'EdgeAlpha', edge_alpha, ...
         'FaceColor', 'white', ...
         'FaceAlpha', 0.2);
        colormap(color_map);
        caxis([min(M_scaled(:)), max(M_scaled(:))]); % 自动调整颜色映射范围
else
    warning('未知的绘图模式，使用默认color模式\n');
    surf(x, y, z_terrain, ...
         'EdgeColor', 'none', ...
         'FaceAlpha', face_alpha, ...
         'FaceColor', 'interp', ...
         'CData', M_scaled);
    colormap(color_map);
end

shading interp;

%% 步骤10: 视角、透视、坐标轴设置
view(3);                             % 三维视角
camproj('perspective');              % 透视投影

axis tight;
axis off;                            % 取消坐标轴
set(gca, 'YDir', 'reverse');         % 让图像坐标与屏幕一致（左上为原点）

% 固定Z轴范围，确保z_scale调节能改变视觉效果
zlim([0, z_axis_max]);               % 固定Z轴范围，使地形在绝对坐标系中变化

% 光照效果（增强三维感）
camlight headlight;
lighting gouraud;
material dull;
%% 步骤11: 显示参数总结
fprintf('========== 参数配置总结 ==========\n');
fprintf('显著性方法......: %s\n', saliency_method);
fprintf('高斯Sigma......: %.1f\n', gaussian_sigma);
fprintf('高度缩放.......: %.2f\n', z_scale);
fprintf('Z轴偏移.......: %.4f\n', z_offset);
fprintf('Z轴显示范围...: [0, %.3f]  (固定)\n', z_axis_max);
fprintf('绘图模式.......: %s\n', plot_mode);
fprintf('表面透明度.....: %.2f\n', face_alpha);
fprintf('手动峰个数.....: %d\n', size(manual_peaks, 1));
fprintf('地形高度范围...: [%.4f, %.4f]\n', min(z_terrain(:)), max(z_terrain(:)));
fprintf('==================================\n');
fprintf('💡 改变z_scale会改变地形在[0,%.3f]范围内的高度\n', z_axis_max);
fprintf('💡 改变z_axis_max会改变整体高度感（更小→更陡峭）\n\n');

%========== 辅助函数 ==========

function saliency = compute_saliency(img, method)
    % 计算图像视觉显著性
    % 
    % 输入:
    %   img (H × W × C) - 输入图像（.png/.jpg等）
    %   method (str)    - 显著性计算方法
    %
    % 输出:
    %   saliency (H × W) - 显著性地图（值越大表示显著性越高）
    
    if size(img, 3) > 1
        img_gray = rgb2gray(img);
    else
        img_gray = img;
    end
    
    switch method
        case 'contrast'
            % 基于梯度对比度（最常用，效果好）
            [Gx, Gy] = imgradient(img_gray);
            saliency = sqrt(Gx.^2 + Gy.^2);
            
        case 'edges'
            % 基于边缘检测（Canny算子）
            saliency = edge(img_gray, 'Canny');
            saliency = imgaussfilt(double(saliency), 3);
            
        case 'intensity'
            % 基于局部强度方差（检测纹理区域）
            local_var = stdfilt(img_gray, ones(9));
            saliency = local_var;
            
        case 'frequency'
            % 基于频域显著性（复杂方法，更能捕捉周期性）
            F = fft2(img_gray);
            F_shifted = fftshift(F);
            A = abs(F_shifted);
            P = angle(F_shifted);
            R = A .^ 0.4;          % 幂律变换（AC分量突出）
            F_reconstructed = (R .* exp(1i * P));
            saliency = abs(ifft2(ifftshift(F_reconstructed)));
            
        otherwise
            % 默认使用对比度方法
            [Gx, Gy] = imgradient(img_gray);
            saliency = sqrt(Gx.^2 + Gy.^2);
    end
    
    % log 变换：避免极端值主导，增强中间值
    saliency = log(saliency + 1);
end

function M = add_manual_peaks(M, peaks_spec, H, W)
    % 在指定位置添加高斯锥形（峰）
    %
    % 输入:
    %   M (H × W)        - 原始地形矩阵
    %   peaks_spec (N × 4) - 峰的规格矩阵
    %                        第1列: pixel_x (1~W)
    %                        第2列: pixel_y (1~H)
    %                        第3列: height_scale (相对高度，推荐 0.5~3.0)
    %                        第4列: local_sigma (高斯宽度，越大越宽泛)
    %   H, W              - 图像高度与宽度
    %
    % 输出:
    %   M (H × W)        - 添加峰后的地形矩阵
    %
    % 示例:
    %   peaks_spec = [100, 150, 1.5, 10;    % 在(100,150)处添加高度1.5、宽度10的峰
    %                 200, 250, 2.0, 15];   % 在(200,250)处添加高度2.0、宽度15的峰
    
    [xx, yy] = meshgrid(1:W, 1:H);
    
    for i = 1:size(peaks_spec, 1)
        peak_x = peaks_spec(i, 1);
        peak_y = peaks_spec(i, 2);
        peak_height = peaks_spec(i, 3);
        peak_sigma = peaks_spec(i, 4);
        
        % 验证坐标范围
        if peak_x < 1 || peak_x > W || peak_y < 1 || peak_y > H
            warning('峰#%d 坐标超出范围 (%.0f, %.0f)，已跳过', i, peak_x, peak_y);
            continue;
        end
        
        % 构造二维高斯函数（钟形山峰）
        gaussian = peak_height * exp(-((xx - peak_x).^2 + (yy - peak_y).^2) / (2 * peak_sigma^2));
        
        % 添加到地形
        M = M + gaussian;
        
        fprintf('  已添加峰#%d: 位置(%.0f, %.0f), 高度%.2f, 宽度%.1f\n', ...
                i, peak_x, peak_y, peak_height, peak_sigma);
    end
end
