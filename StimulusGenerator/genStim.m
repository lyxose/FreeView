
function stimulus = genStim(winRect, ut, bgContrast, tgContrast, tgCenter, GaborSF, GaborWidth, GaborOrient, bgWidth, seed)
    % 生成刺激图像（在圆形1/f粉噪声背景上叠加Gabor目标）
    % Generate stimulus image (Gabor target on circular 1/f pink noise background)
    %
    % -----Parameters-----
    % winRect: 窗口尺寸 [x0, y0, width, height]，或 [width, height]
    %          Window size [x0, y0, width, height], or [width, height]
    % ut: UT类实例，用于单位转换 / UT class instance for unit conversion
    % bgContrast: 背景对比度 (0-1) / Background contrast (0-1)
    % tgContrast: 目标对比度 (0-1) / Target contrast (0-1)
    % tgCenter: 目标中心位置 [x, y]，像素单位（图像坐标系）
    %           Target center [x, y] in pixels (image coordinate system)
    % GaborSF: Gabor空间频率 (cycles per degree)
    % GaborWidth: Gabor宽度 (degree, FWHM)
    % GaborOrient: Gabor方向 (degree)
    % bgWidth: 背景圆形区域宽度 (degree)
    % seed: 随机种子 / Random seed
    %
    % -----Returns-----
    % stimulus: 刺激图像矩阵，范围 [0, 1]
    %           Stimulus image matrix, range [0, 1]
    %
    % Yuxin Lu, IPCAS
    % luyx@psych.ac.cn
    % 2025.1.7

    if length(winRect) == 2
        scWidth = winRect(1);
        scHeight = winRect(2);
    else
        scWidth  = winRect(3);  
        scHeight = winRect(4);
    end
    bgCenter = [scWidth/2, scHeight/2];
    
    % 生成全屏1/f粉噪声背景
    % Generate full-screen 1/f pink noise background
    background = tPinkNoise(scWidth, seed, bgContrast);
    background = background(1:scHeight, 1:scWidth);
    
    % 计算Gabor参数并生成光栅
    % Calculate Gabor parameters and generate grating
    lambda = ut.deg2pix(1/GaborSF);
    if lambda ~= 0
        Texture = grating(size(background), tgCenter, ...
                            1/lambda, GaborOrient, tgContrast);
        Texture = winOverlap(background, Texture, ut.deg2pix(GaborWidth), ...
                              tgCenter, 'cos'); 
    else
        fprintf('空间频率过大 (%.4f)，无法在当前屏幕分辨率下生成光栅\n', GaborSF);
        fprintf('Spatial frequency too large (%.4f), cannot generate grating at current screen resolution\n', GaborSF);
        Texture = background;
    end
    
    % 裁剪为圆形背景区域
    % Crop to circular background area
    stimulus = winOverlap(zeros([scHeight, scWidth]) + 0.5, Texture, ...
                          ut.deg2pix(bgWidth), bgCenter, 'hard');
end
