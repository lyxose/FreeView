function grating = grating(sizePix, centerLoc, freqPix, orientDeg, contrast, phaseRad, show)
    % 生成正弦光栅图像
    % Generate a sinusoidal grating image
    %
    % -----Parameters-----
    % sizePix: 图像尺寸 [height, width] (像素)
    %          Image size [height, width] (pixels)
    % centerLoc: Gabor中心位置 [x, y] (像素), [0,0]为左上角
    %            Gabor center location [x, y] (pixels), [0,0] is upper-left
    % freqPix: 空间频率 (cycles per pixel)
    %          Spatial frequency (cycles per pixel)
    % orientDeg: 方向 (度数)
    %            Orientation (degrees)
    % contrast: 对比度，灰度振幅。0.5限制值在[0.25, 0.75]之间
    %           Contrast, gray value amplitude. 0.5 restricts values to [0.25, 0.75]
    % phaseRad: 正弦波相位 (弧度), 默认0
    %           Sine wave phase (radians), default 0
    % show: 是否显示图像, 默认false
    %       Whether to display image, default false
    %
    % -----Returns-----
    % grating: 光栅图像矩阵，范围[0, 1]
    %          Grating image matrix, range [0, 1]
    %
    % Yuxin Lu, IPCAS
    % luyx@psych.ac.cn
    % 2025.1.7

    if nargin < 7
        show = false;
    end

    if nargin < 6
        phaseRad = 0;
    end

    % 转换方向为弧度
    % Convert orientation to radians
    gaborOrientation = deg2rad(orientDeg);

    % 生成光栅
    % Generate grating
    bgSize = flip(sizePix); % [x, y]
    
    X = -centerLoc(1) + 1 : 1 : bgSize(1) - centerLoc(1); % 中心为0
    Y = -centerLoc(2) + 1 : 1 : bgSize(2) - centerLoc(2); % 中心为0
    
    [Xm, Ym] = meshgrid(X, Y);
    
    % 计算正交方向的距离
    % Calculate distance in orthogonal direction
    d_O = Xm * cos(gaborOrientation) + Ym * sin(gaborOrientation);
    
    % 生成全屏光栅
    % Generate full-screen grating
    grating = (1 + contrast * sin(2*pi*freqPix*d_O + phaseRad)) / 2;

    % 可选显示结果
    % Optionally show the result
    if show
        figure;
        imshow(grating, []);
        title('Grating (光栅)');
        colorbar;
    end

end
