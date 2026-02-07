function img = winOverlap(background, source, widthPix, centerLoc, windowType, show)
    % 生成圆形窗口并将源图像与背景线性混合
    % Generate a circular window and linearly blend source with background
    %
    % 在中心处窗口值为1，在指定宽度处逐渐降至0。然后，将加窗的源矩阵
    % 与背景线性混合（即，窗口中心的源像素具有更高的窗口值，使其更可见）。
    %
    % The window has a value of 1 at the center, gradually decreasing to 0 
    % at the specified width. Then, linearly blend the windowed source matrix 
    % with the background (i.e., pixels at the window center of the source 
    % have a higher window value, making them more visible).
    %
    % -----Parameters-----
    % background: 背景图像矩阵（外围部分）
    %             Background image matrix (periphery part)
    % source: 源图像矩阵（中心部分）。源图像对齐到背景的左上角。
    %         如果源图像小于背景，则用中灰值(0.5)填充到相同大小。
    %         Source image matrix (central part). Aligned to upper-left corner.
    %         If smaller than background, padded with middle gray (0.5).
    % widthPix: 窗口宽度，由半高全宽(FWHM)确定
    %           Window width, determined by full width at half maximum (FWHM)
    % centerLoc: 窗口中心位置 [x, y] (像素)。[1,1]表示背景的左上角。
    %            默认为背景中心。
    %            Window center location [x, y] (pixels). [1,1] is upper-left.
    %            Defaults to background center.
    % windowType: 窗口类型: 'cos', 'Gaussian', 'linear', 'hard'。默认'cos'。
    %             Window type: 'cos', 'Gaussian', 'linear', 'hard'. Default 'cos'.
    % show: 是否显示生成的图像，默认false
    %       Whether to display generated image, default false
    %
    % -----Returns-----
    % img: 混合后的图像矩阵
    %      Blended image matrix
    %
    % Yuxin Lu, IPCAS
    % luyx@psych.ac.cn
    % 2025.1.5

    if nargin < 6
        show = false;
    end
    if nargin < 5
        windowType = 'cos';
    end

    bgSize = flip(size(background)); % [x, y]
    if nargin < 4
        centerLoc = round(bgSize/2);
    end
    
    % 裁剪并填充源图像以适应背景矩阵的大小
    % Cut and pad source to fit the size of background matrix
    ss = size(source);
    bs = size(background);
    ss = min(ss, bs);
    source = source(1:ss(1), 1:ss(2));
    if any(bs > ss)
        % 用0.5填充源图像以匹配背景大小
        % Pad source with 0.5 to match the size of background
        source = padarray(source, bs - ss, 0.5, 'post');
    end

    X = -centerLoc(1) + 1 : 1 : bgSize(1) - centerLoc(1); % 中心为0
    Y = -centerLoc(2) + 1 : 1 : bgSize(2) - centerLoc(2); % 中心为0
    
    [Xm, Ym] = meshgrid(X, Y);
    
    % 定义窗口
    % Define the window
    d_E = sqrt(Xm.^2 + Ym.^2); % 欧几里得空间中的距离 / Euclidean distance
    
    if strcmp(windowType, 'cos')
        % 余弦窗口 / Cosine window
        window = cos(pi/3 * 2/widthPix * d_E);
        window(d_E >= 1.5*widthPix/2) = 0; % 只取中心半周期
    elseif strcmp(windowType, 'Gaussian')
        % 高斯窗口 / Gaussian window
        sigma = widthPix / 2.355; % FWHM to sigma conversion
        window = exp(-(d_E.^2) / (2 * sigma^2));
        window = (window - min(window(:))) / (max(window(:)) - min(window(:)));
    elseif strcmp(windowType, 'linear')
        % 线性窗口 / Linear window
        window = -abs(d_E);
        window(d_E >= widthPix) = -widthPix;
        window = (window - min(window(:))) / (max(window(:)) - min(window(:)));
    elseif strcmp(windowType, 'hard')
        % 硬边窗口 / Hard-edge window
        window = d_E <= widthPix/2;
    else
        error('未知的窗口类型 / Unknown window type specified.');
    end

    % 混合背景与源图像
    % Mix background with source
    img = (1 - window) .* background + window .* source;

    % 可选显示结果
    % Optionally show the result
    if show
        figure;
        imshow(img, []);
        title('Overlay (叠加)');
        colorbar;
    end

end
