function img = winOverlap(background, source, widthPix, centerLoc, windowType, show)
    % Generate a circular window with a value of 1 at the center, gradually
    % decreasing to 0.5 at the specified width. Then, linearly blend the 
    % windowed source matrix with the background (i.e., pixels at the 
    % window center of the source have a higher window value, making them
    % more visible).
    %
    % -----Parameters-----
    % background: The gray-level background image matrix (periphery part).
    % source: The gray-level source image matrix (central part). The source
    %   is aligned to the upper-left corner of the background. If the 
    %   source is smaller than the background, it is padded to the same 
    %   size using a middle gray value (0.5).
    % widthPix: The width of the window, determined by the full width at 
    %   half maximum (FWHM).
    % centerLoc: The center location of the window, array of 2 pixels. 
    %   [1,1] denotes the upper-left corner of the background. Defaults to 
    %   the center of the background.
    % windowType: Specifies the type of window to use: 'cos', 'Gaussian', 
    %   'linear' or 'hard'. Defaults to 'cos'.
    % show: Displays a figure of the generated image if set to true.
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

    bgSize = flip(size(background)); % [x,y]
    if nargin < 4
        centerLoc = round(bgSize/2);
    end
    % Cut and pad source to fit the size of background matrix
    ss = size(source);
    bs = size(background);
    ss = min(ss,bs);
    source = source(1:ss(1),1:ss(2));
    if any(bs > ss)
        % Pad source with 0.5 to match the size of background
        source = padarray(source, bs - ss, 0.5, 'post');
    end

    X = -centerLoc(1)+1 : 1 : bgSize(1)-centerLoc(1); % 0 at center
    Y = -centerLoc(2)+1 : 1 : bgSize(2)-centerLoc(2); % 0 at center
    
    [Xm, Ym] = meshgrid(X, Y);
    
    % Define the window
    d_E = sqrt(Xm.^2+Ym.^2); % distance in Euclidean space
    if strcmp(windowType, 'cos')
        window = cos(pi/3 * 2/widthPix * d_E);
        window(d_E >= 1.5*widthPix/2) = 0; % only the central half-cycle should be taken
    elseif strcmp(windowType, 'Gaussian')
        sigma = widthPix / 2.355; % FWHM to sigma conversion for Gaussian window
        window = exp(-(d_E.^2) / (2 * sigma^2));
        window = (window-min(window(:))) /(max(window(:))-min(window(:))); % normalize to 0~1
    elseif strcmp(windowType, 'linear')
        window = -abs(d_E);
        window(d_E >= widthPix) = -widthPix; % only the area within 2*width should be taken
        window = (window-min(window(:))) /(max(window(:))-min(window(:))); % normalize to 0~1
    elseif strcmp(windowType, 'hard')
        window = d_E <= widthPix/2;
    else
        error('Unknown window type specified.');
    end

    % Mix background with Gabor
    img = (1-window).*background + window.*source;

    % Optionally show the result
    if show
        figure;
        imshow(img, []);
        title('Overlay');
    end

end