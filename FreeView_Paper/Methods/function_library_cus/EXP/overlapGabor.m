function img = overlapGabor(background, centerLoc, widthPix, freqPix, orientDeg, contrast, phaseRad, windowType, show)
    % Draw a Gabor in a given gray-level image texture, with transparency of 
    % grating as 1 on center and gradually decreasing to 0.5 on width.
    %
    % -----Parameter-----
    % background: original texture or img matrix in gray level
    % centerLoc: center location of gabor [x, y], by pixel, [0,0] for
    %   upleft corner of background
    % width: width of Gabor, determined by full width at half maxima (FWHM) 
    % freqPix: spatial frequency of Gabor (cycles per unit length)
    % orientDeg: orientation of Gabor in degrees
    % contrast: gray value amplitude, defaulted to be 1.
    % phaseRad: sine wave phase at the centerLoc
    % windowType: specify the window to be 'cos' or 'Gaussian' or 'linear',
    %   defaulted to be 'cos'
    % show: bool, to show a figure of the generated image if setted as true. 
    % 
    % Yuxin Lu, IPCAS
    % luyx@psych.ac.cn
    % 2025.1.5

    if nargin < 9
        show = false;
    end
    if nargin < 8
        windowType = 'cos';
    end
    if nargin < 7
        phaseRad = 0;
    end

    debug = false;
    if debug
        background= zeros(255,255)+0.5;
        centerLoc = [50,90];
        widthPix = 15;
        freqPix = 1/8;
        orientDeg = 30;
        phaseRad = 0;
        windowType = 'Gaussian';
        show = true;
        contrast = 1;
    end
    % Gabor parameters
    gaborOrientation = deg2rad(orientDeg); % Convert orientation to radians

    % Generate grating
    bgSize = flip(size(background)); % [x,y]
    
    X = -centerLoc(1)+1 : 1 : bgSize(1)-centerLoc(1); % 0 at center
    Y = -centerLoc(2)+1 : 1 : bgSize(2)-centerLoc(2); % 0 at center
    
    [Xm, Ym] = meshgrid(X, Y);
    
    d_O = Xm*cos(gaborOrientation) + Ym*sin(gaborOrientation); % distance in orthogonal direction
    grating = (1 + sin(2*pi*freqPix*d_O + phaseRad))/2; % a full screen grating
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
        window(d_E >= widthPix/2) = -widthPix/2; % only the area within 2*width should be taken
        window = (window-min(window(:))) /(max(window(:))-min(window(:))); % normalize to 0~1
    else
        error('Unknown window type specified.');
    end

    % Mix background with Gabor
    img = (1-window).*background + window.*((grating-0.5)*contrast+0.5);

    % Optionally show the result
    if show
        figure;
        imshow(img, []);
        title('Gabor Overlay on Background');
    end

end