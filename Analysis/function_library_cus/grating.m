function grating = grating(sizePix, centerLoc, freqPix, orientDeg, contrast, phaseRad, show)
    % Draw a grating in a given sizePix
    %
    % -----Parameter-----
    % size: [x,y] for horizontal and vertical pixel number, respectively.
    % centerLoc: center location of gabor [x, y], by pixel, [0,0] for
    %   upleft corner of background
    % freqPix: spatial frequency of Gabor (cycles per unit length)
    % orientDeg: orientation of Gabor in degrees
    % contrast: gray value amplitude. 0.5 restricts the value between [0.25,0.75]
    % phaseRad: sine wave phase at the centerLoc
    % show: bool, to show a figure of the generated image if setted as true. 
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


    debug = false;
    if debug
        sizePix = [255,180];
        centerLoc = [50,90];
        widthPix = 15;
        freqPix = 1/8;
        orientDeg = 30;
        phaseRad = 0;
        windowType = 'Gaussian';
        show = true;
        contrast = 0.5;
    end
    % Gabor parameters
    gaborOrientation = deg2rad(orientDeg); % Convert orientation to radians

    % Generate grating
    bgSize = flip(sizePix); % [x,y]
    
    X = -centerLoc(1)+1 : 1 : bgSize(1)-centerLoc(1); % 0 at center
    Y = -centerLoc(2)+1 : 1 : bgSize(2)-centerLoc(2); % 0 at center
    
    [Xm, Ym] = meshgrid(X, Y);
    
    d_O = Xm*cos(gaborOrientation) + Ym*sin(gaborOrientation); % distance in orthogonal direction
    grating = (1 + contrast * sin(2*pi*freqPix*d_O + phaseRad))/2; % a full screen grating

    % Optionally show the result
    if show
        figure;
        imshow(grating, []);
        title('Grating');
    end

end