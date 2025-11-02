function [texture, CentrRect] = drawCentImg(wptr, imgpath, zoom)
    % Draw the given image at the center of the screen by the given zoom way.
    % 
    % Parameters:
    % wptr:    window pointer (returned by Screen('OpenWindow'))
    % imgpath: path of image file (any format supported by imread)
    % zoom:    how to zoom the picture, can be numeric or a string
    %        - numeric:   scale factor (e.g., 0.5 to reduce size by half)
    %        - 'fill':    scale the image to fill the entire screen while maintaining aspect ratio
    %        - 'fit':     scale the image to fit within the screen while maintaining aspect ratio
    %        - 'stretch': scale the image to fill the entire screen, ignoring aspect ratio
    %        - 'tile':    display the image at its original size, tiling it if necessary

    % Load the image
    [img,cmap] = imread(imgpath); % Replace with your image file name
    if ~isempty(cmap)     % gray level img (single channel)
        img = ind2rgb(img,cmap);
        img = im2uint8(img);
    end

    % Check and handle EXIF metadata if necessary
    % This step is optional and depends on the specific image file
    % For example, you can use the 'exifread' function to read EXIF metadata
    % and adjust the image orientation accordingly
    % Read image metadata using imfinfo
    info = imfinfo(imgpath);
    try
        orientation = info.Orientation;
        
        % Adjust image orientation based on EXIF metadata
        switch orientation
            case 2 % Horizontal flip
                img = fliplr(img);
            case 3 % Rotate 180 degrees
                img = rot90(img, 2);
            case 4 % Vertical flip
                img = flipud(img);
            case 5 % Transpose and horizontal flip
                img = rot90(fliplr(img));
            case 6 % Rotate 90 degrees clockwise
                img = rot90(img, -1);
            case 7 % Transpose and vertical flip
                img = rot90(flipud(img));
            case 8 % Rotate 90 degrees counterclockwise
                img = rot90(img);
        end
    catch me
        disp(me)
        disp('There may not be any orientation information in this img file')
    end

    % Create a texture object from the image
    texture = Screen('MakeTexture', wptr, img);
    
    % Calculate the presentation size
    [texH, texW, ~] = size(img);
    [scrW, scrH] = Screen('WindowSize', wptr);
    
    if isnumeric(zoom)
        % Scale the image by the specified factor
        newW = texW * zoom;
        newH = texH * zoom;
    elseif strcmp(zoom, 'fill')
        % Fill mode: scale the image to fill the entire screen while maintaining aspect ratio
        scale = max(scrW / texW, scrH / texH);
        newW = texW * scale;
        newH = texH * scale;
    elseif strcmp(zoom, 'fit')
        % Fit mode: scale the image to fit within the screen while maintaining aspect ratio
        scale = min(scrW / texW, scrH / texH);
        newW = texW * scale;
        newH = texH * scale;
    elseif strcmp(zoom, 'stretch')
        % Stretch mode: scale the image to fill the entire screen, ignoring aspect ratio
        newW = scrW;
        newH = scrH;
    elseif strcmp(zoom, 'tile')
        % Tile mode: display the image at its original size, tiling it if necessary
        newW = texW;
        newH = texH;
    else
        error('Invalid zoom mode');
    end
    
    % Calculate the center position
    xc = (scrW - newW) / 2;
    yc = (scrH - newH) / 2;
    
    % Draw the texture at the calculated position and size
    CentrRect = [xc, yc, xc + newW, yc + newH];
    Screen('DrawTexture', wptr, texture, [], CentrRect);
    
    if strcmp(zoom, 'tile')
        % Calculate the number of tiles needed in each direction
        numTilesX = ceil(scrW / newW);
        numTilesY = ceil(scrH / newH);
        
        % Draw the texture in a tiled pattern
        for i = -ceil(numTilesX / 2):ceil(numTilesX / 2)
            for j = -ceil(numTilesY / 2):ceil(numTilesY / 2)
                x = xc + i * newW;
                y = yc + j * newH;
                Screen('DrawTexture', wptr, texture, [], [x, y, x + newW, y + newH]);
            end
        end
    end
end
    
   