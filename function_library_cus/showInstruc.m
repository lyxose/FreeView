function [oper, texture] = showInstruc(wptr, instName, instFolder, nextKey, backKey, minSlideT)
    % Show instruction images in sequence based on filename order.
    % 
    % This function reads instruction image files with names containing 
    % instName from the specified folder, presents them in sequence, and 
    % listens for keyboard input to navigate through the images. The 
    % function supports navigation using specified next and back keys.
    % 
    % Parameters:
    % wptr        : Window pointer (returned by Screen('OpenWindow'))
    % instName    : Partial name to match instruction image files
    % instFolder  : Folder path containing the instruction image files
    % nextKey     : Key name for navigating to the next image (e.g., 'space')
    % backKey     : Key name for navigating to the previous image (e.g., 'backspace')
    % minSlideT   : Minimal reading duration (sec) of each slide, default 0
    %               (except the last one when there are more than 1)
    % 
    % Returns:
    % oper        : Operation result (1 for next, -1 for back, 0 for exit)
    % texture     : Texture object of the currently displayed image

    % Read the list of image file names from the specified folder
    afiles = dir(fullfile(instFolder, 'EmptyMatchedStruct'));
    for ftype={'*.png', '*.jpg', '*.jpeg'}
        files = dir(fullfile(instFolder, ftype{1})); % Assuming JPEG images
        afiles = [afiles; files];
    end
    fileNames = {afiles.name};
    
    % Match the file path based on the given name
    matchedFiles = fileNames(contains(fileNames, instName));
    
    % Check if there is at least one match
    if isempty(matchedFiles)
        error('No file matched found for the instruction image.');
    end
    
    % Get the first matched image file path and draw it
    imgPointer = 1;
    imgPath = fullfile(instFolder, matchedFiles{imgPointer});
    texture = drawCentImg(wptr, imgPath, 'fit');
    t0 = Screen('Flip', wptr);
    
    % Wait for a key press
    while 1
        oper = 0;
        while nargin>5 && GetSecs-t0<minSlideT
            WaitSecs(0.1);  % wait until reach the minimal time of each slide
        end
        [keyIsDown, ~, keyCode] = KbCheck;
        if keyIsDown
            if keyCode(KbName(nextKey))
                oper = 1; % Next key pressed
            elseif keyCode(KbName(backKey))
                oper = -1; % Back key pressed
            end
            imgPointer = imgPointer+oper;
            if imgPointer>0 && imgPointer<=length(matchedFiles)
                imgPath = fullfile(instFolder, matchedFiles{imgPointer});
                texture = drawCentImg(wptr, imgPath, 'fit');
                t0 = Screen('Flip', wptr);
            else
                break
            end
        end
        checkend;
        while keyIsDown % wait key release
            keyIsDown = KbCheck;
            WaitSecs(0.1);
        end
        WaitSecs(0.1);
    end
end
