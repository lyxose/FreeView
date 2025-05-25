function [oper, texture] = showInstruc(wptr, instName, instFolder, nextKey, backKey)
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
    Screen('Flip', wptr);
    
    % Prevent skipping the first page (wait the 'nextKey' released)
    [~, ~, keyCode] = KbCheck;
    while keyCode(KbName(nextKey))
        [~, ~, keyCode] = KbCheck;
    end
    % Wait for a key press
    oper = 0;
    while 1
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
                Screen('Flip', wptr);
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
