function quit = disp_rest(w,xc,yc,rest_time)
% DISPLAY REST INSTRUCTION SCREEN
% Usage: quit = disp_rest(w,xc,yc,rest_time)
% Inputs:   
%   w         - window pointer
%   xc, yc    - center coordinates
%   rest_time - duration of the rest period in seconds
% Outputs:
%   quit      - flag indicating if the user chose to quit (1) or continue (0)
% Example:
%   quit = disp_rest(w, xc, yc, 10); % Display rest screen for 10 seconds
% Note: == Will quit immediately when reach rest_time, with a time resolution of 0.05 seconds ==
start=GetSecs;
instFolder = './instructions/';
afiles = dir(fullfile(instFolder, 'EmptyMatchedStruct'));
for ftype={'*.png', '*.jpg', '*.jpeg'}
    files = dir(fullfile(instFolder, ftype{1})); % Assuming JPEG images
    afiles = [afiles; files];
end
fileNames = {afiles.name};

% Match the file path based on the given name
matchedFiles = fileNames(contains(fileNames, 'Rest'));
imgPath = fullfile(instFolder, matchedFiles{1});

[texture, rect] = drawCentImg(w, imgPath, 'fit');

while 1
    timeinfo=round(rest_time-(GetSecs-start));
    % drawcenteredtext_dot(w,'Take a Rest, press "c" to jump or "Esc" to quit',xc,yc,0,30);
    Screen('DrawTexture', w, texture, [], rect);
    drawcenteredtext_dot(w,num2str(timeinfo),xc,yc+60,0,30);
    Screen('Flip',w);
    [~,~,kc]=KbCheck;
    if kc(KbName('c')) || timeinfo<=0
        quit = 0;
        break;
    elseif kc(KbName('escape'))
        quit = 1;
        break
    end
    
    WaitSecs(0.05);
end
