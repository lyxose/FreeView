function quit = disp_rest(w,xc,yc,rest_time)

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
    
    WaitSecs(0.2);
end
