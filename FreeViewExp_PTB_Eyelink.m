%% FreeViewExp_PTB_Eyelink - EyeLink eye tracker experiment script
% FIXED VERSION - Resolves calibration blocking issues
% This script contains all EyeLink specific eye tracking code
% Called from FreeView_main.m
% Based on EyeLink_SimplePicture.m example

%% EyeLink-specific parameters
eyelinkFreq = 1000; % EyeLink sampling rate (Hz), can be 250, 500, 1000, or 2000 depending on model
% Initialize variables
monWidth = 51.1;  % Monitor width in cm 
monHeight = 28.7; % Monitor height in cm 
headDist = default_distance;
% Toggle optional pre-practice free-view block
ENABLE_FREEVIEW = true;      % set false to disable
FV_TRIALS = 20;              % number of free-view trials
FV_DURATION_SEC = 10;        % background duration per trial (s)
DoDriftCorrect = true;        % enable drift correction instead of show_fix

%% STEP 1: Initialize EyeLink connection and open EDF file
dummymode = 0; % 0 for real connection, 1 for dummy mode (testing without eye tracker)

% Initialize EyeLink connection
EyelinkInit(dummymode);
status = Eyelink('IsConnected');
if status < 1 % If EyeLink not connected
    dummymode = 1;
end

% Create EDF file name (max 8 characters)
edfFile = sprintf('S%dS%d', subjID, session);
if length(edfFile) > 8
    edfFile = edfFile(1:8); % Truncate if too long
end

% Open EDF file on Host PC
failOpen = Eyelink('OpenFile', edfFile);
if failOpen ~= 0
    error('Cannot create EDF file %s', edfFile);
end

% Get EyeLink tracker version
ELsoftwareVersion = 0;
[ver, versionstring] = Eyelink('GetTrackerVersion');
if dummymode == 0
    [~, vnumcell] = regexp(versionstring,'.*?(\d)\.\d*?','Match','Tokens');
    ELsoftwareVersion = str2double(vnumcell{1}{1});
    fprintf('Running experiment on %s version %d\n', versionstring, ver);
end

% Add preamble text to EDF file
preambleText = sprintf('RECORDED BY FreeView Experiment Subject %d Session %d', subjID, session);
Eyelink('Command', 'add_file_preamble_text "%s"', preambleText);

%% STEP 2: Configure sample/event data recording
Eyelink('Command', 'file_event_filter = LEFT,RIGHT,FIXATION,SACCADE,BLINK,MESSAGE,BUTTON,INPUT');
Eyelink('Command', 'link_event_filter = LEFT,RIGHT,FIXATION,SACCADE,BLINK,BUTTON,FIXUPDATE,INPUT');

if ELsoftwareVersion > 3
    Eyelink('Command', 'file_sample_data  = LEFT,RIGHT,GAZE,HREF,RAW,AREA,HTARGET,GAZERES,BUTTON,STATUS,INPUT');
    Eyelink('Command', 'link_sample_data  = LEFT,RIGHT,GAZE,GAZERES,AREA,HTARGET,STATUS,INPUT');
else
    Eyelink('Command', 'file_sample_data  = LEFT,RIGHT,GAZE,HREF,RAW,AREA,GAZERES,BUTTON,STATUS,INPUT');
    Eyelink('Command', 'link_sample_data  = LEFT,RIGHT,GAZE,GAZERES,AREA,STATUS,INPUT');
end

if DEBUGlevel>1
    % make screen partially transparent on OSX and windows vista or
    % higher, so we can debug.
    PsychDebugWindowConfiguration;
end
if DEBUGlevel
    Screen('Preference', 'SkipSyncTests', 1);
    Screen('Preference', 'Verbosity', 4);
else
    Screen('Preference', 'Verbosity', 2);
end

%% STEP 3: Open PTB window
if scr>1
    [wpnt,winRect] = PsychImaging('OpenWindow', 1, bgClr, [], [], [], [], 4);
    win_main = [0,0,768,432];
    win_edge = 30;
    [wpnt_main,winRect_main] = Screen('OpenWindow',0,[240,240,240],win_main+win_edge);
    Screen('BlendFunction', wpnt_main, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    alphaChannel = 180 .* ones(flip(win_main(3:4)), 'uint8');
else
    [wpnt,winRect] = PsychImaging('OpenWindow', 0, bgClr, [], [], [], [], 4);
end
hz=Screen('NominalFrameRate', wpnt);
Priority(1);
Screen('BlendFunction', wpnt, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
Screen('Preference', 'TextAlphaBlending', 1);
Screen('Preference', 'TextAntiAliasing', 2);
Screen('Preference', 'TextRenderer', 1);
KbName('UnifyKeyNames');

% Get max color value for EyeLink Host PC integration
colorMaxVal = Screen('ColorRange', wpnt);
[width, height] = Screen('WindowSize', wpnt);

% Ensure a clean first frame and clear any stale key presses before instructions
Screen('FillRect', wpnt, bgClr);
Screen('Flip', wpnt);
KbReleaseWait;

showInstruc(wpnt, 'welcome', instFolder, 'space', 'BackSpace');
showInstruc(wpnt, 'Calib', instFolder, 'space', 'BackSpace');

%% STEP 4: Setup EyeLink calibration
try
    ListenChar(-1);
catch ME
    ListenChar(2);
end

% Configure EyeLink display coordinates
Eyelink('Command','screen_pixel_coords = %ld %ld %ld %ld', 0, 0, width-1, height-1);
Eyelink('Message', 'DISPLAY_COORDS %ld %ld %ld %ld', 0, 0, width-1, height-1);

% EyeLink defaults
el = EyelinkInitDefaults(wpnt);
el.calibrationtargetsize = 3;
el.calibrationtargetwidth = 0.7;
% Use PTB-indexed colors to ensure visibility and correct color range
el.backgroundcolour = repmat(GrayIndex(wpnt),1,3);
el.calibrationtargetcolour = repmat(WhiteIndex(wpnt),1,3);
el.msgfontcolour = repmat(WhiteIndex(wpnt),1,3);

% Use custom calibration target if available
% Use default bull's eye calibration target to avoid missing image issues
% If you want to use a custom image, uncomment the lines below and ensure
% the image path is valid and accessible from the working directory.
% if exist('fixTarget.jpg', 'file')
%     el.calTargetType = 'image';
%     el.calImageTargetFilename = [pwd filesep 'fixTarget.jpg'];
% end

el.targetbeep = 1;
el.feedbackbeep = 1;

% Initialize PsychSound for calibration feedback
InitializePsychSound();
pamaster = PsychPortAudio('Open', [], 8+1);
PsychPortAudio('Start', pamaster);
pahandle = PsychPortAudio('OpenSlave', pamaster, 1);
el.ppa_pahandle = pahandle;

EyelinkUpdateDefaults(el);

% Set calibration type
Eyelink('Command', 'calibration_type = HV9');
Eyelink('Command', 'button_function 5 "accept_target_fixation"');

HideCursor(wpnt);
ListenChar(-1);
Eyelink('Command', 'clear_screen 0');

% 【FIX #1】Perform calibration - properly structured
if DEBUGlevel == 0
    % Ensure our window is frontmost and clean before entering setup
    Screen('Flip', wpnt);
    % Perform calibration
    EyelinkDoTrackerSetup(el);
end

% 【FIX #2】Restore keyboard listening to experiment window
ListenChar(0);

% 【IMPORTANT】Clear the offline mode and prepare for drift correction
Eyelink('SetOfflineMode');
Eyelink('Command', 'clear_screen 0');

% Optional free-view block before practice
if ENABLE_FREEVIEW
    run_freeview_eyelink(wpnt, winRect, el, edfFile, subjID, session, bgCenter, monWidth, monHeight, bgWidth, bgContrast, FV_TRIALS, FV_DURATION_SEC);
end


if length(bgContrast)~=1
    pre_bgContrast = mean(bgContrast);
else
    pre_bgContrast = bgContrast;
end

if threshold==0
    showInstruc(wpnt, 'Practice', instFolder, 'space', 'BackSpace');
    CtrGrad = 10.^ linspace(log10(0.6), log10(0.2), 10);
    passed = 0;
else
    passed = 1;
    ut = UT(monWidth, scWidth, headDist);
end

reCali = false;
while ~passed
    preresults = results(1:10,:);
    preresults.ECC = sqrt((bgWidth/2)^2 .* rand(10,1));
    preresults.Orient = rand(10,1).*360;
    preresults.Xtarg = preresults.ECC .* cosd(preresults.Orient);
    preresults.Ytarg = preresults.ECC .* sind(preresults.Orient);
    preresults.oriF = rand(10,1)*360;
    preresults.seed = randi(1000,10,1);
    
    for pretrial = 1:10
        ut = UT(monWidth, scWidth, headDist);
        fixCenter = ut.Pol2Rect([rFix,preresults.oriF(pretrial)]).*[1,-1]+bgCenter;

        % Start recording BEFORE fixation (show_fix needs live gaze)
        Eyelink('SetOfflineMode');
        Eyelink('Command', 'clear_screen 0');
        Eyelink('StartRecording');
        WaitSecs(0.1);
        Eyelink('Message', 'FIX ON Pre');
        if DoDriftCorrect
            EyelinkDoDriftCorrect(el, fixCenter(1), fixCenter(2), 1, 1); % last parameter = 1: allow recalibration
        else
            [startT, headDist] = show_fix(wpnt, fixCenter(1), fixCenter(2), fixTime, 0, winRect, MaxErr, ...
                                        'eyeTrackerType', 'EyeLink', 'el', el, ...
                                        'monWidth', monWidth, 'monHeight', monHeight);
        end

        % Recompute geometry with updated head distance
        ut = UT(monWidth, scWidth, headDist);
        tgCenter = ut.deg2pix([preresults.Xtarg(pretrial), preresults.Ytarg(pretrial)]);
        stimulus = genStim(winRect, ut, pre_bgContrast, CtrGrad(pretrial), ...
            tgCenter.*[1,-1]+bgCenter, GaborSF, GaborWidth, GaborOrient, bgWidth, preresults.seed(pretrial));
        stiTex = Screen('MakeTexture', wpnt, cat(3,stimulus,stimulus,stimulus)*255);
        
        Screen('Drawtexture',wpnt,stiTex);
        imgT = Screen('Flip',wpnt);
        Eyelink('Message', sprintf('STIM ON: P. trial-%.0f  ECC-%.2f  Ori-%.0f  tgContrast-%.3f  bgContrast-%.2f', ...
            pretrial, preresults.ECC(pretrial), preresults.Orient(pretrial), CtrGrad(pretrial), pre_bgContrast));
        
        timeCost = 0;
        
        % Check key response
        while 1
            checkend;
            [keyIsDown, secs, keyCode] = KbCheck;
            if keyCode(KbName(key1)) || timeCost>=maxTrialDur
                while 1
                    checkend;
                    [keyIsDown2, secs2, keyCode2] = KbCheck;
                    if keyCode2(KbName(key2)) || timeCost>=maxTrialDur
                        % 【FIX #5】Properly extract gaze data while recording
                        lastFixPix_ = bgCenter;
                        if Eyelink('NewFloatSampleAvailable') > 0
                            evt = Eyelink('NewestFloatSample');
                            % Prefer the eye selected by tracker, otherwise use any valid eye
                            eyeUsed = Eyelink('EyeAvailable');
                            if eyeUsed == el.BINOCULAR
                                eyeUsed = el.LEFT_EYE;
                            end
                            % gx/gy are 1-based indexed for left/right eyes
                            if eyeUsed == el.LEFT_EYE && evt.gx(1) ~= -32768 && evt.gy(1) ~= -32768
                                lastFixPix_ = [evt.gx(1), evt.gy(1)];
                            elseif eyeUsed == el.RIGHT_EYE && evt.gx(2) ~= -32768 && evt.gy(2) ~= -32768
                                lastFixPix_ = [evt.gx(2), evt.gy(2)];
                            else
                                % fallback: any valid eye
                                if evt.gx(1) ~= -32768 && evt.gy(1) ~= -32768
                                    lastFixPix_ = [evt.gx(1), evt.gy(1)];
                                elseif evt.gx(2) ~= -32768 && evt.gy(2) ~= -32768
                                    lastFixPix_ = [evt.gx(2), evt.gy(2)];
                                end
                            end
                        end
                        FixNum = 1; % Simplified for practice
                        
                        lastFixPix = lastFixPix_ - bgCenter;
                        lastFixPix(2) = -lastFixPix(2);
                        lastFixDeg = ut.pix2deg(lastFixPix);
                        break
                    end
                    timeCost = WaitSecs(0.01)-imgT;
                end
                if keyIsDown2 || timeCost>=maxTrialDur
                    break
                end
            end
            timeCost = WaitSecs(0.01)-imgT;
        end
        
        % Feedback
        err = norm(lastFixDeg - ut.pix2deg(tgCenter));
        if err<MaxErr
            pointColor = [0, 0, 255];
        else
            pointColor = [255, 255, 0];
        end
        
        Screen('Drawtexture',wpnt, stiTex);
        Screen('DrawDots', wpnt, lastFixPix_, 14, pointColor, [], 3);
        fbT = Screen('Flip',wpnt);
        WaitSecs(0.3);
        
        endT = Screen('Flip',wpnt);
        Screen('Close',stiTex);
        Eyelink('Message', sprintf('STIM OFF: FixNum-%.0f  Err-%.3f  judge-%.0f RT-%.3f', ...
            FixNum, err, err<MaxErr, secs-imgT));
        
        Eyelink('StopRecording');
        WaitSecs(rand(1)*0.4);
        
        % If error, stop early - start the next 10 trials until all correct
        if err>=MaxErr
            break
        end
        if pretrial == 10
            passed = 1;
        end
    end
    
    oper = showInstruc(wpnt, 'Check', instFolder, 'space', 'BackSpace');
    if oper==-1
        % redo calibration then redo practice
        Eyelink('SetOfflineMode');
        Eyelink('Command', 'clear_screen 0');
        EyelinkDoTrackerSetup(el);
        passed = 0;
        continue;
    else
        reCali = false; % passed
    end
end

WaitSecs(0.5);
oper = showInstruc(wpnt, 'T2F', instFolder, 'space', 'BackSpace');

%% Main experiment loop
for trial = 1:trialNum
    if trial ==1
        pThreshold = learnP;
        tGuessSd=1;
        tGuess = log10(tgContrast);
        grain = 0.01;
        beta=3;
        delta=0.01;
        gamma=2^2/bgWidth^2;
        range = 3;
        for i = 1:4
            q{i} = QuestCreate(tGuess,tGuessSd,pThreshold,beta,delta,gamma,grain,range);
            q{i}.normalizePdf = 1;
            q{i} = QuestUpdate(q{i}, tGuess,1);
            q{i} = QuestUpdate(q{i}, log10(0.005),0);
        end
    end
    if trial == learnTNum
        tGuess = tTest;
        for i = 1:4
            range = 2;
            q{i} = QuestCreate(tGuess,tGuessSd,pThreshold,beta,delta,gamma,grain,range);
            q{i}.normalizePdf = 1;
            q{i} = QuestUpdate(q{i}, tGuess,1);
            q{i} = QuestUpdate(q{i}, log10(0.005),0);
        end
    end
    if trial > learnTNum && trial <= learnTNum+connectTNum
        for i = 1:4
            q{i}.grain=0.02;
            q{i}.pThreshold = q{i}.pThreshold - (learnP-testP)/connectTNum;
            q{i} = QuestRecompute(q{i});
        end
    end
    
    % 【FIX #6】Proper message formatting with color calibration
    % Write trial start message
    Eyelink('Message', 'TRIALID %d', trial);
    % Format colors correctly (0-255 range)
    bgClr_scaled = round(bgClr/colorMaxVal*255);
    Eyelink('Message', '!V CLEAR %d %d %d', bgClr_scaled, bgClr_scaled, bgClr_scaled);
    Eyelink('Command', 'record_status_message "TRIAL %d/%d"', trial, trialNum);
    
    % Start recording
    % Start recording (stay in recording while sampling)
    Eyelink('SetOfflineMode');
    Eyelink('Command', 'clear_screen 0');
    Eyelink('StartRecording');
    WaitSecs(0.1);
    
    % Show fixation using unified function
    fixCenter = ut.Pol2Rect([rFix,results.oriF(trial)]).*[1,-1]+bgCenter;
    Eyelink('Message', 'FIX ON');
    if DoDriftCorrect
        EyelinkDoDriftCorrect(el, fixCenter(1), fixCenter(2), 1, 1); % last parameter = 1: allow recalibration
    else
        [startT, headDist] = show_fix(wpnt, fixCenter(1), fixCenter(2), fixTime, 0, winRect, MaxErr, ...
                                    'eyeTrackerType', 'EyeLink', 'el', el, ...
                                    'monWidth', monWidth, 'monHeight', monHeight);
    end
    
    ut = UT(monWidth, scWidth, headDist);
    tgCenter = ut.deg2pix([results.Xtarg(trial), results.Ytarg(trial)]);
    
    % Get Quest contrast
    if results.ClusterTags(trial)==0
        iq=4;
    else
        iq=results.ClusterTags(trial);
    end
    tTest=QuestQuantile(q{iq});
    results.tgContrast(trial) = 10.^tTest;
    
    % Generate and show stimulus
    stimulus = genStim(winRect, ut, results.bgContrast(trial), ...
        results.tgContrast(trial), tgCenter.*[1,-1]+bgCenter, GaborSF, GaborWidth, ...
        GaborOrient, bgWidth, results.seed(trial));
    stiTex = Screen('MakeTexture', wpnt, cat(3,stimulus,stimulus,stimulus)*255);
    
    Screen('Drawtexture',wpnt,stiTex);
    imgT = Screen('Flip',wpnt);
    Eyelink('Message', sprintf('STIM ON: F. trial-%.0f  ECC-%.0f  Ori-%.0f  tgContrast-%.3f bgContrast-%.3f', ...
        trial, results.ECC(trial), results.Orient(trial), results.tgContrast(trial), results.bgContrast(trial)));
    
    timeCost = 0;
    
    % Wait for response
    while 1
        checkend;
        [keyIsDown, secs, keyCode] = KbCheck;
        if keyCode(KbName(key1)) || timeCost>=maxTrialDur
            while 1
                checkend;
                [keyIsDown2, secs2, keyCode2] = KbCheck;
                if keyCode2(KbName(key2)) || timeCost>=maxTrialDur
                    % Get gaze data from EyeLink while recording
                    lastFixPix_ = [0, 0];
                    if Eyelink('NewFloatSampleAvailable') > 0
                        evt = Eyelink('NewestFloatSample');
                        eyeUsed = Eyelink('EyeAvailable');
                        if eyeUsed == el.BINOCULAR
                            eyeUsed = el.LEFT_EYE;
                        end
                        if eyeUsed == el.LEFT_EYE && evt.gx(1) ~= -32768 && evt.gy(1) ~= -32768
                            lastFixPix_ = [evt.gx(1), evt.gy(1)];
                        elseif eyeUsed == el.RIGHT_EYE && evt.gx(2) ~= -32768 && evt.gy(2) ~= -32768
                            lastFixPix_ = [evt.gx(2), evt.gy(2)];
                        else
                            if evt.gx(1) ~= -32768 && evt.gy(1) ~= -32768
                                lastFixPix_ = [evt.gx(1), evt.gy(1)];
                            elseif evt.gx(2) ~= -32768 && evt.gy(2) ~= -32768
                                lastFixPix_ = [evt.gx(2), evt.gy(2)];
                            end
                        end
                    end
                    FixNum = 1; % Simplified
                    
                    lastFixPix = lastFixPix_ - bgCenter;
                    lastFixPix(2) = -lastFixPix(2);
                    lastFixDeg = ut.pix2deg(lastFixPix);
                    break
                end
                timeCost = WaitSecs(0.01)-imgT;
            end
            if keyIsDown2 || timeCost>=maxTrialDur
                break
            end
        end
        timeCost = WaitSecs(0.01)-imgT;
    end
    
    % Feedback
    err = norm(lastFixDeg - [results.Xtarg(trial),results.Ytarg(trial)]);
    if err<MaxErr
        judgement = 1;
        pointColor = [0, 0, 255];
    else
        judgement = 0;
        pointColor = [255, 255, 0];
    end
    
    Screen('Drawtexture',wpnt, stiTex);
    Screen('DrawDots', wpnt, lastFixPix_, 14, pointColor, [], 3);
    Screen('Flip',wpnt);
    
    if scr>1
        tfirn = min(FixNum,firstn);
        earlyFixs = [lastFixPix_(1), lastFixPix_(2)]./(winRect(3)/win_main(3));
        textSize = 16;
        earlyFixs = earlyFixs-[textSize/4,textSize/2];
        Screen('TextSize', wpnt_main, textSize);
        if trial>1
            Screen('DrawTexture', wpnt_main, lastFrameTexture);
        end
        if judgement ==1
            tgColor = [0,0,255];
        else
            tgColor = [237,145,33];
        end
        Screen('DrawText', wpnt_main, '1', earlyFixs(1), earlyFixs(2), tgColor);
        Screen('DrawDots', wpnt_main, (tgCenter.*[1,-1]+bgCenter)./(winRect(3)/win_main(3)), 8, tgColor, [], 3);
        Screen('Flip',wpnt_main);
        lastScreenShot = Screen('GetImage', wpnt_main);
        lastScreenShot = cat(3, lastScreenShot, alphaChannel);
        lastFrameTexture = Screen('MakeTexture', wpnt_main, lastScreenShot);
    end
    
    WaitSecs(0.3);
    
    % Record data
    results.XeyeFix(trial) = lastFixDeg(1);
    results.YeyeFix(trial) = lastFixDeg(2);
    if keyCode(KbName(key1))
        results.key1RT(trial) = secs-imgT;
    end
    if keyCode2(KbName(key2))
        results.key2RT(trial)= secs2-imgT;
    end
    results.Err(trial) = err;
    results.judge(trial) = judgement;
    results.FixNum(trial) = FixNum;
    results.headDist(trial) = headDist;
    
    endT = Screen('Flip',wpnt);
    Screen('Close',stiTex);
    Eyelink('Message', sprintf('STIM OFF: F. FixNum-%.0f  Err-%.3f RT-%.3f  Recent10judge-%s ', ...
        FixNum, err, results.key1RT(trial), num2str(results.judge(max(trial-10,1):trial))));
    
    Eyelink('StopRecording');
    WaitSecs(rand(1)*0.4);
    
    % Rest break
    if mod(trial,blockSize)==0 && trial ~=trialNum
        quit = disp_rest(wpnt, bgCenter(1), bgCenter(2), restTime);
        if quit
            break
        end
    end
    
    % Update Quest
    if judgement || ~isnan(results.key2RT(trial))
        questJudge = 1;
    else
        questJudge = 0;
    end
    q{iq}=QuestUpdate(q{iq},tTest,questJudge);
end

%% STEP 6: Close EDF file and transfer to display PC
Eyelink('SetOfflineMode');
Eyelink('Command', 'clear_screen 0');
WaitSecs(0.5);
Eyelink('CloseFile');

% Transfer EDF file
try
    fprintf('Receiving data file ''%s''\n', edfFile);
    status = Eyelink('ReceiveFile');
    if status > 0
        fprintf('ReceiveFile status %d\n', status);
    end
    if exist([edfFile '.edf'], 'file')
        fprintf('Data file ''%s'' can be found in ''%s''\n', edfFile, pwd);
        % Rename with full experiment info
        newName = sprintf('./Data/Formal/%s_Sub%d_Ses%d_%s.edf', edfFile, subjID, session, DTstr);
        movefile([edfFile '.edf'], newName);
    end
catch
    fprintf('Problem receiving data file ''%s''\n', edfFile);
end

% Cleanup
PsychPortAudio('Close', pahandle);
PsychPortAudio('Close', pamaster);
Eyelink('Shutdown');