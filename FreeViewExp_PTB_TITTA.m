


% 
% if threshold == 0
%     grain = 0.01;  % stepsize
%     tGuess = log10(0.1);  % prior threshold estimate of central position
% else
%     grain = 0.002;
%     tGuess = log10(threshold);  % prior threshold estimate of central position
% end    
% beta=1.5;delta=0.01;gamma=0.25; % Be careful!
% range = 3;  % the intensity difference between the largest and smallest intensity,  log10([0.01 1]) =-2   0
% % create the initial 'q'
% q = QuestCreate(tGuess,tGuessSd,pThreshold,beta,delta,gamma,grain,range);
% q.normalizePdf = 1;
% q = QuestUpdate(q, log10(0.5),1);
% q = QuestUpdate(q, log10(0.005),0);

if DEBUGlevel>1
    % make screen partially transparent on OSX and windows vista or
    % higher, so we can debug.
    PsychDebugWindowConfiguration;
end
if DEBUGlevel
    Screen('Preference', 'SkipSyncTests', 1);
    % Be pretty verbose about information and hints to optimize your code and system.
    Screen('Preference', 'Verbosity', 4);
else
    % Only output critical errors and warnings.
    Screen('Preference', 'Verbosity', 2);
end

% Open PTB window
[wpnt,winRect] = PsychImaging('OpenWindow', scr, bgClr, [], [], [], [], 4);
hz=Screen('NominalFrameRate', wpnt);
Priority(1);
Screen('BlendFunction', wpnt, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
Screen('Preference', 'TextAlphaBlending', 1);
Screen('Preference', 'TextAntiAliasing', 2);
% This preference setting selects the high quality text renderer on
% each operating system: It is not really needed, as the high quality
% renderer is the default on all operating systems, so this is more of
% a "better safe than sorry" setting.
Screen('Preference', 'TextRenderer', 1);
KbName('UnifyKeyNames');    % for correct operation of the setup/calibration interface, calling this is required
% if task==1
showInstruc(wpnt, 'welcome', instFolder, 'space', 'BackSpace');
showInstruc(wpnt, 'Calib', instFolder, 'space', 'BackSpace');
% end
% do calibration
try
    ListenChar(-1);
catch ME
    % old PTBs don't have mode -1, use 2 instead which also supresses
    % keypresses from leaking through to matlab
    ListenChar(2);
end

% calibrate and practice
reCali = true;
while reCali
% Titta initiation
% get setup struct (can edit that of course):
settings = Titta.getDefaults(tobiiMod);
settings.freq = tobiiFreq;
% request some debug output to command window, can skip for normal use
settings.debugMode      = true;
% customize colors of setup and calibration interface (colors of
% everything can be set, so there is a lot here).
% 1. setup screen
settings.UI.setup.bgColor       = bgClr;
settings.UI.setup.instruct.color= fixClrs(1);
settings.UI.setup.fixBackColor  = fixClrs(1);
settings.UI.setup.fixFrontColor = fixClrs(2);
% 2. calibration display
if useAnimatedCalibration
    % custom calibration drawer
    calViz                      = AnimatedCalibrationDisplay();
    settings.cal.drawFunction   = @calViz.doDraw;
    calViz.bgColor              = bgClr;
    calViz.fixBackColor         = fixClrs(1);
    calViz.fixFrontColor        = fixClrs(2);
else
    % set color of built-in fixation points
    settings.cal.bgColor        = bgClr;
    settings.cal.fixBackColor   = fixClrs(1);
    settings.cal.fixFrontColor  = fixClrs(2);
end
% callback function for completion of each calibration point
settings.cal.pointNotifyFunction = @demoCalCompletionFun;
% 3. validation result screen
settings.UI.val.bgColor                 = bgClr;
settings.UI.val.avg.text.color          = fixClrs(1);
settings.UI.val.fixBackColor            = fixClrs(1);
settings.UI.val.fixFrontColor           = fixClrs(2);
settings.UI.val.onlineGaze.fixBackColor = fixClrs(1);
settings.UI.val.onlineGaze.fixFrontColor= fixClrs(2);

% init tobbi
EThndl          = Titta(settings);
% EThndl          = EThndl.setDummyMode();    % just for internal testing, enabling dummy mode for this readme makes little sense as a demo
EThndl.init();
if DEBUGlevel==0 
    if doBimonocularCalibration 
        % do sequential monocular calibrations for the two eyes
        settings                = EThndl.getOptions();
        settings.calibrateEye   = 'left';
        settings.UI.button.setup.cal.string = 'calibrate left eye (<i>spacebar<i>)';
        str = settings.UI.button.val.continue.string;
        settings.UI.button.val.continue.string = 'calibrate other eye (<i>spacebar<i>)';
        EThndl.setOptions(settings);
        try
            tobii.calVal{1}         = EThndl.calibrate(wpnt,1,tobii.calVal{1});
        catch
            tobii.calVal{1}         = EThndl.calibrate(wpnt,1);            
        end
        if ~tobii.calVal{1}.wasSkipped
            settings.calibrateEye   = 'right';
            settings.UI.button.setup.cal.string = 'calibrate right eye (<i>spacebar<i>)';
            settings.UI.button.val.continue.string = str;
            EThndl.setOptions(settings);
            try 
                tobii.calVal{2}         = EThndl.calibrate(wpnt,2,tobii.calVal{2});
            catch
                tobii.calVal{2}         = EThndl.calibrate(wpnt,2);                
            end
        end
    else
        % do binocular calibration
        try
            tobii.calVal{1} = EThndl.calibrate(wpnt,3,tobii.calVal{1});%,threDat.calibration{1});            
        catch
            tobii.calVal{1} = EThndl.calibrate(wpnt,3);%,threDat.calibration{1});            
        end
    end
end

ListenChar(0);

EThndl.buffer.start('gaze');
WaitSecs(.8);   % wait for eye tracker to start and gaze to be picked up
monWidth = EThndl.geom.displayArea.width/10;    % in cm
monHeight = EThndl.geom.displayArea.height/10;  % in cm
headDist = default_distance;

if length(bgContrast)~=1
    pre_bgContrast = mean(bgContrast);
else
    pre_bgContrast = bgContrast;
end

if threshold==0
    showInstruc(wpnt, 'Practice', instFolder, 'space', 'BackSpace');
    CtrGrad = 10.^ linspace(log10(0.6), log10(0.2), 10); % contrast gradient in pretrials, from 0.6 to 0.15
    passed = 0;
else 
    passed = 1;
    ut = UT(monWidth, scWidth, headDist);
end
reCali = false;
while ~passed && ~reCali
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
        [startT, headDist] = show_fix(wpnt, fixCenter(1), fixCenter(2), fixTime, fixClrs, winRect, true, tobiiFreq, EThndl, monWidth, monHeight);
        EThndl.sendMessage('FIX ON Pre',startT);
%         tgCenter_= [sqrt(R_min^2 + (R_max^2 - R_min^2) * rand()), rand() * 360];  % random polor coordinate
        tgCenter = ut.deg2pix([preresults.Xtarg(pretrial), preresults.Ytarg(pretrial)]);
%         tgCenter = ut.Pol2Rect(tgCenter_); % central rectangular coordinate in pixel
        % prep stimuli 
        stimulus = genStim(winRect, ut, pre_bgContrast, CtrGrad(pretrial), ...
            tgCenter.*[1,-1]+bgCenter, GaborSF, GaborWidth, GaborOrient, bgWidth, preresults.seed(pretrial));
        stiTex = Screen('MakeTexture', wpnt, cat(3,stimulus,stimulus,stimulus)*255);
        % show stimulus        
        Screen('Drawtexture',wpnt,stiTex);
        imgT = Screen('Flip',wpnt,startT+fixTime-1/hz/2);                   % bit of slack to make sure requested presentation time can be achieved
        SysImgT = EThndl.getTimeAsSystemTime(imgT);
        EThndl.sendMessage(sprintf('STIM ON: P. P. trial-%.0f  ECC-%.2f  Ori-%.0f  tgContrast-%.3f  bgContrast-%.2f', ...
                                                pretrial,      preresults.ECC(pretrial), preresults.Orient(pretrial), CtrGrad(pretrial), pre_bgContrast), imgT);
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
                        gazeData = EThndl.buffer.peekTimeRange('gaze',SysImgT);
                        [lastFixPix_,tdat] = getLastFix(gazeData, monWidth, monHeight, tobiiFreq, headDist, winRect(3:4), SysImgT, 15, 60);       
                        FixNum = length(tdat.fix.dur);
                        % lastFixPix is in centimeter, [0,0] for upperleft corner
                        lastFixPix = lastFixPix_ - bgCenter; % Change the coordinate origin
                        lastFixPix(2) = -lastFixPix(2);     % Change the y-axis direction
                        lastFixDeg = ut.pix2deg(lastFixPix);
                        break
                    end
                    timeCost = WaitSecs(0.01)-imgT; % reduce the sampling rate to reduce the pressure of CPU
                end
                if keyIsDown2 || timeCost>=maxTrialDur
                    break 
                end
            end
            timeCost = WaitSecs(0.01)-imgT;
        end

        % Color feedback
        err = norm(lastFixDeg - ut.pix2deg(tgCenter));
        if err<MaxErr
            pointColor = [0, 0, 255]; % blue indicates right fixation
        else
            pointColor = [255, 255, 0]; % yellow indicates wrong fixation
        end
        Screen('Drawtexture',wpnt, stiTex);
        Screen('DrawDots', wpnt, lastFixPix_, 14, pointColor, [], 3);
        fbT = Screen('Flip',wpnt);
        WaitSecs(0.3);

        endT = Screen('Flip',wpnt); % finiT+1-1/hz/2 
        Screen('Close',stiTex);
        EThndl.sendMessage(sprintf('STIM OFF: FixNum-%.0f  Err-%.3f  judge-%.0f RT-%.3f', ...
                                              FixNum,      err,      err<MaxErr, secs-imgT),endT);
        WaitSecs(rand(1)*0.4); % to prevent any long term rhythm.

        if err>=MaxErr  % start the next 10 trials, untill all correct
            break
        end
        if pretrial ==10
            passed = 1;
        end
    end
    oper = showInstruc(wpnt, 'Check', instFolder, 'space', 'BackSpace');
    if oper==-1
        reCali = true; % re-calibrate
        EThndl.deInit();
    else 
        reCali = false; % passed
    end
end
end
WaitSecs(0.5);
% if threshold==0
%     oper = showInstruc(wpnt, 'T1F', instFolder, 'space', 'BackSpace');
% else
    oper = showInstruc(wpnt, 'T2F', instFolder, 'space', 'BackSpace');
% end    
for trial = 1:trialNum
    if trial ==1  % learning stage
        % Prepare parameters for Quest
        pThreshold = 0.85;
        %Threshold "t" is measured on an abstract "intensity" scale, which usually corresponds to log10 contrast.
        tGuessSd=1;  % the standard deviation you assign to that guess. Be generous!
        tGuess = log10(0.25);
        grain = 0.01;
        beta=1.5;delta=0.01;gamma=0.25; % Be careful!
        range = 3;  % the intensity difference between the largest and smallest intensity,  log10([0.01 1]) =-2   0
        q = QuestCreate(tGuess,tGuessSd,pThreshold,beta,delta,gamma,grain,range);
        q.normalizePdf = 1;
        q = QuestUpdate(q, log10(0.5),1);
        q = QuestUpdate(q, log10(0.005),0);
    end
    if trial == learnTNum
        % Prepare parameters for Quest
        pThreshold = 0.4;
        %Threshold "t" is measured on an abstract "intensity" scale, which usually corresponds to log10 contrast.
        tGuessSd=1;  % the standard deviation you assign to that guess. Be generous!
        tGuess = log10(results.tgContrast(trial-1));
        grain = 0.001;
        beta=1.5;delta=0.01;gamma=0.25; % Be careful!
        range = 3;  % the intensity difference between the largest and smallest intensity,  log10([0.01 1]) =-2   0
        q = QuestCreate(tGuess,tGuessSd,pThreshold,beta,delta,gamma,grain,range);
        q.normalizePdf = 1;
        for ii = 1:10
        q = QuestUpdate(q, tGuess,1);
        q = QuestUpdate(q, log10(0.005),0);
        end
    end
    % First draw a fixation point
    fixCenter = ut.Pol2Rect([rFix,results.oriF(trial)]).*[1,-1]+bgCenter;
    [startT, headDist] = show_fix(wpnt, fixCenter(1), fixCenter(2), fixTime, fixClrs, winRect, true, tobiiFreq, EThndl, monWidth, monHeight);

    EThndl.sendMessage('FIX ON',startT);

    ut = UT(monWidth, scWidth, headDist);
    tgCenter = ut.deg2pix([results.Xtarg(trial), results.Ytarg(trial)]);
    % prep stimuli 
    % Get recommended Contrast level by Pelli (1987)
    tTest=QuestQuantile(q);
    results.tgContrast(trial) = 10.^tTest; % 

    stimulus = genStim(winRect, ut, results.bgContrast(trial), ...
        results.tgContrast(trial), tgCenter.*[1,-1]+bgCenter, GaborSF, GaborWidth, ...
        GaborOrient, bgWidth, results.seed(trial));
    stiTex = Screen('MakeTexture', wpnt, cat(3,stimulus,stimulus,stimulus)*255);
                    
    % show on screen and log when it was shown in eye-tracker time.
    % NB: by setting a deadline for the flip, we ensure that the previous
    % screen (fixation point) stays visible for the indic m m m m m m m m ated amount of
    % time. See PsychToolbox demos for further elaboration on this way of
    % timing your script.
    Screen('Drawtexture',wpnt,stiTex);
    imgT = Screen('Flip',wpnt,startT+fixTime-1/hz/2);                   % bit of slack to make sure requested presentation time can be achieved
    SysImgT = EThndl.getTimeAsSystemTime(imgT);
    EThndl.sendMessage(sprintf('STIM ON: F. trial-%.0f  ECC-%.0f  Ori-%.0f  tgContrast-%.3f bgContrast-%.3f', ...
                                            trial,      results.ECC(trial), results.Orient(trial), results.tgContrast(trial), results.bgContrast(trial)),imgT);
    timeCost = 0;

    while 1 
        checkend;
        [keyIsDown, secs, keyCode] = KbCheck;
        if keyCode(KbName(key1)) || timeCost>=maxTrialDur
            while 1
                checkend;
                [keyIsDown2, secs2, keyCode2] = KbCheck;
                if keyCode2(KbName(key2)) || timeCost>=maxTrialDur
                    gazeData = EThndl.buffer.peekTimeRange('gaze',SysImgT);
                    [lastFixPix_,tdat] = getLastFix(gazeData, monWidth, monHeight, tobiiFreq, headDist, winRect(3:4), SysImgT, 15, 60);       
                    FixNum = length(tdat.fix.dur);
                    % lastFixPix is in centimeter, [0,0] for upperleft corner
                    lastFixPix = lastFixPix_ - bgCenter; % Change the coordinate origin
                    lastFixPix(2) = -lastFixPix(2);     % Change the y-axis direction
                    lastFixDeg = ut.pix2deg(lastFixPix);
                    break
                end
                timeCost = WaitSecs(0.01)-imgT; % reduce the sampling rate to reduce the pressure of CPU
            end
            if keyIsDown2 || timeCost>=maxTrialDur
                break 
            end
        end
        timeCost = WaitSecs(0.01)-imgT;
    end
    

    % action feedback
    err = norm(lastFixDeg - [results.Xtarg(trial),results.Ytarg(trial)]);
    if err<MaxErr
        judgement = 1;
        pointColor = [0, 0, 255]; % blue indicates right fixation
    else
        judgement = 0;
        pointColor = [255, 255, 0]; % yellow indicates wrong fixation
    end
    Screen('Drawtexture',wpnt, stiTex);
    Screen('DrawDots', wpnt, lastFixPix_, 14, pointColor, [], 3);
    fbT = Screen('Flip',wpnt);
    WaitSecs(0.3); 

    % record data
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
    
    endT = Screen('Flip',wpnt); % finiT+1-1/hz/2 
    Screen('Close',stiTex);
    EThndl.sendMessage(sprintf('STIM OFF: F. FixNum-%.0f  Err-%.3f RT-%.3f  Recent10judge-%s ', ...
                                             FixNum,      err,     results.key1RT(trial),  num2str(results.judge(max(trial-10,1):trial))),endT);
    WaitSecs(rand(1)*0.4); % to prevent any long term rhythm.
    % take a break every blockSize trials
    if mod(trial,blockSize)==0 && trial ~=trialNum
       quit = disp_rest(wpnt, bgCenter(1), bgCenter(2), restTime);
       if quit
           break
       end
    end
    q=QuestUpdate(q,tTest,judgement);  % Add the new datum (actual test intensity and observer response) to the database.
end