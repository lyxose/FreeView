% this demo code is part of Titta, a toolbox providing convenient access to
% eye tracking functionality using Tobii eye trackers
%
% Titta can be found at https://github.com/dcnieho/Titta. Check there for
% the latest version.
% When using Titta, please cite the following paper:
%
% Niehorster, D.C., Andersson, R. & Nystrom, M., (2020). Titta: A toolbox
% for creating Psychtoolbox and Psychopy experiments with Tobii eye
% trackers. Behavior Research Methods.
% doi: https://doi.org/10.3758/s13428-020-01358-8


clear all
sca
if ~exist('./Data', 'dir')
    mkdir('./Data');
end
if ~exist('./Data/Interrupted','dir')
    mkdir('./Data/Interrupted')
end
if ~exist('./Data/Threshold','dir')
    mkdir('./Data/Threshold')
end
DTstr = datestr(datetime, 'yyyymmddTHHMM');

%% Import functions
addpath(genpath('function_library'));
addpath('function_library_cus');
instFolder = './Instructions';
%% Parameters
DEBUGlevel              = 0;
saveRaw                 = true;    % ~200MB for 72 trials, 10min
fixClrs                 = [0 255];
bgClr                   = 127;
useAnimatedCalibration  = true;
doBimonocularCalibration= false;
% task parameters
fixTime                 = 1; % should be more than 0.1
key1 = 'space';
key2 = 'm'; % 
restTime = 15;   % second
blockSize = 96; % trials

[subjID, session, location, subjName, subjGender, subjAge, threshold, threDat] = InformationBox;

tobiiFreq = 250; % hz
tobiiMod = 'Tobii Pro Fusion';

bgWidth = 15;       % background width, in degree
GaborSF = 6;        % cycle per degree
GaborCyc = 2;       % target width (full width at half maxima, FWHM), n cycle
GaborWidth = GaborCyc/GaborSF; % target width in degree
GaborOrient = -45;   % Orientation of Garbor
difficulty = 1/2; 

Eccent = [2, 4, 6]; % Eccentricity of location
nOrient = 8;        % number of location orientation, start by positive x-axis;
Orient = linspace(0,360*(nOrient-1)/nOrient,nOrient); % angle of location to center, counter-clockwise from positive x-axis;
tgContrast = threshold/difficulty; % contrast of target center (Gabor)
bgContrast = 0.2;             % maximum contrast of background texture
repeatTimes = 20;              % repeat n times for each condition

%% Create the result matrix
% the 1nd  column denotes to the Eccent of target;
% the 2rd  column denotes to the angle counter-clockwise from positive x-axis;
% the 3th  column denotes to contrast of target (Gabor);
% the 4th  column denotes to background contrast (1/f noise);
% the 5th  column denotes to repeat id;
% the 6th  column denotes to reaction time of press the first key;
% the 7th  column denotes to reaction time of press the second key;
% the 8th  column denotes to horizontal coordinate of eye fixation;
% the 9th  column denotes to horizontal coordinate of target center;
% the 10th column denotes to vertical coordinate of eye fixation;
% the 11th column denotes to vertical coordinate of target center;
% the 12th column denotes to response judgement;

[~,results] = expRun.generateTrialList('ECC',Eccent,'Orient',Orient, ...
    'tgContrast',tgContrast,'bgContrast',bgContrast,'repeatID',1:repeatTimes, ...
    'key1RT',nan,'key2RT',nan,'XeyeFix',nan,'Xtarg',nan, ...
    'YeyeFix',nan,'Ytarg',nan,'FixNum',nan,'Err',nan,'judge',nan, ...
    'headDist',nan);

trialNum = height(results);
results.seed = randi(10*trialNum,[trialNum,1]);
% shuffle
trialOrder = Shuffle(1:trialNum);
results = results(trialOrder,:);

% ATTENTION!! 
% Take the centre of the screen as the coordinate origin, right and up as 
% the positive direction
results.Xtarg = results.ECC .* cosd(results.Orient);
results.Ytarg = results.ECC .* sind(results.Orient);
taskSpace= results(:,{'Xtarg','Ytarg'});
[~, ia, ic] = unique(taskSpace,"rows");
taskSpace = taskSpace(ia, :);
judger = spatialJudge(taskSpace);

%% Initiation of stimulation

% fix the hardware bugs
[ keyIsDown, ~, keyCode ] = KbCheck;
keyCode = find(keyCode, 1);
if keyIsDown
    ignoreKey=keyCode;
    DisableKeysForKbCheck(ignoreKey);
end
Screen('Preference', 'SyncTestSettings', 0.002); % the systems are a little noisy, give the test a little more leeway
KbName('UnifyKeyNames');

% PTB parameters
scr = max(Screen('Screens'));
screens      = Screen('Screens');
screenNumber = max(screens);
scWidth = Screen('Resolution',screenNumber).width; % width of screen resolution, in pixel
scHeight = Screen('Resolution',screenNumber).height; % width of screen resolution, in pixel
bgCenter = round([scWidth/2, scHeight/2]);

% Titta initiation
try
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
    
    % init
    EThndl          = Titta(settings);
    % EThndl          = EThndl.setDummyMode();    % just for internal testing, enabling dummy mode for this readme makes little sense as a demo
    EThndl.init();

    if DEBUGlevel>1
        % make screen partially transparent on OSX and windows vista or
        % higher, so we can debug.
        PsychDebugWindowConfiguration;
    end
    if DEBUGlevel
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
    
    oper = showInstruc(wpnt, 'welcome', instFolder, 'space', 'BackSpace');

    % do calibration
    try
        ListenChar(-1);
    catch ME
        % old PTBs don't have mode -1, use 2 instead which also supresses
        % keypresses from leaking through to matlab
        ListenChar(2);
    end
    if DEBUGlevel==0
        if doBimonocularCalibration
            % do sequential monocular calibrations for the two eyes
            settings                = EThndl.getOptions();
            settings.calibrateEye   = 'left';
            settings.UI.button.setup.cal.string = 'calibrate left eye (<i>spacebar<i>)';
            str = settings.UI.button.val.continue.string;
            settings.UI.button.val.continue.string = 'calibrate other eye (<i>spacebar<i>)';
            EThndl.setOptions(settings);
            tobii.calVal{1}         = EThndl.calibrate(wpnt,1);
            if ~tobii.calVal{1}.wasSkipped
                settings.calibrateEye   = 'right';
                settings.UI.button.setup.cal.string = 'calibrate right eye (<i>spacebar<i>)';
                settings.UI.button.val.continue.string = str;
                EThndl.setOptions(settings);
                tobii.calVal{2}         = EThndl.calibrate(wpnt,2);
            end
        else
            % do binocular calibration
            tobii.calVal{1} = EThndl.calibrate(wpnt,3,threDat.calibration{1});
        end
    end
    ListenChar(0);
    showInstruc(wpnt, 'T2P', instFolder, 'space', 'BackSpace');
    % later:
    EThndl.buffer.start('gaze');
    WaitSecs(.8);   % wait for eye tracker to start and gaze to be picked up
    monWidth = EThndl.geom.displayArea.width/10;    % in cm
    monHeight = EThndl.geom.displayArea.height/10;  % in cm
    
    passed = 0;
    CtrGrad = 10.^ linspace(log10(0.75), log10(0.2), 10); % contrast gradient in pretrials, from 0.6 to 0.12
    if length(bgContrast)~=1
        pre_bgContrast = mean(bgContrast);
    else
        pre_bgContrast = bgContrast;
    end
    while ~passed
        for pretrial = 1:10
            [startT, headDist] = show_fix(wpnt, fixTime, fixClrs, winRect, true, tobiiFreq, EThndl, monWidth, monHeight);
            EThndl.sendMessage('FIX ON Pre',startT);
            ut = UT(monWidth, scWidth, headDist);
            tgCenter_= rand([1,2]).*[7,360];  % random polor coordinate
            tgCenter = ut.Pol2Rect(tgCenter_); % central rectangular coordinate in pixel
            % prep stimuli 
            stimulus = genStim(winRect, ut, pre_bgContrast, CtrGrad(pretrial), ...
                tgCenter, GaborSF, GaborWidth, GaborOrient, bgWidth, randi(1000));
            stiTex = Screen('MakeTexture', wpnt, cat(3,stimulus,stimulus,stimulus)*255);
            % show stimulus        
            Screen('Drawtexture',wpnt,stiTex);
            imgT = Screen('Flip',wpnt,startT+fixTime-1/hz/2);                   % bit of slack to make sure requested presentation time can be achieved
            SysImgT = EThndl.getTimeAsSystemTime(imgT);
            EThndl.sendMessage(sprintf('STIM ON: P. ECC-%.0f  Ori-%.0f  tgContrast-%.3f  bgContrast-%.2f', ...
                   tgCenter_(1), tgCenter_(2), CtrGrad(pretrial), pre_bgContrast), imgT);

            % Check key response
            while 1 
                checkend;
                [keyIsDown, secs, keyCode] = KbCheck;
                if keyCode(KbName(key1))
                    while 1
                        checkend;
                        [keyIsDown2, secs2, keyCode2] = KbCheck;
                        if keyCode2(KbName(key2))
                            gazeData = EThndl.buffer.peekTimeRange('gaze',SysImgT);
                            [lastFixPix_,tdat] = getLastFix(gazeData, monWidth, monHeight, tobiiFreq, headDist, winRect(3:4), SysImgT, 15, 60);       
                            FixNum = length(tdat.fix.dur);
                            % lastFixPix is in centimeter, [0,0] for upperleft corner
                            lastFixPix = lastFixPix_ - bgCenter; % Change the coordinate origin
                            lastFixPix(2) = -lastFixPix(2);     % Change the y-axis direction
                            lastFixDeg = ut.pix2deg(lastFixPix);
                            break
                        end
                        WaitSecs(0.05); % reduce the sampling rate to reduce the pressure of CPU
                    end
                    if keyIsDown2
                        break 
                    end
                    WaitSecs(0.05);
                end
            end

            % Color feedback
            err = norm(lastFixDeg - ut.pix2deg(tgCenter));
            if err<1
                pointColor = [0, 0, 255]; % blue indicates right fixation
            else
                pointColor = [255, 255, 0]; % yellow indicates wrong fixation
            end
            Screen('Drawtexture',wpnt, stiTex);
            Screen('DrawDots', wpnt, lastFixPix_, 63, pointColor);
            fbT = Screen('Flip',wpnt);
            WaitSecs(0.3);

            endT = Screen('Flip',wpnt); % finiT+1-1/hz/2 
            Screen('Close',stiTex);
            EThndl.sendMessage(sprintf('STIM OFF: P. trial-%.0f  FixNum-%.0f  Err-%.3f  judge-%.0f', ...
                                                  pretrial,      FixNum,      err,      err<1),endT);
            WaitSecs(rand(1)*0.4); % to prevent any long term rhythm.

            if err>=1  % start the next 10 trials, untill all correct
                break
            end
            if pretrial ==10
                passed = 1;
            end
        end
        oper = showInstruc(wpnt, 'T2Check', instFolder, 'space', 'BackSpace');
        if oper==-1
            passed = 0;
            tobii.calVal{1} = EThndl.calibrate(wpnt,3,tobii.calVal{1}); 
            EThndl.buffer.start('gaze');
            WaitSecs(.8);   % wait for eye tracker to start and gaze to be picked up
        end
    end
    
    oper = showInstruc(wpnt, 'T2F', instFolder, 'space', 'BackSpace');
    
    for trial = 1:trialNum
        % First draw a fixation point
        [startT, headDist] = show_fix(wpnt, fixTime, fixClrs, winRect, true, tobiiFreq, EThndl, monWidth, monHeight);

        EThndl.sendMessage('FIX ON',startT);

        ut = UT(monWidth, scWidth, headDist);
        tgCenter = ut.deg2pix([results.Xtarg(trial), results.Ytarg(trial)]);
        % prep stimuli 
        stimulus = genStim(winRect, ut, results.bgContrast(trial), ...
            results.tgContrast(trial), tgCenter, GaborSF, GaborWidth, ...
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
        EThndl.sendMessage(sprintf('STIM ON: F_%0.f_%0.f ECC-%.0f  Ori-%.0f  tgContrast-%.3f  bgContrast-%.2f', ...
               results.ECC(trial), results.Orient(trial), results.ECC(trial), results.Orient(trial), results.tgContrast(trial), results.bgContrast(trial)),imgT);
        while 1 
            checkend;
            [keyIsDown, secs, keyCode] = KbCheck;
            if keyCode(KbName(key1))
                while 1
                    checkend;
                    [keyIsDown2, secs2, keyCode2] = KbCheck;
                    if keyCode2(KbName(key2))
                        gazeData = EThndl.buffer.peekTimeRange('gaze',SysImgT);
                        [lastFixPix_,tdat] = getLastFix(gazeData, monWidth, monHeight, tobiiFreq, headDist, winRect(3:4), SysImgT, 15, 60);       
                        FixNum = length(tdat.fix.dur);
                        % lastFixPix is in centimeter, [0,0] for upperleft corner
                        lastFixPix = lastFixPix_ - bgCenter; % Change the coordinate origin
                        lastFixPix(2) = -lastFixPix(2);     % Change the y-axis direction
                        lastFixDeg = ut.pix2deg(lastFixPix);
                        break
                    end
                    WaitSecs(0.05); % reduce the sampling rate to reduce the pressure of CPU
                end
                if keyIsDown2
                    break 
                end
                WaitSecs(0.05);
            end
        end
        

        % action feedback
        [judgement, err]= judger.judge(lastFixDeg, [results.Xtarg(trial),results.Ytarg(trial)]);
        if judgement == 1
            pointColor = [0, 0, 255]; % blue indicates right fixation
        else
            pointColor = [255, 255, 0]; % yellow indicates wrong fixation
        end
        Screen('Drawtexture',wpnt, stiTex);
        Screen('DrawDots', wpnt, lastFixPix_, 63, pointColor);
        fbT = Screen('Flip',wpnt);
        WaitSecs(0.3); 

        % record data
        results.XeyeFix(trial) = lastFixDeg(1);
        results.YeyeFix(trial) = lastFixDeg(2);
        results.key1RT(trial) = secs-imgT;
        results.key2RT(trial)= secs2-imgT;
        results.Err(trial) = err;
        results.judge(trial) = judgement;
        results.FixNum(trial) = FixNum;
        results.headDist(trial) = headDist;
        
        endT = Screen('Flip',wpnt); % finiT+1-1/hz/2 
        Screen('Close',stiTex);
        EThndl.sendMessage(sprintf('STIM OFF: F_%.0f_%.0f trial-%.0f  FixNum-%.0f  Err-%.3f  judge-%.0f', ...
               results.ECC(trial), results.Orient(trial), trial,      FixNum,      err,      judgement),endT);
        WaitSecs(rand(1)*0.4); % to prevent any long term rhythm.
        % take a break every blockSize trials
        if mod(trial,blockSize)==0 && trial ~=trialNum
           quit = disp_rest(wpnt, bgCenter(1), bgCenter(2), restTime);
           if quit
               break
           end
        end
    end
    dat = EThndl.collectSessionData();
    dat.expt.restTime    = restTime;
    dat.expt.blockSize   = blockSize;
    dat.expt.winRect     = winRect;
    dat.expt.bgWidth     = bgWidth;     
    dat.expt.GaborSF     = GaborSF;     
    dat.expt.GaborCyc    = GaborCyc;      
    dat.expt.GaborWidth  = GaborWidth;        
    dat.expt.GaborOrient = GaborOrient;  
    dat.expt.difficulty  = difficulty;

    EThndl.saveData(dat, fullfile(cd,sprintf('./Data/Dat_Sub%.0f_Ses%.0f',subjID, session)), true);

    save(sprintf('./Data/Result_Sub%.0f_Ses%.0f_%s_%s_%s',subjID, session, location, subjName, DTstr),"results")
    if saveRaw
        save(sprintf('./Data/EXP_Sub%.0f_Ses%.0f_%s_%s_%s',subjID, session, location, subjName, DTstr))
    end
    % if you want to (also) save the data to Apache Parquet and json files
    % that can easily be read in Python (Apache Parquet files are supported
    % by Pandas), use:
    % EThndl.saveDataToParquet(dat, fullfile(cd,'t'), true);
    % All gaze data columns and messages can be dumped to tsv files using:
    % EThndl.saveGazeDataToTSV(dat, fullfile(cd,'t'), true);
    
    % shut down
    EThndl.deInit();
    showInstruc(wpnt, 'End', instFolder, 'space', 'BackSpace')
catch me
    dat = EThndl.collectSessionData();
    dat.expt.restTime    = restTime;
    dat.expt.blockSize   = blockSize;
    dat.expt.winRect     = winRect;
    dat.expt.bgWidth     = bgWidth;     
    dat.expt.GaborSF     = GaborSF;     
    dat.expt.GaborCyc    = GaborCyc;      
    dat.expt.GaborWidth  = GaborWidth;        
    dat.expt.GaborOrient = GaborOrient;  
    dat.expt.difficulty  = difficulty;

    save(sprintf('./Data/Interrupted/EXPINT_Sub%.0f_Ses%.0f_%s_%s_%s',subjID, session, location, subjName, DTstr))
    sca
    ListenChar(0);
    rethrow(me)
end
sca

% 

% to visualize current data
% img.data = cat(3,stimulus,stimulus,stimulus);
% img.x = 1:scWidth;
% img.y = 1:scHeight;
% drawFix(tdat,tdat.fix,[tdat.I2MCopt.xres tdat.I2MCopt.yres],img,[tdat.I2MCopt.missingx tdat.I2MCopt.missingy],sprintf('test'));

% record x seconds of data, then clear screen. Indicate stimulus
% removed, clean up
% finiT = Screen('Flip',wpnt);


