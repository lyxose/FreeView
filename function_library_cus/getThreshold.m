% Get the threshold by Quest
%
function [t,dat] = getThreshold(trialNum, subjID, session, location, subjName)
% clear all                                                                                           
% sca
% trialNum=60;
    DTstr = datestr(datetime, 'yyyymmddTHHMM');
    instFolder = './Instructions';

    DEBUGlevel     = 0;
    saveRaw        = true;    % ~200MB for 72 trials, 10min
    fixClrs        = [0 255];
    bgClr          = 127;
    % task parameters
    fixTime        = 1;    % second
    key1           = ',<';
    key2           = '.>'; 
    key12          = {key1,key2};
    stimT          = 0.25; % duration of showing stimulus, second
    interT         = 0.5;  % duration between two stimulus, second
    
    bgWidth = 15;       % background width, in degree
    GaborSF = 6;        % cycle per degree
    GaborCyc = 2;       % target width (full width at half maxima, FWHM), n cycle
    GaborWidth = GaborCyc/GaborSF; % target width in degree
    GaborOrient = -45;  % Orientation of Garbor
    
    MaxEcc = 1;         % Eccentricity range of random dots
    bgContrast = 0.2;   % maximum contrast of background texture
    %% Titta
    if isempty(TittaMex().findAllEyeTrackers())
        withTobii=false;
    else
        withTobii=true;        
    end
    if withTobii
        tobiiFreq = 250; % hz
        tobiiMod = 'Tobii Pro Fusion';
        useAnimatedCalibration  = true;
        doBimonocularCalibration= false;
    else
        tobiiFreq = [];
        monWidth = input('Screen width, in cm: (54.47853 for Alien of 348-2 room)');
        headDist = input('Enter the distance between head and screen, in cm: ');
        dat = struct();
    end
    
    %% Import functions
    addpath(genpath('function_library'));
    addpath('function_library_cus');
    
    %% Create the result matrix
    % the 1nd  column denotes to trial number
    % the 2rd  column denotes to the 1st or 2nd interval with target
    % the 3th  column denotes to background contrast
    % the 4th  column denotes to target contrast
    % the 5th  column denotes to reaction time
    % the 6th  column denotes to X coordinates of target
    % the 7th  column denotes to Y coordinates of target
    % the 8th  column denotes to judgement of decision
    [~,results] = expRun.generateTrialList('repeat',1:trialNum/2,'target',1:2, ...
                    'bgContrast',bgContrast,'tgContrast',nan,'Xtarg',nan, ...
                    'Ytarg',nan,'RT',nan,'judge',nan, ...
                    'seed1',randi(10*trialNum,[trialNum,1]), ...
                    'seed2',randi(10*trialNum,[trialNum,1]));
    trialOrder = Shuffle(1:trialNum);
    results = results(trialOrder,:);
    
    % Prepare parameters for Quest
    pThreshold = 0.7;
    tGuess = log10(0.05);  % prior threshold estimate of central position
    %Threshold "t" is measured on an abstract "intensity" scale, which usually corresponds to log10 contrast.
    tGuessSd=1;  % the standard deviation you assign to that guess. Be generous!
    beta=3.5;delta=0.01;gamma=0.25; % Be careful!
    grain = 0.01;  % stepsize
    range = 3;  % the intensity difference between the largest and smallest intensity,  log10([0.01 1]) =-2   0
    wrongRight={'wrong','right'};
    % create the initial 'q'
    q = QuestCreate(tGuess,tGuessSd,pThreshold,beta,delta,gamma,grain,range);
    q.normalizePdf = 1;
    q = QuestUpdate(q, log10(1),1);
    q = QuestUpdate(q, log10(0.01),0);
    
    % ATTENTION!! 
    % Take the centre of the screen as the coordinate origin, right and up as 
    % the positive direction
    es = rand(trialNum,1)*MaxEcc; % Random Eccentricities (uniform distribution)
    os = rand(trialNum,1)*2*pi;   % Random Orientations (uniform distribution)
    results.Xtarg = es.*cos(os);
    results.Ytarg = es.*sin(os);
    
    %% Initiation of stimulation
    
    % fix the hardware bugs
    [ keyIsDown, ~, keyCode ] = KbCheck;
    keyCode = find(keyCode, 1);
    if keyIsDown
        ignoreKey=keyCode;
        DisableKeysForKbCheck(ignoreKey);
    end
    Screen('Preference', 'SyncTestSettings', 0.002);    % the systems are a little noisy, give the test a little more leeway

    % PTB parameters
    scr = max(Screen('Screens'));
    screens      = Screen('Screens');
    screenNumber = max(screens);
    scWidth = Screen('Resolution',screenNumber).width; % width of screen resolution, in pixel
    scHeight = Screen('Resolution',screenNumber).height; % width of screen resolution, in pixel
    bgCenter = round([scWidth/2, scHeight/2]);
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
    try
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
        showInstruc(wpnt, 'welcome', instFolder, 'space', 'BackSpace');
        showInstruc(wpnt, 'Calib', instFolder, 'space', 'BackSpace');
        if withTobii
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
       
            % init
            EThndl          = Titta(settings);
            % EThndl          = EThndl.setDummyMode();    % just for internal testing, enabling dummy mode for this readme makes little sense as a demo
            EThndl.init();
            
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
                    tobii.calVal{1}         = EThndl.calibrate(wpnt);
                end
            end
            ListenChar(0);
            showInstruc(wpnt, 'T1P', instFolder, 'space', 'BackSpace');
            % later:
            EThndl.buffer.start('gaze');
            WaitSecs(.8);   % wait for eye tracker to start and gaze to be picked up
            monWidth = EThndl.geom.displayArea.width/10; % in cm
            monHeight = EThndl.geom.displayArea.height/10;% in cm
        end
        
        passed = 0;
        CtrGrad = 10.^ linspace(log10(0.6), log10(0.12), 10); % contrast gradient in pretrials, from 0.6 to 0.12
        if length(bgContrast)~=1
            pre_bgContrast = mean(bgContrast);
        else
            pre_bgContrast = bgContrast;
        end
        while ~passed
            for pretrial = 1:10
                if withTobii
                    [startT, headDist] = show_fix(wpnt, fixTime, fixClrs, winRect, withTobii, tobiiFreq, EThndl, monWidth, monHeight);
                    EThndl.sendMessage('FIX ON Pre',startT);
                else
                    [startT, headDist] = show_fix(wpnt, fixTime, fixClrs, winRect, withTobii);
                end
                ut = UT(monWidth, scWidth, headDist);
                tgCenter = ut.pix2deg(ut.Pol2Rect([rand(1)*MaxEcc,rand(1)*360]));
                % prep stimuli 
                stiTex=[1,0];   % [withTarget, withoutTarget]
                for i = 1:2
                    % there is no target when GaborWidth set to 0
                    stimulus = genStim(winRect, ut, pre_bgContrast, ...
                                       CtrGrad(pretrial), tgCenter, GaborSF, ...
                                       GaborWidth*stiTex(i), GaborOrient, ...
                                       bgWidth, randi(1000));
                    stiTex(i) = Screen('MakeTexture', wpnt, cat(3,stimulus,stimulus,stimulus)*255); % first stimulus is with target
                end
                targetAt = randi(2);

                % Show stimulus 
                % First interval
                Screen('Drawtexture',wpnt,stiTex(targetAt));
                imgT1 = Screen('Flip',wpnt,startT+fixTime-1/hz/2);                   % bit of slack to make sure requested presentation time can be achieved
                if withTobii
                    EThndl.sendMessage(sprintf('STIM ON: First Target-%.0f  tgContrast-%.2f  bgContrast-%.2f', ...
                           targetAt, CtrGrad(pretrial), pre_bgContrast),imgT1);
                end
                imgT1e = Screen('Flip', wpnt, imgT1+stimT-1/hz/2);
                % Second interval
                Screen('Drawtexture',wpnt,stiTex(3-targetAt));
                imgT2 = Screen('Flip',wpnt,imgT1e+interT-1/hz/2);                   % bit of slack to make sure requested presentation time can be achieved
                if withTobii
                    EThndl.sendMessage(sprintf('STIM ON: Second Target-%.0f  tgContrast-%.2f  bgContrast-%.2f', ...
                           targetAt, CtrGrad(pretrial), pre_bgContrast),imgT1);
                end
                imgT2e = Screen('Flip',wpnt,imgT2+stimT-1/hz/2);
        
                % Check key response
                while 1 
                    checkend;
                    [~, secs, keyCode] = KbCheck;
                    if keyCode(KbName(key1)) + keyCode(KbName(key2))==1 % pressing 2 key is not allowed
                        expKey = key12(targetAt);
                        judgement = keyCode(KbName(expKey));
                        break;
                    end
                    WaitSecs(0.05);
                end

                % Color feedback
                if judgement
                    pointColor = [0, 0, 255]; % blue indicates right fixation
                else
                    pointColor = [255, 255, 0]; % yellow indicates wrong fixation
                end
                Screen('DrawDots', wpnt, bgCenter, 63, pointColor);
                fbT = Screen('Flip',wpnt);
                WaitSecs(0.3);
                if ~judgement  % start the next 10 trials, untill all correct
                    break
                end
                if pretrial ==10
                    passed = 1;
                end
            end 
        end

        showInstruc(wpnt, 'T1F', instFolder, 'space', 'BackSpace');

        for trial = 1:trialNum
            % First draw a fixation point
            if withTobii
                [startT, headDist] = show_fix(wpnt, fixTime, fixClrs, winRect, withTobii, tobiiFreq, EThndl, monWidth, monHeight);
                EThndl.sendMessage('FIX ON',startT);
            else
                [startT, headDist] = show_fix(wpnt, fixTime, fixClrs, winRect, withTobii);
            end
            ut = UT(monWidth, scWidth, headDist);
            tgCenter = ut.deg2pix([results.Xtarg(trial), results.Ytarg(trial)]);

            % Get recommended Contrast level by Pelli (1987)
            tTest=QuestQuantile(q);
            results.tgContrast(trial) = 10.^tTest; % the 3rd column denotes the contrast intensity
    
            stiTex=[1,0];   % [withTarget, withoutTarget]
            seeds = [results.seed1(trial),results.seed2(trial)];
            for i = 1:2
                % there is no target when GaborWidth set to 0
                stimulus = genStim(winRect, ut, results.bgContrast(trial), ...
                                   results.tgContrast(trial), tgCenter, GaborSF, ...
                                   GaborWidth*stiTex(i), GaborOrient, ...
                                   bgWidth,seeds(i));
                stiTex(i) = Screen('MakeTexture', wpnt, cat(3,stimulus,stimulus,stimulus)*255); % first stimulus is with target
            end
            % show on screen and log when it was shown in eye-tracker time.
            % NB: by setting a deadline for the flip, we ensure that the previous
            % screen (fixation point) stays visible for the indic m m m m m m m m ated amount of
            % time. See PsychToolbox demos for further elaboration on this way of
            % timing your script.
            
            % Show stimulus 
            % First interval
            Screen('Drawtexture',wpnt,stiTex(results.target(trial)));
            imgT1 = Screen('Flip',wpnt,startT+fixTime-1/hz/2);                   % bit of slack to make sure requested presentation time can be achieved
            if withTobii
                EThndl.sendMessage(sprintf('STIM ON: First Target-%.0f  tgContrast-%.3f  bgContrast-%.2f', ...
                       results.target(trial), results.tgContrast(trial), results.bgContrast(trial)),imgT1);
            end
            imgT1e = Screen('Flip', wpnt, imgT1+stimT-1/hz/2);
            % Second interval
            Screen('Drawtexture',wpnt,stiTex(3-results.target(trial)));
            imgT2 = Screen('Flip',wpnt,imgT1e+interT-1/hz/2);                   % bit of slack to make sure requested presentation time can be achieved
            if withTobii
                EThndl.sendMessage(sprintf('STIM ON: Second Target-%.0f  tgContrast-%.3f  bgContrast-%.2f', ...
                       results.target(trial), results.tgContrast(trial), results.bgContrast(trial)),imgT1);
            end
            imgT2e = Screen('Flip',wpnt,imgT2+stimT-1/hz/2);
    
            % Check key response
            while 1 
                checkend;
                [~, secs, keyCode] = KbCheck;
                if keyCode(KbName(key1)) + keyCode(KbName(key2))==1 % pressing 2 key is not allowed
                    expKey = key12(results.target(trial));
                    results.judge(trial)= keyCode(KbName(expKey));
                    break;
                end
                WaitSecs(0.05);
            end

            % Color feedback
            if results.judge(trial)
                pointColor = [0, 0, 255]; % blue indicates right fixation
            else
                pointColor = [255, 255, 0]; % yellow indicates wrong fixation
            end
            Screen('DrawDots', wpnt, bgCenter, 63, pointColor);
            fbT = Screen('Flip',wpnt);
            WaitSecs(0.3); 
    
            % record data
            results.RT(trial) = secs-imgT2e;
            
            % record x seconds of data, then clear screen. Indicate stimulus
            % removed, clean up
            % finiT = Screen('Flip',wpnt);
            endT = Screen('Flip',wpnt); % finiT+1-1/hz/2 
            Screen('Close',stiTex(1));
            Screen('Close',stiTex(2));
            q=QuestUpdate(q,tTest,results.judge(trial));  % Add the new datum (actual test intensity and observer response) to the database.
    
            if withTobii
                EThndl.sendMessage(sprintf('STIM OFF: trial-%.0f  RT-%.3f  judge-%.0f', ...
                                                  trial,results.RT(trial),results.judge(trial)),endT);
            end
            % slightly less precise ISI is fine..., about 1s give or take a frame
            % or continue immediately if the above upload actions took longer than
            % a second.  ??????? How ?????????
            WaitSecs(1-(GetSecs-endT));
        end
        
        if withTobii 
            % keep the gaze data in buffer for debuging
            gazeData = EThndl.buffer.peekTimeRange('gaze');
            % stop recording
            % save data to mat file, adding info about the experiment
            dat = EThndl.collectSessionData();
            dat.expt.resolution = winRect(3:4);
            dat.expt.stim       = stimulus;
            if ~exist('./Data', 'dir')
                mkdir('./Data');
            end
            EThndl.saveData(dat, fullfile(cd,sprintf('./Data/Threshold/EYEDat_Sub%.0f_Ses%.0f_%s_%s_%s',subjID, session, location, subjName, datestr(datetime,'yyyymmdd'))), true);
            % shut down
            EThndl.deInit();
        end
        t = 10.^QuestMean(q);
        if nargin == 5
            if saveRaw
                save(sprintf('./Data/Threshold/THREXP_Sub%.0f_Ses%.0f_%s_%s_%s',subjID, session, location, subjName, DTstr))
            end
            save(sprintf('./Data/Threshold/THR_Sub%.0f_Ses%.0f_%s_%s_%s',subjID, session, location, subjName, datestr(datetime,'yyyymmdd')),"results")
        end
        % if you want to (also) save the data to Apache Parquet and json files
        % that can easily be read in Python (Apache Parquet files are supported
        % by Pandas), use:
        % EThndl.saveDataToParquet(dat, fullfile(cd,'t'), true);
        % All gaze data columns and messages can be dumped to tsv files using:
        % EThndl.saveGazeDataToTSV(dat, fullfile(cd,'t'), true);
        
    catch me
        if nargin == 5
            save(sprintf('./Data/Threshold/THRINT_Sub%.0f_Ses%.0f_%s_%s_%s',subjID, session, location, subjName, DTstr))
        else
            save(sprintf('./Data/Threshold/THRINT_%s', DTstr))
        end            
        sca
        ListenChar(0);
        rethrow(me)
    end
    sca

% 



