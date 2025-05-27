% Get the threshold by Quest
%
function [threshold,dat] = getThreshold(wpnt, winRect, scFreq, EThndl, taskPar, repeatTimes, maxTrialDur, subjID, session, location, subjName)
% clear all                                                                                           
% sca
% trialNum=60;
instFolder = './Instructions';
if ~exist('./Data/Threshold','dir')
    mkdir('./Data/Threshold')
end
if ~exist('./Data/Threshold/Interrupted','dir')
    mkdir('./Data/Threshold/Interrupted')
end
DTstr = datestr(datetime, 'yyyymmddTHHMM');
%% Parameters
saveRaw                 = true;    % ~200MB for 72 trials, 10min

% task parameters
fixTime   = 1; % should be more than 0.1
fixClrs   = taskPar.fixClrs;
key1      = taskPar.key1;
key2      = taskPar.key2; % 
restTime  = taskPar.restTime;   % second
blockSize = taskPar.blockSize; % trials

bgWidth     = taskPar.bgWidth;       % background width, in degree
GaborSF     = taskPar.GaborSF;        % cycle per degree
GaborCyc    = taskPar.GaborCyc;       % target width (full width at half maxima, FWHM), n cycle
GaborWidth  = GaborCyc/GaborSF; % target width in degree
GaborOrient = taskPar.GaborOrient;   % Orientation of Garbor

Eccent  = taskPar.Eccent; % Eccentricity of location
Orient  = taskPar.Orient; % angle of location to center, counter-clockwise from positive x-axis;
bgContrast = taskPar.bgContrast;             % maximum contrast of background texture
tobiiFreq = EThndl.frequency;
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
    'tgContrast',0,'bgContrast',bgContrast,'repeatID',1:repeatTimes, ...
    'key1RT',nan,'key2RT',nan,'XeyeFix',nan,'Xtarg',nan, ...
    'YeyeFix',nan,'Ytarg',nan,'FixNum',nan,'Err',nan,'judge',nan, ...
    'headDist',nan);

trialNum = height(results);
results.seed = randi(10*trialNum,[trialNum,1]);
% shuffle
trialOrder = Shuffle(1:trialNum);
results = results(trialOrder,:);

% Prepare parameters for Quest
pThreshold = 0.5;
tGuess = log10(0.05);  % prior threshold estimate of central position
%Threshold "t" is measured on an abstract "intensity" scale, which usually corresponds to log10 contrast.
tGuessSd=1;  % the standard deviation you assign to that guess. Be generous!
beta=3.5;   % slope
delta=0.01; % 
gamma=0.25; % 
grain = 0.01;  % stepsize
range = 3;  % the intensity difference between the largest and smallest intensity,  log10([0.01 1]) =-2   0
% create the initial 'q'
q = QuestCreate(tGuess,tGuessSd,pThreshold,beta,delta,gamma,grain,range);
q.normalizePdf = 1;
q = QuestUpdate(q, log10(0.2),1);
q = QuestUpdate(q, log10(0.01),0);

% ATTENTION!! 
% Take the centre of the screen as the coordinate origin, right and up as 
% the positive direction
results.Xtarg = results.ECC .* cosd(results.Orient);
results.Ytarg = results.ECC .* sind(results.Orient);
taskSpace= results(:,{'Xtarg','Ytarg'});
[~, ia, ~] = unique(taskSpace,"rows");
taskSpace = taskSpace(ia, :);
judger = spatialJudge(taskSpace);

%% Initiation of stimulation


try
    ListenChar(0);
    showInstruc(wpnt, 'TP', instFolder, 'space', 'BackSpace');
    
    % later:
    EThndl.buffer.start('gaze');
    WaitSecs(.8);   % wait for eye tracker to start and gaze to be picked up
    monWidth = EThndl.geom.displayArea.width/10;    % in cm
    monHeight = EThndl.geom.displayArea.height/10;  % in cm
    scWidth = winRect(3);  % width of screen resolution, in pixel
    scHeight = winRect(4); % height of screen resolution, in pixel
    bgCenter = round([scWidth/2, scHeight/2]);

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
            imgT = Screen('Flip',wpnt,startT+fixTime-1/scFreq/2);                   % bit of slack to make sure requested presentation time can be achieved
            SysImgT = EThndl.getTimeAsSystemTime(imgT);
            EThndl.sendMessage(sprintf('STIM ON: P. ECC-%.0f  Ori-%.0f  tgContrast-%.3f  bgContrast-%.2f', ...
                   tgCenter_(1), tgCenter_(2), CtrGrad(pretrial), pre_bgContrast), imgT);

            % Check key response
            timeCost = 0;
            while 1 
                checkend;
                [~, secs, keyCode] = KbCheck;
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
            if err<1
                pointColor = [0, 0, 255]; % blue indicates right fixation
            else
                pointColor = [255, 255, 0]; % yellow indicates wrong fixation
            end
            Screen('Drawtexture',wpnt, stiTex);
            Screen('DrawDots', wpnt, lastFixPix_, 14, pointColor,[],3);
            Screen('Flip',wpnt);
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
        oper = showInstruc(wpnt, 'Check', instFolder, 'space', 'BackSpace');
        if oper==-1
            passed = 0;
            tobii.calVal{1} = EThndl.calibrate(wpnt,3,tobii.calVal{1}); 
            EThndl.buffer.start('gaze');
            WaitSecs(.8);   % wait for eye tracker to start and gaze to be picked up
        end
    end
    
    showInstruc(wpnt, 'T1F', instFolder, 'space', 'BackSpace');
    
    for trial = 1:trialNum
        % First draw a fixation point
        [startT, headDist] = show_fix(wpnt, fixTime, fixClrs, winRect, true, tobiiFreq, EThndl, monWidth, monHeight);

        EThndl.sendMessage('FIX ON',startT);

        ut = UT(monWidth, scWidth, headDist);
        tgCenter = ut.deg2pix([results.Xtarg(trial), results.Ytarg(trial)]);
        
        % prep stimuli 
        
        % Get recommended Contrast level by Pelli (1987)
        tTest=QuestQuantile(q);
        results.tgContrast(trial) = 10.^tTest; % the 3rd column denotes the contrast intensity

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
        imgT = Screen('Flip',wpnt,startT+fixTime-1/scFreq/2);                   % bit of slack to make sure requested presentation time can be achieved
        SysImgT = EThndl.getTimeAsSystemTime(imgT);
        EThndl.sendMessage(sprintf('STIM ON: F_%0.f_%0.f ECC-%.0f  Ori-%.0f  tgContrast-%.3f  bgContrast-%.2f', ...
               results.ECC(trial), results.Orient(trial), results.ECC(trial), results.Orient(trial), results.tgContrast(trial), results.bgContrast(trial)),imgT);
        timeCost = 0;
        while 1 
            checkend;
            [~, secs, keyCode] = KbCheck;
            if keyCode(KbName(key1))|| timeCost>=maxTrialDur
                while 1
                    checkend;
                    [keyIsDown2, secs2, keyCode2] = KbCheck;
                    if keyCode2(KbName(key2))|| timeCost>=maxTrialDur
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
        [judgement, err]= judger.judge(lastFixDeg, [results.Xtarg(trial),results.Ytarg(trial)]);
        if judgement == 1
            pointColor = [0, 0, 255]; % blue indicates right fixation
        else
            pointColor = [255, 255, 0]; % yellow indicates wrong fixation
        end
        Screen('Drawtexture',wpnt, stiTex);
        Screen('DrawDots', wpnt, lastFixPix_, 14, pointColor,[],3);
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
        
        % 2.3s as expected searching time
        q=QuestUpdate(q,tTest,results.key2RT(trial)<2.4);  % Add the new datum (actual test intensity and observer response) to the database.
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

    EThndl.saveData(dat, fullfile(cd,sprintf('./Data/Threshold/Dat_Sub%.0f_Ses%.0f',subjID, session)), true);

    save(sprintf('./Data/Threshold/Result_Sub%.0f_Ses%.0f_%s_%s_%s',subjID, session, location, subjName, DTstr),"results")
    if saveRaw
        save(sprintf('./Data/Threshold/EXP_Sub%.0f_Ses%.0f_%s_%s_%s',subjID, session, location, subjName, DTstr))
    end
    % if you want to (also) save the data to Apache Parquet and json files
    % that can easily be read in Python (Apache Parquet files are supported
    % by Pandas), use:
    % EThndl.saveDataToParquet(dat, fullfile(cd,'t'), true);
    % All gaze data columns and messages can be dumped to tsv files using:
    % EThndl.saveGazeDataToTSV(dat, fullfile(cd,'t'), true);
    
    % shut down
    EThndl.buffer.stop('gaze');

    % update the SubjInfo.csv
    infoFilePath = './Data/SubjInfo.csv';
    SubjInfo = readtable(infoFilePath);
    rowIdx = find(SubjInfo.subjID == subjID,1);
    threshold = 10.^QuestQuantile(q);
    SubjInfo.threshold(rowIdx) = threshold;
    writetable(SubjInfo,infoFilePath);


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

    save(sprintf('./Data/Threshold/Interrupted/EXPINT_Sub%.0f_Ses%.0f_%s_%s_%s',subjID, session, location, subjName, DTstr))
    sca
    ListenChar(0);
    rethrow(me)
end

% 
% 



