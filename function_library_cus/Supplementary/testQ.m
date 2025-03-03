clear; clc;
addpath('../')

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
restTime = 10;   % second
blockSize = 20; % trials
maxTrialDur = 4; % second

monWidth = 54.4785;
monHeight = 30.6442;

bgWidth = 15;       % background width, in degree
GaborSF = 6;        % cycle per degree
GaborCyc = 2;       % target width (full width at half maxima, FWHM), n cycle
GaborWidth = GaborCyc/GaborSF; % target width in degree
GaborOrient = -45;   % Orientation of Garbor
bgContrast = 0.2;             % maximum contrast of background texture
trialNum = 80;              % repeat n times for each condition

R_max = 7;
R_min = 2;
tgSeed = randi(10000);
original_rng = rng;    % original seed
rng(tgSeed);           % use customized seed for each subject
Orient = rand(trialNum, 1) * 360;  
Eccent = sqrt(R_min^2 + (R_max^2 - R_min^2) * rand(trialNum, 1)); % 概率积分变换定理
rng(original_rng);     % 

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

[~,results] = expRun.generateTrialList('ECC',nan,'Orient',nan, ...
    'tgContrast',nan,'bgContrast',bgContrast,'repeatID',1:trialNum, ...
    'key1RT',nan,'key2RT',nan,'XeyeFix',nan,'Xtarg',nan, ...
    'YeyeFix',nan,'Ytarg',nan,'FixNum',nan,'Err',nan,'judge',nan, ...
    'headDist',nan);
results.ECC = Eccent;
results.Orient = Orient;
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
Screen('Preference', 'SkipSyncTests', 1);
scr = max(Screen('Screens'));
screens      = Screen('Screens');
screenNumber = max(screens);
scWidth = Screen('Resolution',screenNumber).width; % width of screen resolution, in pixel
scHeight = Screen('Resolution',screenNumber).height; % width of screen resolution, in pixel
bgCenter = round([scWidth/2, scHeight/2]);
%%    

% Prepare parameters for Quest
pThreshold = 0.5;
tGuess = log10(0.05);  % prior threshold of contrast
%Threshold "t" is measured on an abstract "intensity" scale, which usually corresponds to log10 contrast.
tGuessSd=1;  % the standard deviation you assign to that guess. Be generous!
beta=1.5;delta=0.01;gamma=0; % Be careful!
grain = 0.01;  % stepsize
range = 3;  % the intensity difference between the largest and smallest intensity,  log10([0.01 1]) =-2   0
wrongRight={'wrong','right'};
% create the initial 'q'
q = QuestCreate(tGuess,tGuessSd,pThreshold,beta,delta,gamma,grain,range);
q.normalizePdf = 1;
q = QuestUpdate(q, log10(1),1);
q = QuestUpdate(q, log10(0.01),0);

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

% oper = showInstruc(wpnt, 'welcome', instFolder, 'space', 'BackSpace');

for trial = 1:80
    % First draw a fixation point

    [startT, headDist] = show_fix(wpnt, fixTime, fixClrs, winRect, false);


    ut = UT(monWidth, scWidth, headDist);
    tgCenter = ut.deg2pix([results.Xtarg(trial), results.Ytarg(trial)]);
    % prep stimuli 
    tTest=QuestQuantile(q);
    results.tgContrast(trial) = 10.^tTest;
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
    timeCost = 0;
    while 1 
        checkend;
        [keyIsDown, secs, keyCode] = KbCheck;
        if keyCode(KbName(key1)) || timeCost>=maxTrialDur
            while 1
                checkend;
                [keyIsDown2, secs2, keyCode2] = KbCheck;
                if keyCode2(KbName(key2)) 
                    results.judge(trial)= 1;
                    results.key1RT(trial) = secs-imgT;
                    results.key2RT(trial) = secs2-imgT;
                    break
                elseif timeCost>=maxTrialDur
                    results.judge(trial)= 0;
                    break
                end
                timeCost = WaitSecs(0.01)-imgT; % reduce the sampling rate to reduce the pressure of CPU
            end
            if keyCode2(KbName(key2)) 
                break
            elseif timeCost>=maxTrialDur
                results.judge(trial)= 0;
                break 
            end
        end
        timeCost = WaitSecs(0.01)-imgT;
    end
    

    % quest    
    q=QuestUpdate(q,tTest,results.judge(trial));  % Add the new datum (actual test intensity and observer response) to the database.
    if mod(trial,20)==0
        Screen('Flip',wpnt);
        WaitSecs(10);
    end
end

%%
plot(results.tgContrast(1:80))
xlabel('trial')
ylabel('contrast')