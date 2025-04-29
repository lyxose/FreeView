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
if ~exist('./Data/Threshold','dir')
    mkdir('./Data/Threshold')
end
if ~exist('./Data/Threshold/Interrupted','dir')
    mkdir('./Data/Threshold/Interrupted')
end
if ~exist('./Data/Formal/Interrupted','dir')
    mkdir('./Data/Formal/Interrupted')
end


DTstr = datestr(datetime, 'yyyymmddTHHMM');
%% Import functions
addpath(genpath('function_library'));
addpath('function_library_cus');
instFolder = './Instructions';
%% Parameters
DEBUGlevel              = 0;
trialNum                = 500;  % total trial number
learnTNum               = 150;  % quest to 0.85 for statistical learning 
learnP                  = 0.85;
tTest                   = 0.4;
connectTNum             = 100;
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
blockSize = 20;  % trials
maxTrialDur = 5; % second
MaxErr = 1;      % max distance of correct judgement in degree
[subjID, session, location, subjName, subjGender, subjAge, threshold] = InformationBox;

tobiiFreq = 250; % hz
tobiiMod = 'Tobii Pro Fusion';

bgWidth = 20;       % background width, in degree
GaborSF = 6;        % cycle per degree
GaborCyc = 2;       % target width (full width at half maxima, FWHM), n cycle
GaborWidth = GaborCyc/GaborSF; % target width in degree
GaborOrient = -45;  % Orientation of Garbor
bgContrast = 0.2;   % maximum contrast of background texture
noiseP= 0.08;       % probability of target out of ROI
rot_ang = -35;       % randomly rotate the rect to avoid influence from rect orientation, in degree
rFix = bgWidth/2+2; % radius of fixations
tgContrast = threshold;
tgSeed = randi(10000);
default_distance = 68;
% R_max = 7;
% R_min = 2;
firstn = 3;% number of fixations to show on experiementer's screen
%% PTB parameters
scr = max(Screen('Screens'));
screens      = Screen('Screens');
screenNumber = max(screens);
scWidth = Screen('Resolution',screenNumber).width; % width of screen resolution, in pixel
scHeight = Screen('Resolution',screenNumber).height; % width of screen resolution, in pixel
bgCenter = round([scWidth/2, scHeight/2]);

%% TASK
try
%% Generate the task space
% A full rectangel
[Eccent, Orient, Ximg, Yimg] = TaskSpace_gapRect([-6.5 9.5], [1.5 1.6], [0 3], trialNum, noiseP, scWidth, scHeight, bgWidth, rot_ang, tgSeed);

%%

% for task = 1:2
%     if threshold==0
%     else
%         tgContrast = threshold; % contrast of target center (Gabor)
%         trialNum = 360;              % formal exp
%     end
%     Orient = rand(trialNum, 1) * 360;  
%     Eccent = sqrt(R_min^2 + (R_max^2 - R_min^2) * rand(trialNum, 1)); % 概率积分变换定理
figure()
scatter(Ximg, Yimg, 30, 'k', 'filled',...
    'MarkerFaceAlpha', 0.3,...
    'MarkerEdgeColor', 'none');
axis([0 scWidth 0 scHeight]);
set(gca, 'YDir', 'reverse', 'Color', [1 1 1]); % 坐标系匹配图像
title(sprintf('%d个空间采样点分布', trialNum));
%% 
% Screen('Preference', 'SkipSyncTests', 1);

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
    results.oriF = rand(trialNum,1)*360;
    
    % ATTENTION!! 
    % Take the centre of the screen as the coordinate origin, right and up as 
    % the positive direction
%     results.Xtarg_img = Ximg;
%     results.Ytarg_img = Yimg;
    results.Xtarg = results.ECC .* cosd(results.Orient);
    results.Ytarg = results.ECC .* sind(results.Orient);
%     taskSpace= results(:,{'Xtarg','Ytarg'});
%     [~, ia, ic] = unique(taskSpace,"rows");
%     taskSpace = taskSpace(ia, :);
%     judger = spatialJudge(taskSpace);
    
    %% Initiation of stimulation
    
    % fix the hardware bugs
    [ keyIsDown, ~, keyCode ] = KbCheck;
    keyCode = find(keyCode, 1);
    if keyIsDown
        ignoreKey=keyCode;
        DisableKeysForKbCheck(ignoreKey);
        fprintf('KEY "%s" was disabled!',KbName(keyCode));
    end
    Screen('Preference', 'SyncTestSettings', 0.002); % the systems are a little noisy, give the test a little more leeway
    KbName('UnifyKeyNames');
    
    FreeViewExp_PTB_TITTA;

    dat = EThndl.collectSessionData();
    dat.expt.restTime    = restTime;
    dat.expt.blockSize   = blockSize;
    dat.expt.winRect     = winRect;
    dat.expt.bgWidth     = bgWidth;     
    dat.expt.GaborSF     = GaborSF;     
    dat.expt.GaborCyc    = GaborCyc;      
    dat.expt.GaborWidth  = GaborWidth;        
    dat.expt.GaborOrient = GaborOrient;  
    dat.expt.tgSeed      = tgSeed;
    dat.expt.rFix        = rFix;
%     if threshold==0
%         EThndl.saveData(dat, fullfile(cd,sprintf('./Data/Threshold/Dat_Sub%.0f_Ses%.0f',subjID, session)), true);
%         save(sprintf('./Data/Threshold/Result_Sub%.0f_Ses%.0f_%s_%s_%s',subjID, session, location, subjName, DTstr),"results")
%         if saveRaw
%             save(sprintf('./Data/Threshold/EXP_Sub%.0f_Ses%.0f_%s_%s_%s',subjID, session, location, subjName, DTstr))
%         end
    EThndl.saveData(dat, fullfile(cd,sprintf('./Data/Formal/Dat_Sub%.0f_Ses%.0f',subjID, session)), true);
    save(sprintf('./Data/Formal/Result_Sub%.0f_Ses%.0f_%s_%s_%s',subjID, session, location, subjName, DTstr),"results")
    if saveRaw
        save(sprintf('./Data/Formal/EXP_Sub%.0f_Ses%.0f_%s_%s_%s',subjID, session, location, subjName, DTstr))
    end

    % update table
    SubjInfo = readtable('./Data/SubjInfo.csv');
    rowIdx = find(SubjInfo.subjID == subjID,1);
    SubjInfo(rowIdx,'threshold') = {results.tgContrast(end)};
    writetable(SubjInfo,'./Data/SubjInfo.csv');
    threshold = results.tgContrast(end);
    % show quest
    figure;
    plot(results.tgContrast);
    xlabel('trial');
    ylabel('contrast');
%     input('Continue? (press ENTER)');
    exportgraphics(gcf,sprintf('./Data/Threshold/Quest_Sub%.0f_Ses%.0f_%s_%s_%s.png',subjID, session, location, subjName, DTstr),'Resolution',300)
%     else
%         EThndl.saveData(dat, fullfile(cd,sprintf('./Data/Formal/Dat_Sub%.0f_Ses%.0f',subjID, session)), true);
%         save(sprintf('./Data/Formal/Resu   lt_Sub%.0f_Ses%.0f_%s_%s_%s',subjID, session, location, subjName, DTstr),"results")
%         if saveRaw
%             save(sprintf('./Data/Formal/EXP_Sub%.0f_Ses%.0f_%s_%s_%s',subjID, session, location, subjName, DTstr))
%         end
%     end

    % if you want to (also) save the data to Apache Parquet and json files
    % that can easily be read in Python (Apache Parquet files are supported
    % by Pandas), use:
    % EThndl.saveDataToParquet(dat, fullfile(cd,'t'), true);
    % All gaze data columns and messages can be dumped to tsv files using:
    % EThndl.saveGazeDataToTSV(dat, fullfile(cd,'t'), true);
    
    % shut down
EThndl.deInit();
% end
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
    dat.expt.tgSeed      = tgSeed;
%     if threshold==0
%         save(sprintf('./Data/Threshold/Interrupted/EXPINT_Sub%.0f_Ses%.0f_%s_%s_%s',subjID, session, location, subjName, DTstr))
%     else
        save(sprintf('./Data/Formal/Interrupted/EXPINT_Sub%.0f_Ses%.0f_%s_%s_%s',subjID, session, location, subjName, DTstr))
%     end
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


