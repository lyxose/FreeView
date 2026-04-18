function run_freeview_eyelink(wpnt, winRect, el, edfFile, subjInfo, bgCenter, monWidth, ~, bgWidth, bgContrast, nTrials, durationSec, spotProb, dummyMode)
% run_freeview_eyelink  Present a block of free-view trials and save to a separate EDF
% - After calibration, before practice: shows center fixation until gaze
%   is at center, then displays random background for durationSec seconds.
% - Records EyeLink data to a dedicated EDF file (prefixed FV), then reopens
%   the main experiment EDF.
% - Emits messages 'FV ON'/'FV OFF' with trial/seed info for downstream parsing.
%
% INPUT:
%   subjInfo - struct with optional fields: subjID, session, location, subjName, subjGender, subjAge
%              Can be empty struct or have partial fields. Dialog will always show for confirmation.
%   spotProb - optional, probability of spot-detection trials (default 0.2)
%   dummyMode - optional, set true to skip EyeLink dependency and run key-only
%               trials; fixation, drift correction, and gaze-dependent metrics
%               are skipped and saved as NaN.

% Standalone mode: if called without arguments, use defaults from FreeView_main
% and initialize PTB + EyeLink here, prompting for subj/session and suffix.
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

standalone = (nargin == 0);
customFvEdf = false;  % whether we override FV EDF host filename here
% fvSaveSuffix = '';  % suffix disabled for standalone FreeView-only runs
if nargin < 14 || isempty(dummyMode)
    dummyMode = false;
end
dummyMode = logical(dummyMode);

% Auto-fallback: if no live EyeLink connection exists, force dummy mode.
if ~dummyMode
    % Initialize EyeLink connection + calibration
    try
        if ~EyelinkInit(0)
            dummyMode = true;
            fprintf('[FV] EyelinkInit failed, auto-switching to dummy mode.\n');
        end
    catch
        dummyMode = true;
        fprintf('[FV] Error occurred while initializing EyeLink, auto-switching to dummy mode.\n');
    end
end

if standalone
    % Try to put project libs on path based on this file location
    try
        thisDir = fileparts(mfilename('fullpath'));
        projRoot = thisDir; % in project root
        addpath(genpath(fullfile(projRoot, 'function_library')));
        addpath(fullfile(projRoot, 'function_library_cus'));
    catch
    end

    % Defaults aligned with FreeView_main / EyeLink version
    bgClr = 127.5;
    monWidth   = 50.92;    % cm (default lab monitor width)
    bgWidth    = 22;      % deg
    bgContrast = 0.2;     % background contrast
    nTrials    = 150;      % free-view trials
    durationSec= 5;      % seconds per trial
    spotProb   = 0.2;     % probability of spot-detection trial

    % Initialize keyboard
    KbName('UnifyKeyNames');
    
    % Release keyboard listening for dialog input (avoid key accumulation from PTB)
    ListenChar(0);
    
    % Initialize empty subjInfo structure
    subjInfo = struct();
    
    % Get subject information via InformationBox_FV (always shows dialog)
    [subjID, session, location, subjName, ~, ~] = InformationBox_FV(subjInfo);
    
    if isnan(subjID) || isnan(session)
        ListenChar(2);
        error('Invalid Subject ID or Session');
    end

    % Open PTB window
    AssertOpenGL;
    screenId = 1;% max(Screen('Screens'));
    [wpnt, winRect] = PsychImaging('OpenWindow', screenId, bgClr, [], [], [], [], 4);
    Screen('BlendFunction', wpnt, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    HideCursor(wpnt);
    bgCenter = round([winRect(3)/2, winRect(4)/2]);

    if ~dummyMode
        % Suffix disabled: use fixed EDF naming without user suffix.
        % prompts = {'FV EDF suffix (for FVS*S*_suffix)'};
        % defs    = {'A'};
        % answ = inputdlg(prompts, 'FreeView EyeLink Setup', 1, defs);
        % if isempty(answ)
        %     ListenChar(2);
        %     return
        % end
        % fvSaveSuffix = regexprep(answ{1}, '[^A-Za-z0-9]', '');


        if ~dummyMode
            el = EyelinkInitDefaults(wpnt);
            % 关闭所有声音 
            el.targetbeep = 0;
            el.feedbackbeep = 0; 
            el.calibrationtargetbeep = 0; 
            el.calibrationfeedbackbeep = 0;

            EyelinkUpdateDefaults(el);
            Eyelink('Command','screen_pixel_coords = %ld %ld %ld %ld', 0, 0, winRect(3)-1, winRect(4)-1);
            Eyelink('Message','DISPLAY_COORDS %ld %ld %ld %ld', 0, 0, winRect(3)-1, winRect(4)-1);
            Eyelink('Command', 'calibration_type = HV9');
            Eyelink('Command', 'button_function 5 "accept_target_fixation"');
            EyelinkDoTrackerSetup(el);

            % Compose FV EDF name for Host (<=8 chars): FVS{Sub}S{Ses}
            hostBase = sprintf('FVS%dS%d', subjID, session);
            hostBase = regexprep(hostBase, '[^A-Za-z0-9]', '');
            fvEdf = hostBase(1:min(8, numel(hostBase)));
            customFvEdf = true;

            % No main EDF to reopen in standalone
            edfFile = '';
        else
            el = [];
            fvEdf = '';
            % fvSaveSuffix = '';
            customFvEdf = false;
        end
    else
        el = [];
        fvEdf = '';
        % fvSaveSuffix = '';
        customFvEdf = false;
        fprintf('[FV] Dummy mode enabled: running without EyeLink, fixation, or drift correction.\n');
    end
    
    % Re-engage keyboard listening for experiment
    ListenChar(2);  % Restore keyboard listening (unformatted mode)
end

if nargin < 11
    nTrials = 20;
end
if nargin < 12
    durationSec = 5;
end
if nargin < 13
    spotProb = 0.2;
end

spotProb = min(max(spotProb, 0), 1);
maskRadius_deg = 1.5;     % center circular hole radius (deg)
spotRadius_deg = 1/3;    % colorful spot radius (deg)
spotFadeDurSec = 2.5;    % spot grow/fade duration (sec)
fixtime = 1;
fixJitterSec = 0.6;
restEveryN = 20;          % rest every N trials
restTimeSec = 10;         % default rest time (s)
fixDotRadiusPix = max(4, round(winRect(3)/150));
fixHoleRadiusPix = max(1, round(fixDotRadiusPix * 0.35));
catchRespKey = KbName('RightControl');

% Handle subject information for non-standalone mode
% If subjInfo not provided or incomplete, initialize and prompt
if nargin < 5 || isempty(subjInfo)
    subjInfo = struct();
end

% Prompt for subject/session/location information in non-standalone mode.
if ~standalone
    [subjID, session, location, subjName, ~, ~] = InformationBox_FV(subjInfo);
end

% Compose a secondary EDF name (<=8 chars) only when EyeLink is active.
if ~dummyMode
    if ~customFvEdf
        if nargin < 4 || isempty(edfFile)
            edfFile = 'MAIN';
        end
        fvEdf = ['FV' edfFile];
        if length(fvEdf) > 8
            fvEdf = fvEdf(1:8);
        end
    end
else
    fvEdf = '';
end

preambleText = sprintf('FREEVIEW BLOCK Subject %d Session %d', subjID, session);
if ~dummyMode
    % Close current file if any, then open the free-view EDF
    try
        Eyelink('CloseFile');
    catch
    end
    failOpenFV = Eyelink('OpenFile', fvEdf);
    if failOpenFV ~= 0
        error('Cannot create free-view EDF file %s', fvEdf);
    end

    % Annotate
    Eyelink('Command', 'add_file_preamble_text "%s"', preambleText);
end
% Prepare free-view trial log table (similar to main experiment results table)
fvResults = table((1:nTrials)', ...
    nan(nTrials,1), nan(nTrials,1), nan(nTrials,1), ...
    nan(nTrials,1), nan(nTrials,1), nan(nTrials,1), ...
    false(nTrials,1), ...
    nan(nTrials,1), nan(nTrials,1), ...
    nan(nTrials,1), nan(nTrials,1), ...
    nan(nTrials,1), nan(nTrials,1), ...
    nan(nTrials,1), nan(nTrials,1), ...
    nan(nTrials,1), nan(nTrials,1), nan(nTrials,1), ...
    nan(nTrials,1), nan(nTrials,1), ...
    false(nTrials,1), false(nTrials,1), false(nTrials,1), ...
    nan(nTrials,1), nan(nTrials,1), nan(nTrials,1), ...
    'VariableNames', {'trial','seed','headDistCm','driftStatus', ...
    'trialColorR','trialColorG','trialColorB', ...
    'spotTrial','spotXdeg','spotYdeg','spotXpix','spotYpix', ...
    'catchXdeg','catchYdeg','catchXpix','catchYpix', ...
    'spotColorR','spotColorG','spotColorB','spotOnsetSec','spotRespRTsec', ...
    'spotResponded','spotMissed','falseAlarm','driftErrDeg','catchErrDeg','trialDurSec'});
fvResults.thisfixTime = nan(nTrials,1);

% Build a fixed-count catch schedule so each participant gets exactly
% round(nTrials * spotProb) catch trials, randomly ordered across trials.
nSpotTarget = round(nTrials * spotProb);
nSpotTarget = min(max(nSpotTarget, 0), nTrials);
spotSchedule = false(nTrials, 1);
if nSpotTarget > 0
    spotSchedule(randperm(nTrials, nSpotTarget)) = true;
end

% Guidance and practice block before formal trials.
run_fv_guidance_and_practice(wpnt, winRect, bgCenter, monWidth, bgWidth, ...
    bgContrast, maskRadius_deg, spotRadius_deg, spotFadeDurSec, catchRespKey, dummyMode, el);

% Run nTrials free-view trials
for t = 1:nTrials
    if t > 1 && mod(t-1, restEveryN) == 0
        quit = disp_rest(wpnt, bgCenter(1), bgCenter(2), restTimeSec);
        if quit
            fprintf('[FV] Aborted during rest at trial %d/%d\n', t, nTrials);
            break
        end
    end

    if ~dummyMode
        Eyelink('Command', 'record_status_message "FREEVIEW %d/%d"', t, nTrials);
        Eyelink('Message', sprintf('TRIALID %d', t));
        Eyelink('SetOfflineMode');
        Eyelink('Command', 'clear_screen 0');

        % Use center drift correction like the main EyeLink script and redraw
        % fixation point after key acceptance to prevent immediate disappearance.
        Screen('FillRect', wpnt, 127);
        % Draw hollow fixation for drift-correction confirmation.
        Screen('gluDisk', wpnt, 0, bgCenter(1), bgCenter(2), fixDotRadiusPix);
        Screen('gluDisk', wpnt, 127, bgCenter(1), bgCenter(2), fixHoleRadiusPix);
        Screen('Flip', wpnt);
        driftStatus = EyelinkDoDriftCorrection(el, bgCenter(1), bgCenter(2), 0, 1);
        % Drift-correction accepted: switch to solid fixation to indicate pass.
        Screen('FillRect', wpnt, 127);
        Screen('gluDisk', wpnt, 0, bgCenter(1), bgCenter(2), fixDotRadiusPix);
        Screen('Flip', wpnt);

        Eyelink('StartRecording');
        Eyelink('Message', 'FIX ON FV');
        thisfixTime = fixtime + (rand-0.5) * fixJitterSec - 0.3;  % press the space bar take ~0.3s
        WaitSecs(thisfixTime);
    else
        driftStatus = NaN;
        thisfixTime = 0;
    end

    headDist = 62;
    ut = UT(monWidth, winRect(3), headDist);
    if ~dummyMode
        gazeAtFixPix = get_latest_gaze_px(el, bgCenter);
        gazeAtFixFromCenter = gazeAtFixPix - bgCenter;
        gazeAtFixFromCenter(2) = -gazeAtFixFromCenter(2);
        driftErrDeg = norm(ut.pix2deg(gazeAtFixFromCenter));
    else
        driftErrDeg = NaN;
    end

    % Generate background-only stimulus (no target)
    seed = randi(1e6);
    tgCenter = ut.deg2pix([0, 0]);
    stimulus = genStim(winRect, ut, bgContrast, 0, tgCenter.*[1,-1]+bgCenter, 1, 0.01, 0, bgWidth, seed);

    % Dig a circular hole in the center, same strategy as PTB EyeLink main script
    maskRadius_pix = ut.deg2pix(maskRadius_deg);
    stimHeight = size(stimulus, 1);
    stimWidth = size(stimulus, 2);
    [meshX, meshY] = meshgrid(1:stimWidth, 1:stimHeight);
    stimCenter_pix = [stimWidth / 2, stimHeight / 2];
    distFromCenter = sqrt((meshX - stimCenter_pix(1)).^2 + (meshY - stimCenter_pix(2)).^2);
    maskMat = distFromCenter < maskRadius_pix;
    stimulus(maskMat) = 0.5;

    % Apply a single-hue color family per trial. The texture varies mainly in
    % lightness/saturation while staying within one hue band.
    [stimRGB, trialColor] = make_single_hue_texture(stimulus);

    stiTex = Screen('MakeTexture', wpnt, uint8(stimRGB * 255));

    % Spot-detection configuration for this trial (fixed count, random order)
    spotTrial = spotSchedule(t);
    spotDeg = [NaN, NaN];
    spotPix = [NaN, NaN];
    spotColor = [NaN, NaN, NaN];
    spotOnsetSec = NaN;
    if spotTrial
        spotEccMin = maskRadius_deg + 0.5;
        spotEccMax = bgWidth/2 - 0.5;
        spotEcc = spotEccMin + rand * max(spotEccMax - spotEccMin, 0.01);
        spotOri = rand * 360;
        spotDeg = [spotEcc * cosd(spotOri), spotEcc * sind(spotOri)];
        spotPix = ut.deg2pix(spotDeg) .* [1, -1] + bgCenter;
        spotColor = random_vivid_rgb();
        spotOnsetSec = 0.5 + 2.0 * rand;
    end

    fvResults.seed(t) = seed;
    fvResults.headDistCm(t) = headDist;
    fvResults.driftStatus(t) = driftStatus;
    fvResults.trialColorR(t) = trialColor(1);
    fvResults.trialColorG(t) = trialColor(2);
    fvResults.trialColorB(t) = trialColor(3);
    fvResults.spotTrial(t) = spotTrial;
    fvResults.spotXdeg(t) = spotDeg(1);
    fvResults.spotYdeg(t) = spotDeg(2);
    fvResults.spotXpix(t) = spotPix(1);
    fvResults.spotYpix(t) = spotPix(2);
    fvResults.catchXdeg(t) = spotDeg(1);
    fvResults.catchYdeg(t) = spotDeg(2);
    fvResults.catchXpix(t) = spotPix(1);
    fvResults.catchYpix(t) = spotPix(2);
    fvResults.spotColorR(t) = spotColor(1);
    fvResults.spotColorG(t) = spotColor(2);
    fvResults.spotColorB(t) = spotColor(3);
    fvResults.spotOnsetSec(t) = spotOnsetSec;
    fvResults.driftErrDeg(t) = driftErrDeg;
    fvResults.thisfixTime(t) = thisfixTime;

    Screen('DrawTexture', wpnt, stiTex);
    imgT = Screen('Flip', wpnt);
    if ~dummyMode
        Eyelink('Message', sprintf('FV ON: trial-%.0f seed-%.0f bgContrast-%.3f col-[%.2f %.2f %.2f] spot-%.0f', ...
            t, seed, bgContrast, trialColor(1), trialColor(2), trialColor(3), spotTrial));
    end

    % Hold display for durationSec; spot trials can end early with RightControl key
    t0 = GetSecs();
    spotAnnounced = false;
    spotOnsetAbsSec = NaN;
    trialEndedByResponse = false;
    falseAlarmThisTrial = false;
    % KbReleaseWait;
    while GetSecs() - t0 < durationSec
        % Allow emergency exit via experiment's key handling
        try
            checkend;
        catch
        end

        elapsed = GetSecs() - t0;
        Screen('DrawTexture', wpnt, stiTex);

        spotVisible = false;
        if spotTrial && elapsed >= spotOnsetSec
            spotVisible = true;
            if ~spotAnnounced && ~dummyMode
                Eyelink('Message', sprintf('FV SPOT ON: trial-%.0f xDeg-%.3f yDeg-%.3f col-[%.2f %.2f %.2f]', ...
                    t, spotDeg(1), spotDeg(2), spotColor(1), spotColor(2), spotColor(3)));
                spotAnnounced = true;
            end
            fadeP = min(1, (elapsed - spotOnsetSec) / spotFadeDurSec);
            fadeP = max(0, fadeP);
            % Exponential ease-in for both alpha and size: slow start, faster later.
            kEase = 3.0;
            fadeEase = (exp(kEase * fadeP) - 1) / (exp(kEase) - 1);
            spotMaxDiamPix = 2 * max(4, round(ut.deg2pix(spotRadius_deg)));
            spotDrawDiamPix = max(1, round(spotMaxDiamPix * fadeEase));
            alphaVal = round(255 * fadeEase);
            dotColorRGBA = [spotColor * 255, alphaVal];
            Screen('DrawDots', wpnt, spotPix, spotDrawDiamPix, dotColorRGBA, [], 1);
        end

        vbl = Screen('Flip', wpnt);
        if spotVisible && isnan(spotOnsetAbsSec)
            spotOnsetAbsSec = vbl;
        end

        [keyIsDown, keySec, keyCode] = KbCheck;
        respPressed = keyIsDown && any(keyCode(catchRespKey));
        if respPressed
            if spotVisible
                if ~dummyMode
                    gazePix = get_latest_gaze_px(el, bgCenter);
                    gazePixFromCenter = gazePix - bgCenter;
                    gazePixFromCenter(2) = -gazePixFromCenter(2);
                    gazeDeg = ut.pix2deg(gazePixFromCenter);
                    spotErrDeg = norm(gazeDeg - spotDeg);
                else
                    spotErrDeg = NaN;
                end

                if isnan(spotOnsetAbsSec)
                    spotOnsetAbsSec = t0 + spotOnsetSec;
                end
                spotRT = keySec - spotOnsetAbsSec;
                fprintf('[FV] trial %d 注视点目标距离 = %.3f deg, RT(target onset) = %.3f s\n', ...
                    t, spotErrDeg, spotRT);
                if ~dummyMode
                    Eyelink('Message', sprintf('FV SPOT RESP: trial-%.0f errDeg-%.3f RT-%.3f', ...
                        t, spotErrDeg, spotRT));
                end
                fvResults.spotResponded(t) = true;
                fvResults.spotRespRTsec(t) = spotRT;
                fvResults.catchErrDeg(t) = spotErrDeg;
                trialEndedByResponse = true;
                break
            elseif spotTrial
                fprintf('[FV] trial %d RightControl ignored (spot not started yet, elapsed=%.3fs, onset=%.3fs)\n', ...
                    t, elapsed, spotOnsetSec);
                if ~dummyMode
                    Eyelink('Message', sprintf('FV RC EARLY: trial-%.0f elapsed-%.3f onset-%.3f', t, elapsed, spotOnsetSec));
                end
            else
                fprintf('[FV] trial %d RightControl pressed in non-spot trial (monitor only)\n', t);
                if ~dummyMode
                    Eyelink('Message', sprintf('FV RC NONSPOT: trial-%.0f elapsed-%.3f', t, elapsed));
                end
                falseAlarmThisTrial = true;
            end
        end

        WaitSecs(0.01);
    end

    if spotTrial && ~trialEndedByResponse
        fprintf('[FV] trial %d spot missed (no valid RightControl response)\n', t);
        if ~dummyMode
            Eyelink('Message', sprintf('FV SPOT MISS: trial-%.0f', t));
        end
        fvResults.spotMissed(t) = true;
    end

    % End trial
    trialDurSec = GetSecs() - imgT;
    if ~dummyMode
        Eyelink('Message', sprintf('FV OFF: trial-%.0f dur-%.3f endedByResp-%.0f driftErrDeg-%.3f', t, trialDurSec, trialEndedByResponse, driftErrDeg));
    end
    fvResults.trialDurSec(t) = trialDurSec;
    fvResults.falseAlarm(t) = falseAlarmThisTrial;
    Screen('Flip', wpnt);
    Screen('Close', stiTex);
    if ~dummyMode
        Eyelink('StopRecording');
    end
    WaitSecs(0.2 + rand(1)*0.2);
end

% Save trial table and full FV workspace like main script style
if ~exist('./Data', 'dir'), mkdir('./Data'); end
if ~exist('./Data/FV', 'dir'), mkdir('./Data/FV'); end

DTstr = char(datetime('now', 'Format', 'yyyyMMdd''T''HHmm'));
resultFile = sprintf('./Data/FV/Result_FV_Sub%.0f_Ses%.0f_%s_%s_%s.mat', subjID, session, location, subjName, DTstr);
save(resultFile, 'fvResults');

expFile = sprintf('./Data/FV/EXP_FV_Sub%.0f_Ses%.0f_%s_%s_%s.mat', subjID, session, location, subjName, DTstr);
save(expFile);

% End-of-block quality summary for experimenter monitoring
nSpotTrials = sum(fvResults.spotTrial);
nNoSpotTrials = sum(~fvResults.spotTrial);
nMiss = sum(fvResults.spotMissed);
nFalseAlarm = sum(fvResults.falseAlarm);

if nSpotTrials > 0
    missRate = nMiss / nSpotTrials;
else
    missRate = NaN;
end

if nNoSpotTrials > 0
    falseAlarmRate = nFalseAlarm / nNoSpotTrials;
else
    falseAlarmRate = NaN;
end

fprintf('[FV Summary] False alarm: %d/%d (%.2f%%) | Miss: %d/%d (%.2f%%)\n', ...
    nFalseAlarm, nNoSpotTrials, falseAlarmRate*100, ...
    nMiss, nSpotTrials, missRate*100);

% Close and transfer the free-view EDF
if ~dummyMode
    Eyelink('SetOfflineMode');
    Eyelink('Command', 'clear_screen 0');
    WaitSecs(0.5);
    Eyelink('CloseFile');

    try
        fprintf('Receiving free-view data file ''%s''\n', fvEdf);
        status = Eyelink('ReceiveFile');
        if status > 0
            fprintf('ReceiveFile status %d\n', status);
        end
        if exist([fvEdf '.edf'], 'file')
            % Ensure output folder exists
            if ~exist('./Data', 'dir'), mkdir('./Data'); end
            if ~exist('./Data/FV', 'dir'), mkdir('./Data/FV'); end
            fprintf('Free-view file ''%s'' saved to ''%s''\n', fvEdf, pwd);
            DTstr = char(datetime('now', 'Format', 'yyyyMMdd''T''HHmm'));
            if customFvEdf
                % Suffix disabled: save as FVS{Sub}S{Ses}_{DT}.edf
                saveStem = sprintf('FVS%dS%d_%s', subjID, session, DTstr);
            else
                % Backward-compatible naming
                saveStem = sprintf('%s_Sub%d_Ses%d_%s', fvEdf, subjID, session, DTstr);
            end
            newName = fullfile('./Data/FV', [saveStem '.edf']);
            movefile([fvEdf '.edf'], newName);
            fprintf('Free-view EDF moved to %s\n', newName);
        end
    catch
        fprintf('Problem receiving free-view data file ''%s''\n', fvEdf);
    end
end

if standalone
    % Standalone: clean up and close
    if ~dummyMode
        Eyelink('Shutdown');
    end
    try
        ShowCursor(wpnt);
    catch
    end
    sca;
else
    % Reopen the main experiment EDF to continue
    if ~dummyMode
        failOpen = Eyelink('OpenFile', edfFile);
        if failOpen ~= 0
            error('Cannot re-open main EDF file %s after free-view block', edfFile);
        end
        Eyelink('Command', 'add_file_preamble_text "RESUME EXPERIMENT Subject %d Session %d"', subjID, session);
    end
end

end

function [stimRGB, trialColor] = make_single_hue_texture(stimulus)
% Colorize a grayscale texture using one hue per trial.
% The texture varies mainly in value (brightness) and slightly in saturation,
% so the entire pattern stays in one color family.
h = rand();
s = 0.55;
vMin = 0;
vMax = 1;
stimRGB = 0.5 * ones([size(stimulus), 3]);
stimMask = abs(stimulus - 0.5) > 1e-6;
v = vMin + (vMax - vMin) .* (stimulus(stimMask));
hsvImg = cat(3, h * ones(nnz(stimMask), 1), s * ones(nnz(stimMask), 1), min(max(v, 0), 1));
stimRGBMask = hsv2rgb(hsvImg);
stimRGB(repmat(stimMask, [1, 1, 3])) = reshape(stimRGBMask, [], 1);
trialColor = hsv2rgb([h, s, vMax]);
end

function rgb = random_vivid_rgb()
% Random vivid color in RGB space using HSV sampling.
rgb = hsv2rgb([rand, 0.9 + 0.1 * rand, 0.9 + 0.1 * rand]);
end

function gazePix = get_latest_gaze_px(el, fallbackPix)
% Return latest valid gaze in pixel coordinates.
gazePix = fallbackPix;
if Eyelink('NewFloatSampleAvailable') > 0
    evt = Eyelink('NewestFloatSample');
    eyeUsed = Eyelink('EyeAvailable');
    if eyeUsed == el.BINOCULAR
        eyeUsed = el.LEFT_EYE;
    end
    if eyeUsed == el.LEFT_EYE && evt.gx(1) ~= -32768 && evt.gy(1) ~= -32768
        gazePix = [evt.gx(1), evt.gy(1)];
    elseif eyeUsed == el.RIGHT_EYE && evt.gx(2) ~= -32768 && evt.gy(2) ~= -32768
        gazePix = [evt.gx(2), evt.gy(2)];
    else
        if evt.gx(1) ~= -32768 && evt.gy(1) ~= -32768
            gazePix = [evt.gx(1), evt.gy(1)];
        elseif evt.gx(2) ~= -32768 && evt.gy(2) ~= -32768
            gazePix = [evt.gx(2), evt.gy(2)];
        end
    end
end
end

function run_fv_guidance_and_practice(wpnt, winRect, bgCenter, monWidth, bgWidth, bgContrast, maskRadius_deg, spotRadius_deg, spotFadeDurSec, catchRespKey, dummyMode, el)
% Show short task instruction and a mandatory practice loop via images.
spaceKey = KbName('space');
backspaceKey = KbName('BackSpace');
instFolder = './Instructions';



repeatPractice = true;
while repeatPractice
    showInstruc(wpnt, 'FV_Guide_01_TaskAndFix', instFolder, 'space', 'BackSpace');
    % Step 1: fixation preview and confirmation
    run_fixation_demo(wpnt, bgCenter, winRect, dummyMode, el);
    % wait_for_any_key(spaceKey);

    % Step 2-3: practice flow instruction (no-target then with-target)
    showInstruc(wpnt, 'FV_Guide_02_PracticeFlow1', instFolder, 'space', 'BackSpace');
    % wait_for_any_key(spaceKey);

    % Practice A: no-target demo
    run_practice_trial(wpnt, winRect, bgCenter, monWidth, bgWidth, bgContrast, ...
        maskRadius_deg, spotRadius_deg, spotFadeDurSec, catchRespKey, false);
    showInstruc(wpnt, 'FV_Guide_02_PracticeFlow2', instFolder, 'space', 'BackSpace');

    % Practice B: with-target demo
    run_practice_trial(wpnt, winRect, bgCenter, monWidth, bgWidth, bgContrast, ...
        maskRadius_deg, spotRadius_deg, spotFadeDurSec, catchRespKey, true);

    showInstruc(wpnt, 'FV_Guide_03_PracEndChoice', instFolder, 'space', 'BackSpace');
    pressedKey = wait_for_any_key([spaceKey, backspaceKey]);
    repeatPractice = any(pressedKey == backspaceKey);
end

Screen('FillRect', wpnt, 127);
Screen('Flip', wpnt);
KbReleaseWait;
end

function run_practice_trial(wpnt, winRect, bgCenter, monWidth, bgWidth, bgContrast, maskRadius_deg, spotRadius_deg, spotFadeDurSec, catchRespKey, withSpot)
% Present one practice trial with or without dynamic spot.
ut = UT(monWidth, winRect(3), 62);
[stiTex, spotPix, spotColor, spotOnsetSec] = create_practice_texture( ...
    wpnt, winRect, ut, bgCenter, bgContrast, bgWidth, maskRadius_deg, withSpot);

t0 = GetSecs();
maxDurSec = 5;
spotOn = false;
while GetSecs() - t0 < maxDurSec
    elapsed = GetSecs() - t0;
    Screen('DrawTexture', wpnt, stiTex);

    if withSpot && elapsed >= spotOnsetSec
        spotOn = true;
        fadeP = min(1, (elapsed - spotOnsetSec) / spotFadeDurSec);
        fadeP = max(0, fadeP);
        kEase = 3.0;
        fadeEase = (exp(kEase * fadeP) - 1) / (exp(kEase) - 1);
        spotMaxDiamPix = 2 * max(4, round(ut.deg2pix(spotRadius_deg)));
        spotDrawDiamPix = max(1, round(spotMaxDiamPix * fadeEase));
        alphaVal = round(255 * fadeEase);
        dotColorRGBA = [spotColor * 255, alphaVal];
        Screen('DrawDots', wpnt, spotPix, spotDrawDiamPix, dotColorRGBA, [], 1);
    end

    Screen('Flip', wpnt);

    [keyIsDown, ~, keyCode] = KbCheck;
    if keyIsDown
        if withSpot && spotOn && any(keyCode(catchRespKey))
            break
        end
    end
    WaitSecs(0.01);
end

Screen('Close', stiTex);
Screen('FillRect', wpnt, 127);
Screen('Flip', wpnt);
WaitSecs(0.2);
KbReleaseWait;
end

function [stiTex, spotPix, spotColor, spotOnsetSec] = create_practice_texture(wpnt, winRect, ut, bgCenter, bgContrast, bgWidth, maskRadius_deg, withSpot)
% Build a practice texture and optional demo spot at one of 8 fixed directions.
seed = randi(1e6);
tgCenter = ut.deg2pix([0, 0]);
stimulus = genStim(winRect, ut, bgContrast, 0, tgCenter.*[1,-1]+bgCenter, 1, 0.01, 0, bgWidth, seed);

maskRadius_pix = ut.deg2pix(maskRadius_deg);
stimHeight = size(stimulus, 1);
stimWidth = size(stimulus, 2);
[meshX, meshY] = meshgrid(1:stimWidth, 1:stimHeight);
stimCenter_pix = [stimWidth / 2, stimHeight / 2];
distFromCenter = sqrt((meshX - stimCenter_pix(1)).^2 + (meshY - stimCenter_pix(2)).^2);
stimulus(distFromCenter < maskRadius_pix) = 0.5;

[stimRGB, ~] = make_single_hue_texture(stimulus);
stiTex = Screen('MakeTexture', wpnt, uint8(stimRGB * 255));

spotPix = [NaN, NaN];
spotColor = [NaN, NaN, NaN];
spotOnsetSec = NaN;
if withSpot
    dirSet = 22.5 + (0:7) * 45;
    spotOri = dirSet(randi(numel(dirSet)));
    spotEccMin = maskRadius_deg + 0.8;
    spotEccMax = bgWidth/2 - 0.8;
    if spotEccMax <= spotEccMin
        spotEcc = max(maskRadius_deg + 0.2, bgWidth/4);
    else
        spotEcc = min(max(3.0, spotEccMin), spotEccMax);
    end
    spotDeg = [spotEcc * cosd(spotOri), spotEcc * sind(spotOri)];
    spotPix = ut.deg2pix(spotDeg) .* [1, -1] + bgCenter;
    spotColor = random_vivid_rgb();
    spotOnsetSec = 0.8;
end
end

function run_fixation_demo(wpnt, bgCenter, winRect, dummyMode, el)
% Mimic formal fixation sequence: hollow dot -> drift correction -> solid dot.
fixDotRadiusPix = max(4, round(winRect(3)/150));
fixHoleRadiusPix = max(1, round(fixDotRadiusPix * 0.35));
spaceKey = KbName('space');

Screen('FillRect', wpnt, 127);
% Hollow fixation before drift correction
Screen('gluDisk', wpnt, 0, bgCenter(1), bgCenter(2), fixDotRadiusPix);
Screen('gluDisk', wpnt, 127, bgCenter(1), bgCenter(2), fixHoleRadiusPix);
Screen('Flip', wpnt);

if dummyMode
    % In dummy mode, skip eye data checks and only wait for Space.
    wait_for_any_key(spaceKey);
else
    Eyelink('SetOfflineMode');
    Eyelink('Command', 'clear_screen 0');
    EyelinkDoDriftCorrection(el, bgCenter(1), bgCenter(2), 0, 1);
end

% Drift correction accepted: switch to solid fixation.
Screen('FillRect', wpnt, 127);
Screen('gluDisk', wpnt, 0, bgCenter(1), bgCenter(2), fixDotRadiusPix);
Screen('Flip', wpnt);
WaitSecs(0.2);
KbReleaseWait;
end


function pressedKey = wait_for_any_key(keyIndices)
% Wait for any key in keyIndices and return matched index/indices.
pressedKey = [];
while isempty(pressedKey)
    [keyIsDown, ~, keyCode] = KbCheck;
    if keyIsDown
        keyHit = find(keyCode);
        keyHit = keyHit(ismember(keyHit, keyIndices));
        if ~isempty(keyHit)
            pressedKey = keyHit;
        end
    end
    WaitSecs(0.01);
end
KbReleaseWait;
end
