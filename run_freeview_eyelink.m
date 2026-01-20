function run_freeview_eyelink(wpnt, winRect, el, edfFile, subjInfo, bgCenter, monWidth, monHeight, bgWidth, bgContrast, nTrials, durationSec)
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

% Standalone mode: if called without arguments, use defaults from FreeView_main
% and initialize PTB + EyeLink here, prompting for subj/session and suffix.

standalone = (nargin == 0);
customFvEdf = false;  % whether we override FV EDF host filename here
fvSaveSuffix = '';    % user-provided suffix for saved filename
infoCollected = false; % track whether subject info already gathered

if standalone
    % Try to put project libs on path based on this file location
    try
        thisDir = fileparts(mfilename('fullpath'));
        projRoot = fileparts(fileparts(thisDir)); % ../../ → project root
        addpath(genpath(fullfile(projRoot, 'function_library')));
        addpath(fullfile(projRoot, 'function_library_cus'));
    catch
    end

    % Defaults aligned with FreeView_main / EyeLink version
    bgClr = 127;
    monWidth   = 51.1;    % cm (default lab monitor width)
    monHeight  = 28.7;    % cm (unused downstream but kept for consistency)
    bgWidth    = 15;      % deg
    bgContrast = 0.2;     % background contrast
    nTrials    = 20;      % free-view trials
    durationSec= 10;      % seconds per trial

    % Initialize keyboard
    KbName('UnifyKeyNames');
    
    % Release keyboard listening for dialog input (avoid key accumulation from PTB)
    ListenChar(0);
    
    % Initialize empty subjInfo structure
    subjInfo = struct();
    
    % Get subject information via InformationBox_FV (always shows dialog)
    [subjID, session, location, subjName, subjGender, subjAge] = InformationBox_FV(subjInfo);
    infoCollected = true;
    
    % Prompt for FV EDF suffix
    prompts = {'FV EDF suffix (for FVS*S*_suffix)'};
    defs    = {'A'};
    answ = inputdlg(prompts, 'FreeView EyeLink Setup', 1, defs);
    if isempty(answ)
        ListenChar(2);
        return
    end
    fvSaveSuffix = regexprep(answ{1}, '[^A-Za-z0-9]', '');

    if isnan(subjID) || isnan(session)
        ListenChar(2);
        error('Invalid Subject ID or Session');
    end

    % Open PTB window
    AssertOpenGL;
    screenId = max(Screen('Screens'));
    [wpnt, winRect] = PsychImaging('OpenWindow', screenId, bgClr, [], [], [], [], 4);
    Screen('BlendFunction', wpnt, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    HideCursor(wpnt);
    bgCenter = round([winRect(3)/2, winRect(4)/2]);

    % Initialize EyeLink connection + calibration
    dummymode = 0;
    if ~EyelinkInit(dummymode)
        sca; error('EyelinkInit failed.');
    end
    if Eyelink('IsConnected') < 1
        sca; error('EyeLink not connected.');
    end
    el = EyelinkInitDefaults(wpnt);
    EyelinkUpdateDefaults(el);
    Eyelink('Command','screen_pixel_coords = %ld %ld %ld %ld', 0, 0, winRect(3)-1, winRect(4)-1);
    Eyelink('Message','DISPLAY_COORDS %ld %ld %ld %ld', 0, 0, winRect(3)-1, winRect(4)-1);
    Eyelink('Command', 'calibration_type = HV9');
    Eyelink('Command', 'button_function 5 "accept_target_fixation"');
    EyelinkDoTrackerSetup(el);

    % Compose FV EDF name for Host (<=8 chars): FVS{Sub}S{Ses}{suffix}
    hostBase = sprintf('FVS%dS%d%s', subjID, session, fvSaveSuffix);
    hostBase = regexprep(hostBase, '[^A-Za-z0-9]', '');
    fvEdf = hostBase(1:min(8, numel(hostBase)));
    customFvEdf = true;

    % No main EDF to reopen in standalone
    edfFile = '';
    
    % Re-engage keyboard listening for experiment
    ListenChar(2);  % Restore keyboard listening (unformatted mode)
end

if nargin < 11
    nTrials = 20;
end
if nargin < 12
    durationSec = 10;
end

% Handle subject information for non-standalone mode
% If subjInfo not provided or incomplete, initialize and prompt
if nargin < 5 || isempty(subjInfo)
    subjInfo = struct();
end

% Always prompt for confirmation/completion via InformationBox_FV (only once)
if ~infoCollected
    % Release keyboard listening for dialog input
    ListenChar(0);
    [subjID, session, location, subjName, subjGender, subjAge] = InformationBox_FV(subjInfo);
    ListenChar(2);  % Restore keyboard listening
    infoCollected = true;
end

% Compose a secondary EDF name (<=8 chars)
if ~customFvEdf
    fvEdf = ['FV' edfFile];
    if length(fvEdf) > 8
        fvEdf = fvEdf(1:8);
    end
end

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
preambleText = sprintf('FREEVIEW BLOCK Subject %d Session %d', subjID, session);
Eyelink('Command', 'add_file_preamble_text "%s"', preambleText);

% Run nTrials free-view trials
for t = 1:nTrials
    Eyelink('Command', 'record_status_message "FREEVIEW %d/%d"', t, nTrials);
    Eyelink('Message', 'TRIALID %d', t);
    Eyelink('SetOfflineMode');
    Eyelink('Command', 'clear_screen 0');
    Eyelink('StartRecording');
    WaitSecs(0.1);

    % Show fixation until gaze at center
    [~, headDist] = show_fix(wpnt, bgCenter(1), bgCenter(2), 0.8, 0, winRect, 2.0, ...
        'eyeTrackerType', 'EyeLink', 'el', el, 'monWidth', monWidth, 'monHeight', monHeight);
    ut = UT(monWidth, winRect(3), headDist);

    % Generate background-only stimulus (no target)
    seed = randi(1e6);
    tgCenter = ut.deg2pix([0, 0]);
    stimulus = genStim(winRect, ut, bgContrast, 0, tgCenter.*[1,-1]+bgCenter, 1, 0.01, 0, bgWidth, seed);
    
    stiTex = Screen('MakeTexture', wpnt, cat(3, stimulus, stimulus, stimulus)*255);
    Screen('Drawtexture', wpnt, stiTex);
    imgT = Screen('Flip', wpnt);
    Eyelink('Message', sprintf('FV ON: trial-%.0f seed-%.0f bgContrast-%.3f', t, seed, bgContrast));

    % Hold background for durationSec seconds
    t0 = GetSecs();
    while GetSecs() - t0 < durationSec
        % Allow emergency exit via experiment's key handling
        try; checkend(); catch; end
        WaitSecs(0.01);
    end

    % End trial
    Eyelink('Message', sprintf('FV OFF: trial-%.0f dur-%.3f', t, GetSecs() - imgT));
    Screen('Flip', wpnt);
    Screen('Close', stiTex);
    Eyelink('StopRecording');
    WaitSecs(0.2 + rand(1)*0.2);
end

% Close and transfer the free-view EDF
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
        DTstr = datestr(now,'yyyymmddTHHMM');
        if customFvEdf
            % Use requested save pattern: FVS{Sub}S{Ses}_{suffix}_{DT}.edf
            saveStem = sprintf('FVS%dS%d_%s_%s', subjID, session, fvSaveSuffix, DTstr);
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

if standalone
    % Standalone: clean up and close
    Eyelink('Shutdown');
    try
        ShowCursor(wpnt);
    catch
    end
    sca;
else
    % Reopen the main experiment EDF to continue
    failOpen = Eyelink('OpenFile', edfFile);
    if failOpen ~= 0
        error('Cannot re-open main EDF file %s after free-view block', edfFile);
    end
    Eyelink('Command', 'add_file_preamble_text "RESUME EXPERIMENT Subject %d Session %d"', subjID, session);
end

end