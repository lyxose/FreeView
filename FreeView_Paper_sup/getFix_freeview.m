%% getFix_freeview  Batch-convert EyeLink EDF into Dat_Sub*_Ses*.mat (FreeView)
% Processes FreeView EDF files without result tables
% Creates empty trial table and extracts gaze/fixation data via I2MC

%% Setup
clear variables; clear global; close all; fclose('all'); clc
dbstop if error

rootDir = fileparts(mfilename('fullpath'));          % FreeView_Paper_sup
projRoot = fileparts(rootDir);                       % FreeView
rawDir = fullfile(projRoot, 'Data', 'FV');
processedDir = fullfile(projRoot, 'Data', 'FV_Processed');
if ~isfolder(processedDir); mkdir(processedDir); end

addpath(genpath(fullfile(projRoot, 'function_library')));
addpath(genpath(fullfile(projRoot, 'function_library_cus')));

edf2ascExe = 'edf2asc';   % adjust if not on PATH
rewrite = false;          % set true to force re-run even if Dat_ exists

%% Locate EDF files (FreeView files: FVS*S*_Sub*_Ses*.edf or FVS*S*_*.edf)
edfFiles = dir(fullfile(rawDir, 'FVS*.edf'));
if isempty(edfFiles)
    warning('No FV EDF files found under %s', rawDir);
    return;
end

for i = 1:numel(edfFiles)
    fname = edfFiles(i).name;
    fprintf('Processing FreeView EDF: %s\n', fname);
    
    % Parse filename: FVS{?}S{?}_Sub{SubID}_Ses{SesID}_*.edf or FVS{?}S{?}_*.edf
    % Try to extract SubID and SesID from filename
    pat = 'FVS\d+S\d+_Sub(\d+)_Ses(\d+)_';
    tok = regexp(fname, pat, 'tokens', 'once');
    
    if isempty(tok)
        % Alternative pattern without Sub/Ses
        pat_alt = 'FVS(\d+)S(\d+)_';
        tok_alt = regexp(fname, pat_alt, 'tokens', 'once');
        if ~isempty(tok_alt)
            subID = str2double(tok_alt{1});
            sesID = str2double(tok_alt{2});
        else
            warning('Cannot parse SubID/SesID from %s; skipping.', fname);
            continue;
        end
    else
        subID = str2double(tok{1});
        sesID = str2double(tok{2});
    end
    
    stem = sprintf('Sub%d_Ses%d', subID, sesID);
    edfPath = fullfile(rawDir, fname);
    ascPath = strrep(edfPath, '.edf', '.asc');
    datMatPath = fullfile(processedDir, sprintf('Dat_%s.mat', stem));
    
    % Skip if Dat_ already exists and rewrite is false
    if exist(datMatPath, 'file') && ~rewrite
        fprintf('  Dat_%s.mat already exists; skipping.\n', stem);
        continue;
    end
    
    %% Step 1: EDF -> ASC if needed
    if ~exist(ascPath, 'file')
        fprintf('  Converting EDF -> ASC using %s...\n', edf2ascExe);
        cmd = sprintf('"%s" -y -miss -o "%s" "%s"', edf2ascExe, ascPath, edfPath);
        status = system(cmd);
        if status ~= 0 
            warning('edf2asc returned %d.', status);
        end
        if ~exist(ascPath, 'file')
            fprintf('  Failed to generate ASC for %s; skipping.\n', stem);
            continue;
        end
    end
    
    %% Step 2: Read samples/messages
    rawSmp = sup_readSamples_eyelink(ascPath); % [time, xp, yp, pupilSize] - monocular
    if isempty(rawSmp)
        warning('Empty samples for %s; skipping.', stem);
        continue;
    end
    msgs = sup_loadMsgs_eyelink(ascPath);      % {timestamp, text}
    if isempty(msgs)
        warning('Empty messages for %s; skipping.', stem);
        continue;
    end
    
    % Detect which eye (LEFT or RIGHT) from ASC file
    eyeUsed = detect_eye_from_asc(ascPath);
    
    %% Step 3: Build Dat_ (Tobii-compatible gaze struct)
    [expt, geometry, settings, systemInfo] = derive_session_meta();
    ts = rawSmp(:,1);
    xp = rawSmp(:,2); yp = rawSmp(:,3);
    pupilSize = rawSmp(:,4);
    
    data = struct();
    % Map monocular data based on detected eye
    if strcmp(eyeUsed, 'LEFT')
        data.gaze.left.gazePoint.onDisplayArea  = [xp./expt.winRect(3), yp./expt.winRect(4)]';
        data.gaze.left.gazePoint.valid  = ~isnan(xp) & ~isnan(yp);
        data.gaze.left.pupil.valid  = ~isnan(pupilSize);
        data.gaze.left.pupil.diameter  = pupilSize;
        data.gaze.right.gazePoint.onDisplayArea = nan(2, numel(xp));
        data.gaze.right.gazePoint.valid = false(numel(xp), 1);
        data.gaze.right.pupil.valid = false(numel(xp), 1);
        data.gaze.right.pupil.diameter = nan(numel(xp), 1);
    else  % RIGHT eye
        data.gaze.right.gazePoint.onDisplayArea = [xp./expt.winRect(3), yp./expt.winRect(4)]';
        data.gaze.right.gazePoint.valid = ~isnan(xp) & ~isnan(yp);
        data.gaze.right.pupil.valid = ~isnan(pupilSize);
        data.gaze.right.pupil.diameter = pupilSize;
        data.gaze.left.gazePoint.onDisplayArea = nan(2, numel(xp));
        data.gaze.left.gazePoint.valid = false(numel(xp), 1);
        data.gaze.left.pupil.valid = false(numel(xp), 1);
        data.gaze.left.pupil.diameter = nan(numel(xp), 1);
    end
    data.gaze.systemTimeStamp = ts;
    
    % Save Dat_ file
    save(datMatPath, 'data', 'expt', 'geometry', 'systemInfo', 'settings');
    
    %% Step 4: Run I2MC per trial and build expT
    expT = build_expT_with_I2MC_fv(datMatPath, rawSmp, msgs, eyeUsed);
    try
        save(datMatPath, 'expT', '-append');
    catch
        warning('Failed to append expT to %s', datMatPath);
    end
    
    fprintf('  Successfully processed Dat_%s.mat\n', stem);
end

fprintf('FreeView batch conversion finished. Outputs in %s\n', processedDir);

%% Helpers --------------------------------------------------------------
function eyeUsed = detect_eye_from_asc(ascPath)
    % Detect which eye (LEFT or RIGHT) from ASC file header
    % Look for "SAMPLES GAZE    LEFT" or "SAMPLES GAZE    RIGHT"
    % Default to LEFT if not found
    
    eyeUsed = 'LEFT';  % default
    
    try
        fid = fopen(ascPath, 'rt');
        if fid < 1, return; end
        
        cleanup = onCleanup(@() fclose(fid));
        
        while true
            ln = fgetl(fid);
            if ~ischar(ln), break; end
            
            % Look for SAMPLES line with eye designation
            if contains(upper(ln), 'SAMPLES GAZE')
                if contains(upper(ln), ' LEFT')
                    eyeUsed = 'LEFT';
                    break;
                elseif contains(upper(ln), ' RIGHT')
                    eyeUsed = 'RIGHT';
                    break;
                end
            end
        end
    catch
        % If detection fails, use default LEFT
    end
end

function [expt, geometry, settings, systemInfo] = derive_session_meta()
    % Derive screen/geometry/settings with sensible fallbacks
    expt = struct(); geometry = struct(); settings = struct(); systemInfo = struct();
    expt.winRect = [0 0 1920 1080];
    geometry.displayArea.width = 511;   % mm
    geometry.displayArea.height = 287;  % mm
    systemInfo.model = 'EyeLink';
    settings.freq = 1000;               % Default EyeLink sampling rate
end

function f = estimate_freq(ts)
    ts = ts(:);
    if numel(ts) < 2
        f = 1000;
        return;
    end
    dt = diff(ts);
    medDt = median(dt(~isnan(dt)));
    if isempty(medDt) || medDt <= 0
        f = 1000;
    else
        f = round(1000 ./ medDt);
    end
end

function expT = build_expT_with_I2MC_fv(datMatPath, rawSmp, msgs, eyeUsed)
    %% Parse trial markers from messages: 'FV ON' and 'FV OFF'
    [starts, ends_, trials_seed_contrast] = parse_freeview_trials(msgs);
    numTrials = numel(starts);
    
    if numTrials == 0
        warning('No FV trials found in messages.');
        % Create minimal empty table
        expT = table();
        expT.trial = [];
        expT.seed = [];
        expT.bgContrast = [];
        expT.headDist = [];
        expT.dat = [];
        return;
    end
    
    %% Create empty results table with FreeView parameters
    expT = table();
    expT.trial = zeros(numTrials, 1);
    expT.seed = zeros(numTrials, 1);
    expT.bgContrast = nan(numTrials, 1);
    expT.headDist = nan(numTrials, 1);
    
    % Extract trial info from messages
    for k = 1:numTrials
        info = trials_seed_contrast{k};
        expT.trial(k) = info.trial;
        expT.seed(k) = info.seed;
        expT.bgContrast(k) = info.bgContrast;
    end
    
    %% Load pre-processed data
    S = load(datMatPath);
    sess.data = S.data;
    sess.expt = S.expt;
    sess.geometry = S.geometry;
    sess.settings = S.settings;
    if isfield(S, 'systemInfo')
        sess.systemInfo = S.systemInfo;
    else
        sess.systemInfo.model = 'EyeLink';
    end
    
    %% Set up I2MC options
    ts = rawSmp(:,1);
    lx = rawSmp(:,2);
    ly = rawSmp(:,3);
    pupil = rawSmp(:,4);
    
    % Validate monocular data format
    if size(rawSmp, 2) ~= 4
        error('Expected monocular data with 4 columns (time, x, y, pupil), got %d columns.', size(rawSmp, 2));
    end
    
    opt = struct();
    opt.xres = sess.expt.winRect(3);
    opt.yres = sess.expt.winRect(4);
    opt.missingx = nan;
    opt.missingy = nan;
    opt.scrSz = [sess.geometry.displayArea.width sess.geometry.displayArea.height] / 10; % mm -> cm
    opt.freq = sess.settings.freq;
    if opt.freq > 120
        opt.downsamples = [2 5 10];
        opt.chebyOrder = 8;
    elseif opt.freq == 120
        opt.downsamples = [2 3 5];
        opt.chebyOrder = 7;
    else
        opt.downsampFilter = false;
        opt.downsamples = [2 3];
    end
    opt.maxMergeDist = 15;
    opt.minFixDur = 60;
    opt.disttoscreen = 60;  % Default distance in cm
    
    %% Initialize dat field template
    datTemplate = struct('time', [], 'left', struct(), 'right', struct(), 'average', struct(), ...
                         'finalweights', [], 'fix', struct(), 'I2MCopt', struct());
    expT.dat = repmat(datTemplate, numTrials, 1);
    
    %% Process each trial
    for k = 1:numTrials
        t_on = starts(k);
        t_off = ends_(k);
        inTrial = ts >= t_on & ts <= t_off;
        
        if ~any(inTrial)
            warning('Trial %d: no samples in time range [%.0f, %.0f]', k, t_on, t_off);
            continue;
        end
        
        % Extract trial data
        data = [ts(inTrial) lx(inTrial) ly(inTrial)];
        dat = struct();
        dat.time = data(:,1) - double(t_on);
        % All monocular data goes to 'left' field for I2MC
        dat.left = struct('X', data(:,2), 'Y', data(:,3), 'pupil', pupil(inTrial));
        
        % Run I2MC fixation detection
        try
            [fix, datOut] = I2MCfunc(dat, opt);
        catch ME
            warning('I2MC failed on trial %d: %s', k, ME.message);
            fix = struct('startT', [], 'endT', [], 'dur', [], 'xpos', [], 'ypos', []);
            datOut = dat;
        end
        
        % Store result
        datOut.fix = fix;
        datOut.I2MCopt = opt;
        reqFields = fieldnames(datTemplate);
        for fi = 1:numel(reqFields)
            f = reqFields{fi};
            if ~isfield(datOut, f)
                datOut.(f) = datTemplate.(f);
            end
        end
        expT.dat(k) = datOut;
    end
end

function [starts, ends_, trials_info] = parse_freeview_trials(msgs)
    % Parse messages with pattern 'FV ON: trial-X seed-Y bgContrast-Z'
    %                       and 'FV OFF: trial-X dur-Y'
    nMsg = size(msgs, 1);
    starts = nan(nMsg, 1);
    ends_ = nan(nMsg, 1);
    trials_info = cell(nMsg, 1);
    stimOnCount = 0;
    
    for i = 1:nMsg
        ts_msg = double(msgs{i, 1});
        tstr = msgs{i, 2};
        
        if contains(tstr, 'FV ON')
            % Parse: 'FV ON: trial-{trial} seed-{seed} bgContrast-{contrast}'
            tr = parse_number(tstr, 'trial-');
            sd = parse_number(tstr, 'seed-');
            bc = parse_number(tstr, 'bgContrast-');
            
            stimOnCount = stimOnCount + 1;
            starts(stimOnCount) = ts_msg;
            
            info = struct();
            info.trial = tr;
            info.seed = sd;
            info.bgContrast = bc;
            trials_info{stimOnCount} = info;
            
        elseif contains(tstr, 'FV OFF')
            % Parse: 'FV OFF: trial-{trial} dur-{duration}'
            idx_open = find(isnan(ends_(1:stimOnCount)), 1, 'first');
            if isempty(idx_open)
                idx_open = find(isnan(ends_(1:stimOnCount)), 1, 'last');
            end
            if ~isempty(idx_open)
                ends_(idx_open) = ts_msg;
            end
        end
    end
    
    % Trim to valid trials
    starts = starts(1:stimOnCount);
    ends_ = ends_(1:stimOnCount);
    trials_info = trials_info(1:stimOnCount);
    
    % Keep only trials with both ON and OFF messages
    valid = ~isnan(ends_);
    starts = starts(valid);
    ends_ = ends_(valid);
    trials_info = trials_info(valid);
end

function val = parse_number(str, key)
    val = NaN;
    idx = strfind(str, key);
    if isempty(idx)
        return;
    end
    s = str(idx(1)+length(key):end);
    tok = regexp(s, '([-+]?\d*\.?\d+)', 'match');
    if ~isempty(tok)
        val = str2double(tok{1});
    end
end
