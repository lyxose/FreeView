%% getFix_eyelink  Batch-convert EyeLink EDF into Dat_Sub*_Ses*.mat + expT (I2MC)
% 适配asc文件与edf同名的情况，避免误删和误跳过

%% Setup
clear variables; clear global; close all; fclose('all'); clc
dbstop if error

rootDir = fileparts(mfilename('fullpath'));          % FreeView_Paper_sup
projRoot = fileparts(rootDir);                       % FreeView
rawDir = fullfile(projRoot, 'Data', 'Formal');
processedDir = fullfile(rootDir, 'Processed_data');
if ~isfolder(processedDir); mkdir(processedDir); end

addpath(genpath(fullfile(projRoot, 'function_library')));
addpath(genpath(fullfile(projRoot, 'function_library_cus')));

edf2ascExe = 'edf2asc';   % adjust if not on PATH
rewrite = false;          % set true to force re-run even if expT exists

%% Locate EDF files
edfFiles = dir(fullfile(rawDir, 'S*S*_Sub*_Ses*_*.edf'));
pat = '^S\d+S\d+_Sub(\d+)_Ses(\d+)_\d{8}T\d{4}\.edf$';
keep = ~cellfun(@isempty, regexp({edfFiles.name}.', pat, 'tokens'));
edfFiles = edfFiles(keep);
if isempty(edfFiles)
    warning('No EDF files found under %s', rawDir);
end

for i =1:numel(edfFiles)%:-1:1
    toks = regexp(edfFiles(i).name, pat, 'tokens', 'once');
    subID = str2double(toks{1});
    sesID = str2double(toks{2});
    stem = sprintf('Sub%d_Ses%d', subID, sesID);
    fprintf('Processing %s (EDF: %s)\n', stem, edfFiles(i).name);

    edfPath = fullfile(rawDir, edfFiles(i).name);
    ascPath = strrep(edfPath, '.edf', '.asc');
    datMatPath = fullfile(processedDir, sprintf('Dat_%s.mat', stem));

    % Matching EXP/Result mats (copied to Processed_data for consistency)
    expSrc = dir(fullfile(rawDir, sprintf('EXP_Sub%d_Ses%d*.mat', subID, sesID)));
    resSrc = dir(fullfile(rawDir, sprintf('Result_Sub%d_Ses%d*.mat', subID, sesID)));
    if isempty(resSrc)
        warning('Result mat missing for %s — skipping.', stem);
        continue;
    end
    expMatPath = '';
    if ~isempty(expSrc)
        expMatPath = fullfile(processedDir, expSrc(1).name);
        if ~exist(expMatPath, 'file')
            copyfile(fullfile(rawDir, expSrc(1).name), expMatPath);
        end
    end
    resMatPath = fullfile(processedDir, resSrc(1).name);
    if ~exist(resMatPath, 'file')
        copyfile(fullfile(rawDir, resSrc(1).name), resMatPath);
    end

    % Skip if expT already exists and rewrite is false
    if exist(datMatPath, 'file') && ~rewrite
        vars = who('-file', datMatPath);
        if ismember('expT', vars)
            fprintf('  expT already exists for %s; skipping.\n', stem);
            continue;
        end
    end

    %% Step 1: EDF -> ASC if needed (同名asc)
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
    eyeUsed = detect_eye_from_asc(ascPath);  % 'LEFT' or 'RIGHT'

    %% Step 3: Build Dat_ (Tobii-compatible gaze struct)
    [expt, geometry, settings, systemInfo] = derive_session_meta(expMatPath, rawSmp);
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

    % Save Dat_ file (overwrite if rewrite, otherwise create once)
    save(datMatPath, 'data', 'expt', 'geometry', 'systemInfo', 'settings');

    %% Step 4: Run I2MC per trial and append expT
    expT = build_expT_with_I2MC(datMatPath, resMatPath, rawSmp, msgs, eyeUsed);
    try
        save(datMatPath, 'expT', '-append');
    catch
        warning('Failed to append expT to %s', datMatPath);
    end
end

fprintf('EyeLink batch conversion finished. Outputs in %s\n', processedDir);

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

function [expt, geometry, settings, systemInfo] = derive_session_meta(expMatPath, rawSmp)
    % Derive screen/geometry/settings with sensible fallbacks
    expt = struct(); geometry = struct(); settings = struct(); systemInfo = struct();
    expt.winRect = [0 0 1920 1080];
    geometry.displayArea.width = 511;   % mm
    geometry.displayArea.height = 287;  % mm
    systemInfo.model = 'EyeLink';
    settings.freq = estimate_freq(rawSmp(:,1));

    if exist(expMatPath, 'file')
        tmp = load(expMatPath);
        if isfield(tmp, 'exp') && isstruct(tmp.exp)
            if isfield(tmp.exp, 'winRect') && numel(tmp.exp.winRect) == 4
                expt.winRect = double(tmp.exp.winRect(:))';
            end
            if isfield(tmp.exp, 'displayArea')
                if isfield(tmp.exp.displayArea, 'width'); geometry.displayArea.width = double(tmp.exp.displayArea.width); end
                if isfield(tmp.exp.displayArea, 'height'); geometry.displayArea.height = double(tmp.exp.displayArea.height); end
            end
        end
    end
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

function expT = build_expT_with_I2MC(datMatPath, resMatPath, rawSmp, msgs, eyeUsed)
    S = load(datMatPath);
    sess.data = S.data; sess.expt = S.expt; sess.geometry = S.geometry; sess.settings = S.settings;
    if isfield(S, 'systemInfo'); sess.systemInfo = S.systemInfo; else; sess.systemInfo.model = 'EyeLink'; end

    % Parse trial markers
    [starts, ends_, eccs, oris, trials, fixStarts] = parse_trials(msgs);
    numTrials = numel(starts);

    R = load(resMatPath);
    if ~isfield(R, 'results'); error('Result mat missing results table: %s', resMatPath); end
    results = R.results;
    if ~ismember('ECC', results.Properties.VariableNames); results.ECC = zeros(height(results),1); end
    if ~ismember('Orient', results.Properties.VariableNames); results.Orient = zeros(height(results),1); end
    if ~ismember('Xtarg', results.Properties.VariableNames); results.Xtarg = zeros(height(results),1); end
    if ~ismember('Ytarg', results.Properties.VariableNames); results.Ytarg = zeros(height(results),1); end

    for k = 1:numTrials
        trIdx = trials(k);
        if trIdx>=1 && trIdx<=height(results)
            results.ECC(trIdx) = eccs(k);
            results.Orient(trIdx) = oris(k);
            results.Xtarg(trIdx) = results.ECC(trIdx) * cosd(results.Orient(trIdx));
            results.Ytarg(trIdx) = results.ECC(trIdx) * sind(results.Orient(trIdx));
        end
    end

    expT = results;
    datTemplate = struct('time', [], 'left', struct(), 'right', struct(), 'average', struct(), ...
                         'finalweights', [], 'fix', struct(), 'I2MCopt', struct());
    expT.dat = repmat(datTemplate, height(results), 1);

    ts = rawSmp(:,1); lx = rawSmp(:,2); ly = rawSmp(:,3);
    
    % Handle monocular data only (4 columns: time, xp, yp, pupilSize)
    if size(rawSmp, 2) ~= 4
        error('Expected monocular data with 4 columns (time, x, y, pupil), got %d columns.', size(rawSmp, 2));
    end
    
    pupil = rawSmp(:,4);
    
    % Create dat structure based on which eye was detected
    % I2MC expects either dat.left or dat.right with X, Y, pupil fields
    if strcmp(eyeUsed, 'LEFT')
        % Monocular LEFT eye: lx, ly are left eye data
        datFieldName = 'left';
    else
        % Monocular RIGHT eye: lx, ly are right eye data, but I2MC expects them in 'left'
        % (I2MC always uses left field internally)
        datFieldName = 'left';
    end

    opt = struct();
    opt.xres = sess.expt.winRect(3);
    opt.yres = sess.expt.winRect(4);
    opt.missingx = nan;
    opt.missingy = nan;
    opt.scrSz = [sess.geometry.displayArea.width sess.geometry.displayArea.height] / 10; % mm -> cm
    opt.freq = sess.settings.freq;
    if opt.freq > 120
        opt.downsamples = [2 5 10]; opt.chebyOrder = 8;
    elseif opt.freq == 120
        opt.downsamples = [2 3 5]; opt.chebyOrder = 7;
    else
        opt.downsampFilter = false; opt.downsamples = [2 3];
    end
    opt.maxMergeDist = 15; opt.minFixDur = 60;

    for k = 1:numTrials
        trIdx = trials(k);
        t_fix = fixStarts(k);   % Trial data starts from FIX ON
        t_on = starts(k);       % STIM ON is time alignment point (time = 0)
        t_off = ends_(k);       % Trial data ends at STIM OFF
        
        % Extract data from FIX ON to STIM OFF
        if isnan(t_fix)
            % If FIX ON not found, fallback to STIM ON
            t_fix = t_on;
        end
        
        inTrial = ts >= t_fix & ts <= t_off;
        if ~any(inTrial) || trIdx < 1 || trIdx > height(expT)
            continue;
        end

        data = [ts(inTrial) lx(inTrial) ly(inTrial)];
        dat = struct();
        % Align to STIM ON (t_on as zero point)
        % Times before STIM ON will be negative
        dat.time = data(:,1) - double(t_on);
        
        % All monocular data goes to 'left' field for I2MC (which always uses left internally)
        dat.left = struct('X', data(:,2), 'Y', data(:,3), 'pupil', pupil(inTrial));

        if ismember('headDist', expT.Properties.VariableNames) && ~isnan(expT.headDist(trIdx))
            opt.disttoscreen = expT.headDist(trIdx);
        else
            opt.disttoscreen = 62;
        end

        try
            [fix, datOut] = I2MCfunc(dat, opt);
        catch ME
            warning('I2MC failed on trial %d: %s', trIdx, ME.message);
            fix = struct('startT', [], 'endT', [], 'dur', [], 'xpos', [], 'ypos', []);
            datOut = dat;
        end

        datOut.fix = fix;
        datOut.I2MCopt = opt;
        reqFields = fieldnames(datTemplate);
        for fi = 1:numel(reqFields)
            f = reqFields{fi};
            if ~isfield(datOut, f)
                datOut.(f) = datTemplate.(f);
            end
        end
        expT.dat(trIdx) = datOut;
        expT.FixNum(trIdx) = size(fix.startT, 1);
    end
end

function [starts, ends_, eccs, oris, trials, fixStarts] = parse_trials(msgs)
    % Parse trial markers, excluding practice trials (FIX ON Pre)
    % Returns:
    %   starts : STIM ON time (STIM ON为时间对齐的0点)
    %   ends_  : STIM OFF time
    %   fixStarts : FIX ON time (trial数据提取的起点)
    %   eccs, oris, trials : 对应的参数
    
    nMsg = size(msgs,1);
    starts = nan(nMsg,1); ends_ = nan(nMsg,1); eccs = nan(nMsg,1); oris = nan(nMsg,1); trials = nan(nMsg,1); fixStarts = nan(nMsg,1);
    stimOnCount = 0;
    currentFixStart = NaN;  % Track FIX ON for current trial
    
    for i = 1:nMsg
        tstr = msgs{i,2};
        
        % FIX ON (practice): skip without starting trial tracking
        if contains(tstr, 'FIX ON Pre')
            % Practice trial - ignore
            continue;
        end
        
        % FIX ON (formal): start tracking (not FIX ON Pre)
        if contains(tstr, 'FIX ON') && ~contains(tstr, 'FIX ON Pre')
            currentFixStart = double(msgs{i,1});
        end
        
        % STIM ON: begin formal trial record
        if contains(tstr, 'STIM ON')
            ts_on = double(msgs{i,1});
            tr = parse_number(tstr, 'trial-');
            ec = parse_number(tstr, 'ECC-');
            or = parse_number(tstr, 'Ori-');
            stimOnCount = stimOnCount + 1;
            starts(stimOnCount) = ts_on;
            ends_(stimOnCount) = NaN;
            eccs(stimOnCount) = ec;
            oris(stimOnCount) = or;
            trials(stimOnCount) = tr;
            fixStarts(stimOnCount) = currentFixStart;  % Store associated FIX ON
        elseif contains(tstr, 'STIM OFF')
            ts_off = double(msgs{i,1});
            idx_open = find(isnan(ends_(1:stimOnCount)), 1, 'first');
            if isempty(idx_open); idx_open = find(isnan(ends_(1:stimOnCount)), 1, 'last'); end
            if ~isempty(idx_open); ends_(idx_open) = ts_off; end
        end
    end
    
    % Trim to valid trials
    starts = starts(1:stimOnCount); 
    ends_ = ends_(1:stimOnCount); 
    eccs = eccs(1:stimOnCount); 
    oris = oris(1:stimOnCount); 
    trials = trials(1:stimOnCount);
    fixStarts = fixStarts(1:stimOnCount);
    
    % Keep only trials with both ON and OFF messages
    valid = ~isnan(ends_);
    starts = starts(valid); 
    ends_ = ends_(valid); 
    eccs = eccs(valid); 
    oris = oris(valid); 
    trials = trials(valid);
    fixStarts = fixStarts(valid);
end

function val = parse_number(str, key)
    val = NaN;
    idx = strfind(str, key);
    if isempty(idx); return; end
    s = str(idx(1)+length(key):end);
    tok = regexp(s, '([-+]?\d*\.?\d+)', 'match');
    if ~isempty(tok); val = str2double(tok{1}); end
end
