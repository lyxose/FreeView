%% combine_eyelink_runs
% Combine two interrupted EyeLink runs (same subject/session) into one Dat_/expT
% Example input suffixes: 'Sub2_Ses1_20260112T1222' and 'Sub2_Ses1_20260112T1204'.
% The script will find matching EDF/EXP/Result files (ignoring name inserts like _CAU_WuYan),
% run I2MC per trial on each, concatenate expT, renumber trials, and save
% Dat_SubX_SesY.mat + expT under Processed_data.

clear variables; clear global; close all; fclose('all'); clc; dbstop if error;

%% User input ------------------------------------------------------------
suffix_recent = 'Sub2_Ses1_20260112T1222'; % longer run (later timestamp)
suffix_early  = 'Sub2_Ses1_20260112T1204'; % shorter run (earlier timestamp)
% Explicit trial counts to extract from each run (set [] to auto by messages/results)
n_recent_keep = 420;                        % e.g., 420 trials from recent run
n_early_keep  = 60;                         % e.g., 60 trials from early run
rewrite = true;                             % overwrite combined Dat_ if exists

%% Paths -----------------------------------------------------------------
rootDir = fileparts(mfilename('fullpath'));      % FreeView_Paper_sup
projRoot = fileparts(rootDir);                   % FreeView
rawDir = fullfile(projRoot, 'Data', 'Formal');
processedDir = fullfile(rootDir, 'Processed_data');
if ~isfolder(processedDir); mkdir(processedDir); end

addpath(genpath(fullfile(projRoot, 'function_library')));
addpath(genpath(fullfile(projRoot, 'function_library_cus')));

edf2ascExe = 'edf2asc';

%% Process each run ------------------------------------------------------
runEarly  = process_one_run(rawDir, processedDir, edf2ascExe, suffix_early,  n_early_keep);
runRecent = process_one_run(rawDir, processedDir, edf2ascExe, suffix_recent, n_recent_keep);

% Basic sanity
if runEarly.subID ~= runRecent.subID || runEarly.sesID ~= runRecent.sesID
    error('Subject/session mismatch between runs (%d/%d vs %d/%d).', runEarly.subID, runEarly.sesID, runRecent.subID, runRecent.sesID);
end

% Combine expT (append; keep column consistency)
expT = [runEarly.expT; runRecent.expT];
if ismember('trial', expT.Properties.VariableNames)
    expT.trial = (1:height(expT))';
end

% Save combined Dat_ (reuse later meta; data not needed for downstream ANA)
datMatPath = fullfile(processedDir, sprintf('Dat_Sub%d_Ses%d.mat', runEarly.subID, runEarly.sesID));
if exist(datMatPath,'file') && ~rewrite
    warning('Dat file exists and rewrite=false, skipping save: %s', datMatPath);
else
    data = runRecent.data; %#ok<NASGU> retain last run gaze struct for completeness
    expt = runRecent.expt; %#ok<NASGU>
    geometry = runRecent.geometry; %#ok<NASGU>
    settings = runRecent.settings; %#ok<NASGU>
    systemInfo = runRecent.systemInfo; %#ok<NASGU>
    save(datMatPath, 'data','expt','geometry','settings','systemInfo','expT');
    fprintf('Combined Dat saved: %s\n', datMatPath);
end

fprintf('Done. Total trials combined: %d (early %d + recent %d)\n', height(expT), height(runEarly.expT), height(runRecent.expT));

%% -----------------------------------------------------------------------
function out = process_one_run(rawDir, processedDir, edf2ascExe, suffix, maxKeep)
    % Discover files (robust to inserted names like _CAU_WuYan)
    patEDF = sprintf('*%s*.edf', suffix);
    edfList = dir(fullfile(rawDir, patEDF));
    if isempty(edfList); error('EDF not found for suffix %s', suffix); end
    edfPath = fullfile(rawDir, edfList(1).name);

    % Parse sub/ses
    tok = regexp(suffix, 'Sub(\d+)_Ses(\d+)_\d{8}T\d{4}', 'tokens', 'once');
    subID = str2double(tok{1}); sesID = str2double(tok{2});
    tsTok = regexp(suffix, '(\d{8}T\d{4})', 'tokens', 'once');
    runTS = tsTok{1};

    % Match EXP/Result
    expAll = dir(fullfile(rawDir, sprintf('EXP_Sub%d_Ses%d*.mat', subID, sesID)));
    resAll = dir(fullfile(rawDir, sprintf('Result_Sub%d_Ses%d*.mat', subID, sesID)));
    % Prefer those with matching timestamp substring
    expList = expAll(contains({expAll.name}, runTS));
    resList = resAll(contains({resAll.name}, runTS));
    if isempty(expList); expList = expAll; end
    if isempty(resList); error('Result mat missing for %s', suffix); end
    expMatPath = '';
    if ~isempty(expList)
        expMatPath = fullfile(processedDir, expList(1).name);
        if ~exist(expMatPath,'file'); copyfile(fullfile(rawDir, expList(1).name), expMatPath); end
    end
    resMatPath = fullfile(processedDir, resList(1).name);
    if ~exist(resMatPath,'file'); copyfile(fullfile(rawDir, resList(1).name), resMatPath); end

    % ASC path (same stem as EDF)
    [~, edfBase] = fileparts(edfPath);
    ascPath = fullfile(rawDir, [edfBase '.asc']);
    if ~exist(ascPath,'file')
        fprintf('  Converting EDF -> ASC (%s)\n', edfBase);
        cmd = sprintf('"%s" -y -miss -o "%s" "%s"', edf2ascExe, ascPath, edfPath);
        status = system(cmd);
%         if status~=0 || ~exist(ascPath,'file')
%             error('edf2asc failed for %s (status %d)', edfBase, status);
%         end
    end

    % Load samples/messages
    rawSmp = sup_readSamples_eyelink(ascPath);
    msgs = sup_loadMsgs_eyelink(ascPath);
    if isempty(rawSmp) || isempty(msgs)
        error('Empty samples/messages for %s', ascPath);
    end

    % Build Dat-like meta
    [expt, geometry, settings, systemInfo] = derive_session_meta(expMatPath, rawSmp);
    ts = rawSmp(:,1);
    gxL = rawSmp(:,2); gyL = rawSmp(:,3);
    gxR = rawSmp(:,4); gyR = rawSmp(:,5);
    pupilL = rawSmp(:,6); pupilR = rawSmp(:,7);
    data = struct();
    data.gaze.left.gazePoint.onDisplayArea  = [gxL./expt.winRect(3), gyL./expt.winRect(4)]';
    data.gaze.right.gazePoint.onDisplayArea = [gxR./expt.winRect(3), gyR./expt.winRect(4)]';
    data.gaze.left.gazePoint.valid  = ~isnan(gxL) & ~isnan(gyL);
    data.gaze.right.gazePoint.valid = ~isnan(gxR) & ~isnan(gyR);
    data.gaze.left.pupil.valid  = ~isnan(pupilL);
    data.gaze.right.pupil.valid = ~isnan(pupilR);
    data.gaze.left.pupil.diameter  = pupilL;
    data.gaze.right.pupil.diameter = pupilR;
    data.gaze.systemTimeStamp = ts;

    % Build expT with I2MC
    expT = build_expT_with_I2MC(struct('data',data,'expt',expt,'geometry',geometry,'settings',settings,'systemInfo',systemInfo), resMatPath, rawSmp, msgs, maxKeep, suffix);

    out = struct('subID',subID,'sesID',sesID,'expT',expT,'data',data,'expt',expt,'geometry',geometry,'settings',settings,'systemInfo',systemInfo);
end

%% Helpers --------------------------------------------------------------
function [expt, geometry, settings, systemInfo] = derive_session_meta(expMatPath, rawSmp)
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
    if numel(ts) < 2; f = 1000; return; end
    dt = diff(ts); medDt = median(dt(~isnan(dt)));
    if isempty(medDt) || medDt <= 0; f = 1000; else; f = round(1000 ./ medDt); end
end

function expT = build_expT_with_I2MC(sess, resMatPath, rawSmp, msgs, maxKeep, runLabel)
    [starts, ends_, eccs, oris, trials, numOn] = parse_trials(msgs);
    numTrials = numel(starts);

    R = load(resMatPath);
    if ~isfield(R, 'results'); error('Result mat missing results table: %s', resMatPath); end
    results = R.results;
    if ~ismember('ECC', results.Properties.VariableNames); results.ECC = zeros(height(results),1); end
    if ~ismember('Orient', results.Properties.VariableNames); results.Orient = zeros(height(results),1); end
    if ~ismember('Xtarg', results.Properties.VariableNames); results.Xtarg = zeros(height(results),1); end
    if ~ismember('Ytarg', results.Properties.VariableNames); results.Ytarg = zeros(height(results),1); end

    % Align by available trials and requested maxKeep
    keep = min(numTrials, height(results));
    if ~isempty(maxKeep) && isnumeric(maxKeep) && isfinite(maxKeep)
        if maxKeep > numTrials
            warning('Requested %d trials for %s exceeds available %d by messages (possible missing STIM OFF). Using %d.', maxKeep, runLabel, numTrials, numTrials);
        end
        keep = min(keep, maxKeep);
    end
    starts = starts(1:keep); ends_ = ends_(1:keep); eccs = eccs(1:keep); oris = oris(1:keep); trials = trials(1:keep);
    results = results(1:keep,:);

    for k = 1:keep
        trIdx = k;
        results.ECC(trIdx) = eccs(k);
        results.Orient(trIdx) = oris(k);
        results.Xtarg(trIdx) = results.ECC(trIdx) * cosd(results.Orient(trIdx));
        results.Ytarg(trIdx) = results.ECC(trIdx) * sind(results.Orient(trIdx));
    end

    expT = results;
    datTemplate = struct('time', [], 'left', struct(), 'right', struct(), 'average', struct(), ...
                         'finalweights', [], 'fix', struct(), 'I2MCopt', struct());
    expT.dat = repmat(datTemplate, height(results), 1);

    ts = rawSmp(:,1); lx = rawSmp(:,2); ly = rawSmp(:,3); rx = rawSmp(:,4); ry = rawSmp(:,5); pl = rawSmp(:,6); pr = rawSmp(:,7);

    opt = struct();
    opt.xres = sess.expt.winRect(3);
    opt.yres = sess.expt.winRect(4);
    opt.missingx = nan; opt.missingy = nan;
    opt.scrSz = [sess.geometry.displayArea.width sess.geometry.displayArea.height] / 10;
    opt.freq = sess.settings.freq;
    if opt.freq > 120
        opt.downsamples = [2 5 10]; opt.chebyOrder = 8;
    elseif opt.freq == 120
        opt.downsamples = [2 3 5]; opt.chebyOrder = 7;
    else
        opt.downsampFilter = false; opt.downsamples = [2 3];
    end
    opt.maxMergeDist = 15; opt.minFixDur = 60;

    for k = 1:keep
        t_on = starts(k); t_off = ends_(k);
        inTrial = ts >= t_on & ts <= t_off;
        if ~any(inTrial); continue; end

        data = [ts(inTrial) lx(inTrial) ly(inTrial) rx(inTrial) ry(inTrial) pl(inTrial) pr(inTrial)];
        dat = struct();
        dat.time = data(:,1) - double(t_on);
        dat.left = struct('X', data(:,2), 'Y', data(:,3), 'pupil', data(:,6));
        dat.right = struct('X', data(:,4), 'Y', data(:,5), 'pupil', data(:,7));

        if ismember('headDist', expT.Properties.VariableNames) && ~isnan(expT.headDist(k))
            opt.disttoscreen = expT.headDist(k);
        else
            opt.disttoscreen = 60;
        end

        try
            [fix, datOut] = I2MCfunc(dat, opt);
        catch ME
            warning('I2MC failed on trial %d: %s', k, ME.message);
            fix = struct('startT', [], 'endT', [], 'dur', [], 'xpos', [], 'ypos', []);
            datOut = dat;
        end

        datOut.fix = fix; datOut.I2MCopt = opt;
        reqFields = fieldnames(datTemplate);
        for fi = 1:numel(reqFields)
            f = reqFields{fi};
            if ~isfield(datOut, f); datOut.(f) = datTemplate.(f); end
        end
        expT.dat(k) = datOut;
    end
end

function [starts, ends_, eccs, oris, trials, numOn] = parse_trials(msgs)
    nMsg = size(msgs,1);
    starts = nan(nMsg,1); ends_ = nan(nMsg,1); eccs = nan(nMsg,1); oris = nan(nMsg,1); trials = nan(nMsg,1);
    stimOnCount = 0;
    for i = 1:nMsg
        tstr = msgs{i,2};
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
        elseif contains(tstr, 'STIM OFF')
            ts_off = double(msgs{i,1});
            idx_open = find(isnan(ends_(1:stimOnCount)), 1, 'first');
            if isempty(idx_open); idx_open = find(isnan(ends_(1:stimOnCount)), 1, 'last'); end
            if ~isempty(idx_open); ends_(idx_open) = ts_off; end
        end
    end
    starts = starts(1:stimOnCount); ends_ = ends_(1:stimOnCount); eccs = eccs(1:stimOnCount); oris = oris(1:stimOnCount); trials = trials(1:stimOnCount);
    valid = ~isnan(ends_);
    starts = starts(valid); ends_ = ends_(valid); eccs = eccs(valid); oris = oris(valid); trials = trials(valid);
    numOn = stimOnCount;
    if sum(valid) < stimOnCount
        warning('Detected %d STIM ON but only %d paired OFF; trailing trial(s) incomplete and dropped.', stimOnCount, sum(valid));
    end
end

function val = parse_number(str, key)
    val = NaN;
    idx = strfind(str, key);
    if isempty(idx); return; end
    s = str(idx(1)+length(key):end);
    tok = regexp(s, '([-+]?\d*\.?\d+)', 'match');
    if ~isempty(tok); val = str2double(tok{1}); end
end