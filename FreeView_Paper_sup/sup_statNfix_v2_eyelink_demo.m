% Smoke-test pipeline for EyeLink data (v2, 60 trials)
% Assumes getFix_eyelink.m has already processed raw EDF files under
% Data/Formal into Processed_data/Dat_Sub*_Ses*.mat with expT appended.
% This script now only runs the minimal ANA preprocessing (build_fixTable)
% and writes a small CSV to verify downstream compatibility.

%% User config -----------------------------------------------------------
rootDir    = fileparts(mfilename('fullpath')); % FreeView_Paper_sup
projRoot   = fileparts(rootDir);               % FreeView
processedDir = fullfile(rootDir, 'Processed_data');
outDir     = fullfile(rootDir, 'Results_test');
maxTrials  = 480;                              % v2 smoke length
vers       = {'v2_test'};                     % label for outputs
if ~exist(outDir,'dir'); mkdir(outDir); end
if ~exist(processedDir,'dir'); mkdir(processedDir); end

%% Paths for ANA helpers -------------------------------------------------
addpath(genpath(fullfile(projRoot, 'FreeView_Paper', 'Methods', 'function_library_cus', 'ANA')));
addpath(genpath(fullfile(projRoot, 'function_library')));
addpath(genpath(fullfile(projRoot, 'function_library_cus')));

%% Step: Run ANA getFix + stat-like analysis -----------------------------
% Flags to reduce work during smoke test
RUN_DATA_PREP = true;
PLOT_HEATMAP = false; PLOT_BASIC_STATS = false; PLOT_ANG_PROP_SCAN = false;
PLOT_TIME_SERIES = false; PLOT_AXIS_EFFECT = false; PLOT_OBLI_EFFECT = false;
PLOT_TRIAL_SLIDING_GAU = false;

% Locate the first processed Dat_ if user did not specify explicitly
datFiles = dir(fullfile(processedDir, 'Dat_Sub*_Ses*.mat'));
if isempty(datFiles)
    error('No Dat_Sub*_Ses*.mat found in %s. Run getFix_eyelink.m first.', processedDir);
end
tok = regexp(datFiles(2).name, 'Dat_Sub(\d+)_Ses(\d+)', 'tokens', 'once');
subID = str2double(tok{1});
sesID = str2double(tok{2});

if RUN_DATA_PREP
    fprintf('Running minimal build_fixTable on processed mats...\n');
    select_sess = 1;
    exclude_sub = [];
    sub_ses_res = [subID, sesID];
    dirs.mat = processedDir;
    dirs.fix = processedDir;
    resfiles = struct('fname', sprintf('Dat_Sub%d_Ses%d', subID, sesID));
    learn_stage_n = 0;
    last_trial = maxTrials;
    skip_corr = false;
    [fixTable, Nsubj, img_width, img_height, ut, center, digPlace] = build_fixTable(select_sess, exclude_sub, sub_ses_res, dirs, resfiles, learn_stage_n, last_trial, skip_corr); %#ok<ASGLU>
    writetable(fixTable, fullfile(outDir, sprintf('ALL_fixTable_%dSubj.csv', Nsubj)));
end

fprintf('Smoke test setup complete. Check %s for generated mats and %s for outputs.\n', processedDir, outDir);
