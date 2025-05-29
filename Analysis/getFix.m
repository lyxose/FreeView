%% get ready
clear variables; clear global; clear mex; close all; fclose('all'); clc
dbstop if error % for debugging: trigger a debug point when an error occurs

% setup directories
myDir = fileparts(mfilename('fullpath'));
cd(myDir);
cd ..\Data;                        dirs.mat       = cd;
cd ..\Analysis
    if ~isdir('Processed_data') 
        mkdir(fullfile(cd,'Processed_data'));
    end
    cd Processed_data
        if ~isdir('msgs_ophak') %#ok<*ISDIR>
            mkdir(fullfile(cd,'msgs_ophak'));
        end
        cd msgs_ophak;          dirs.msgsO      = cd;
    cd ..;
        if ~isdir('samples_ophak')
            mkdir(fullfile(cd,'samples_ophak'));
        end
        cd samples_ophak;       dirs.samplesO   = cd;
    cd ..; 
        if ~isdir('fixDet') %#ok<*ISDIR>
            mkdir(fullfile(cd,'fixDet'));
        end
        cd fixDet;              dirs.fix        = cd;
    cd ..;
cd ..;  cd function_library;    dirs.funclib    = cd;
cd ..;  cd function_library_cus;dirs.funclib_cus= cd;
cd ..;
addpath(genpath(dirs.funclib));                 % add dirs to path
addpath(genpath(dirs.funclib_cus))

rewrite = false;
%%% I2MC params 
maxMergeDist = 15;
minFixDur    = 60;

%%% get eye tracker files 
% filter so we only get data that matches the filter. uses regexp
[files,~] = FileFromFolder(dirs.mat,[],'mat');
filtstr = '^Dat_Sub(\d+)_Ses(\d+)\.mat$';
matched = regexpi({files.name}.',filtstr,'tokens');
eyefile_idx = ~cellfun(@isempty,matched);
eyefiles   = files(eyefile_idx);
sub_ses = cellfun(@(x) [str2double(x{1}{1}), str2double(x{1}{2})], matched(eyefile_idx), 'UniformOutput', false);
sub_ses = cell2mat(sub_ses(~cellfun(@isempty, matched(eyefile_idx))));
nsess  = length(eyefiles);
%%% get result table
filtstr = '^Result_Sub(\d+)_Ses(\d+)[_a-zA-Z]+\d+T\d+\.mat$';
matched = regexpi({files.name}.',filtstr,'tokens');
resfile_idx = ~cellfun(@isempty,matched);
resfiles   = files(resfile_idx);
sub_ses_res = cellfun(@(x) [str2double(x{1}{1}), str2double(x{1}{2})], matched(resfile_idx), 'UniformOutput', false);
sub_ses_res = cell2mat(sub_ses_res(~cellfun(@isempty, matched(resfile_idx))));
if ~isequal(sub_ses_res,sub_ses)
    unmatched_sub_ses = setdiff(sub_ses, sub_ses_res, 'rows');
    warning('Files not matched... Sub%.0f_Ses%.0f', unmatched_sub_ses(1), unmatched_sub_ses(2))
end

%% 
for p=1:nsess
    fprintf('Subject: %.0f, Session-%.0f\n', sub_ses(p,1), sub_ses(p,2))
    if exist(fullfile(dirs.fix,[files(p).fname '.mat']),"file") && ~rewrite
        disp('Fixation result already exist! Skip to the next data...')
        continue
    end
    % read msgs and data
    sess     = load(fullfile(dirs.mat,eyefiles(p).name));
    try
        scrRes  = sess.expt.winRect(3:4);
    catch
        scrRes  = [1920,1080];
    end
    ts      = sess.data.gaze.systemTimeStamp;
    % the Pro SDK does not guarantee invalid data is nan. Set to nan if
    % invalid
    sess.data.gaze. left.gazePoint.onDisplayArea(:,~sess.data.gaze. left.gazePoint.valid) = nan;
    sess.data.gaze.right.gazePoint.onDisplayArea(:,~sess.data.gaze.right.gazePoint.valid) = nan;
    sess.data.gaze. left.pupil.diameter(~sess.data.gaze. left.pupil.valid) = nan;
    sess.data.gaze.right.pupil.diameter(~sess.data.gaze.right.pupil.valid) = nan;
    % collect data from the file, and turn gaze positions from normalized
    % coordinates into pixels
    samp    = [bsxfun(@times,sess.data.gaze.left.gazePoint.onDisplayArea,scrRes.'); bsxfun(@times,sess.data.gaze.right.gazePoint.onDisplayArea,scrRes.'); sess.data.gaze.left.pupil.diameter; sess.data.gaze.right.pupil.diameter];
    header  = {'t','gaze_point_LX','gaze_point_LY','gaze_point_RX','gaze_point_RY','pupil_diameter_L','pupil_diameter_R'};
    
    % skip threshold stage
    target = '^STIM OFF: F_\d+_\d+ trial-1  FixNum-.*$';  % ^表示字符串开头，\.转义点号
    logical_idx = ~cellfun(@isempty, regexp(sess.messages(:,2), target));
    indices = find(logical_idx);
    if length(indices)==2
        sess.messages = sess.messages(indices(2)-2:end,:);
    end

    [timest,what,msgs] = parseMsgs_Form(sess.messages);
    % load exp result table
    expT = load(fullfile(dirs.mat,resfiles(p).name)).results;
    if ~ismember('headDist', expT.Properties.VariableNames)
        expT.headDist = repmat(68, height(expT), 1);
    end
    numTrials = length(timest.fix);
%     alldat = cell(numTrials,1);
    % for each trial
    for trial=1:numTrials
        qSel = ts>=timest.fix(trial) & ts<=timest.end(trial);
        data = [num2cell(ts(qSel)); num2cell(samp(:,qSel))];
        data = cell2mat(cellfun(@double, data', 'UniformOutput', false));
        
        % event detection
        % make params struct (only have to specify those you want to be
        % different from their defaults)
        opt.xres          = scrRes(1);
        opt.yres          = scrRes(2);
        opt.missingx      = nan;
        opt.missingy      = nan;
        opt.scrSz         = [sess.geometry.displayArea.width sess.geometry.displayArea.height]/10;  % mm -> cm
        opt.disttoscreen  = expT.headDist(trial);
        opt.freq          = sess.settings.freq;
        if opt.freq>120
            opt.downsamples   = [2 5 10];
            opt.chebyOrder    = 8;
        elseif opt.freq==120
            opt.downsamples   = [2 3 5];
            opt.chebyOrder    = 7;
        else
            % 90 Hz, 60 Hz, 30 Hz
            opt.downsampFilter= false;
            opt.downsamples   = [2 3];
        end
        if strcmp(sess.systemInfo.model,'X2-30_Compact')
            if sess.settings.freq==40
                % for some weird reason the X2-30 reports 40Hz even though it is 30
                opt.freq = 30;
            end
        end
        if opt.freq==30
            warning('Be careful about using I2MC with data that is only 30 Hz. In a brief test, this did not appear to work well with the settings in this file.')
        end
        if (~isfield(opt,'downsampFilter') || opt.downsampFilter) && ~exist('cheby1','file')
            warning('By default, I2MC runs a Chebyshev filter over the data as part of its operation. It appears that this filter (the function ''cheby1'' from the signal processing toolbox) is not available in your installation. I am thus disabling the filter.')
            opt.downsampFilter= false;
        end
        opt.maxMergeDist  = maxMergeDist;
        opt.minFixDur     = minFixDur;
        
        % make data struct
        clear dat;
        dat.time        = (data(:,1)-double(timest.start(trial)))./1000; % mu_s to ms, make samples relative to onset of texture
        dat.left.X      = data(:,2);
        dat.left.Y      = data(:,3);
        dat.right.X     = data(:,4);
        dat.right.Y     = data(:,5);
        dat.left.pupil  = data(:,6);    % add pupil data to file. not used by I2MC but good for plotting
        dat.right.pupil = data(:,7);
        [fix,dat]       = I2MCfunc(dat,opt);
        
        % collect info and store
        dat.fix         = fix;
        dat.I2MCopt     = opt;
        
%         alldat{trial} = dat;
        expT.dat(trial) = dat;
    end
    save(fullfile(dirs.fix,[files(p).fname '.mat']),'expT');
end

rmpath(genpath(dirs.funclib));                  % cleanup path

