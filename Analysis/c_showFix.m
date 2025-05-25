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

clear variables; clear global; clear mex; close all; fclose('all'); clc

dbstop if error % for debugging: trigger a debug point when an error occurs

% setup directories
myDir = fileparts(mfilename('fullpath'));
cd(myDir);
                                dirs.home       = cd;
cd data;                        dirs.data       = cd;
        cd samples_ophak;       dirs.samples    = cd;
cd ..;  cd fixDet;              dirs.fix        = cd;
cd ..;  cd behav_results;       dirs.behav      = cd;
cd ..;  cd msgs_ophak;          dirs.msgsO      = cd;
cd ..;  cd mat;                 dirs.mat        = cd;
cd ..;
cd ..;
cd function_library;            dirs.funclib    = cd;
cd ..;
cd results;                     dirs.res        = cd;
        cd 'AOImasks';          dirs.AOImasks   = cd;
cd ..;
cd(dirs.home);
addpath(genpath(dirs.funclib));                 % add dirs to path
addpath('function_library_cus');


%%% get all trials, parse into subject and stimulus
[files,nfiles]  = FileFromFolder(dirs.fix,[],'mat');
files           = parseFileNames(files);

% load result table from main script
res_Table = load([dirs.behav,'\', FileFromFolder(dirs.behav,[],'mat').name]).results;

fhndl   = -1;
lastRead= '';
sessionFileName = sprintf('%s.mat',files(1).subj);
expt = load(fullfile(dirs.mat,sessionFileName),'expt').expt;
geom = load(fullfile(dirs.mat,sessionFileName),'geometry').geometry.displayArea;
if nfiles~=height(res_Table)
    error('some trials are losted')
end
for p=1:nfiles
    % load fix data
    dat  = load(fullfile(dirs.fix,[files(p).fname '.mat'])); dat = dat.dat;
    if isempty(dat.time)
        warning('no data for %s, empty file',files(p).fname);
        continue;
    end
    
    % get msgs
    msgs    = loadMsgs(fullfile(dirs.msgsO,[files(p).fname '.txt']));
    [times,what,msgs] = parseMsgs_Form(msgs);

%     if ~strcmp(lastRead,sessionFileName)
%         lastRead = sessionFileName;
%         fInfo = [sess.expt.stim.fInfo];
%     end
%     qWhich= strcmp({fInfo.name},what{1});
    
    % load img, if only one
%     if ~~exist(fullfile(dirs.AOImasks,what{1}),'file')
%         img.data = imread(fullfile(dirs.AOImasks,what{1}));
%     elseif ~~exist(fullfile(sess.expt.stim(qWhich).fInfo.folder,what{1}),'file')
%         img.data = imread(fullfile(sess.expt.stim(qWhich).fInfo.folder,what{1}));
%     else
%     grayImg = uint8(sess.expt.stim*255);
    ut = UT(geom.width/10, expt.winRect(3), res_Table.headDist(p));
    bgCenter = expt.winRect(3:4)/2;
    tgCenter_ = [res_Table.ECC(p),res_Table.Orient(p)];
    tgCenter = ut.Pol2Rect(tgCenter_);
    
    grayImg = genStim(expt.winRect, ut, res_Table.bgContrast(p), ...
        1, tgCenter, expt.GaborSF, expt.GaborWidth, ...
        expt.GaborOrient, expt.bgWidth, res_Table.seed(p));
    img.data     = uint8 (cat(3,grayImg,grayImg,grayImg)*255);
    
    if ~isempty(img)
        % get position on screen
        stimRect = [0,0,size(img.data,2),size(img.data,1)];
        img.x    = linspace(stimRect(1),stimRect(3),size(img.data,2));
        img.y    = linspace(stimRect(2),stimRect(4),size(img.data,1));
    end
    
    % plot
    if ~ishghandle(fhndl)
        fhndl = figure('Units','normalized','Position',[0 0 1 1]);  % make fullscreen figure
    else
        figure(fhndl);
        clf;
    end
    set(fhndl,'Visible','on');  % assert visibility to bring window to front again after keypress
    drawFix(dat,dat.fix,[dat.I2MCopt.xres dat.I2MCopt.yres],img,[dat.I2MCopt.missingx dat.I2MCopt.missingy],sprintf('subj %s, trial %03d, stim: %s',files(p).subj,files(p).runnr,what{1}));
    pause
    if ~ishghandle(fhndl)
        return;
    end
end
if ishghandle(fhndl)
    close(fhndl);
end

rmpath(genpath(dirs.funclib));                  % cleanup path
