% 2026-03-27 by GitHub Copilot
% Merge PNGs in a folder into a sequential MP4 video.

clear variables; clear global; fclose('all'); clc

rootDir = fileparts(mfilename('fullpath'));
addpath(fullfile(rootDir, '..', 'function_library'));

%% ========================= User Config =========================
% You can set one or more folders to convert.
inputFolders = {
    fullfile(rootDir, 'Results', 'Full360_time_windows', 'v2', 'win_0500ms_step_0100ms')
    fullfile(rootDir, 'Results', 'Full360_time_windows', 'v2', 'win_1000ms_step_0100ms')
};

frameRate = 10;
videoName = 'angle_scan_evolution.mp4';

%% ========================= Merge =========================
for iFolder = 1:numel(inputFolders)
    inDir = inputFolders{iFolder};
    if ~exist(inDir, 'dir')
        warning('Folder not found, skipped: %s', inDir);
        continue;
    end

    outVideo = fullfile(inDir, videoName);
    local_png_folder_to_video(inDir, outVideo, frameRate);
end

fprintf('Done.\n');

%% ========================= Local Function =========================
function local_png_folder_to_video(inDir, outVideo, frameRate)
    files = dir(fullfile(inDir, '*.png'));
    if isempty(files)
        warning('No PNG files found in folder: %s', inDir);
        return;
    end

    names = {files.name};
    if exist('natsortfiles', 'file') == 2
        names = natsortfiles(names);
    else
        names = sort(names);
    end

    firstImg = imread(fullfile(inDir, names{1}));
    if ismatrix(firstImg)
        firstImg = repmat(firstImg, [1 1 3]);
    end
    targetH = size(firstImg, 1);
    targetW = size(firstImg, 2);

    vw = VideoWriter(outVideo, 'MPEG-4');
    vw.FrameRate = frameRate;
    open(vw);

    for i = 1:numel(names)
        img = imread(fullfile(inDir, names{i}));
        if ismatrix(img)
            img = repmat(img, [1 1 3]);
        end
        if size(img, 1) ~= targetH || size(img, 2) ~= targetW
            img = imresize(img, [targetH, targetW]);
        end
        writeVideo(vw, img);
    end

    close(vw);
    fprintf('Video saved: %s (frames=%d, fps=%g)\n', outVideo, numel(names), frameRate);
end
