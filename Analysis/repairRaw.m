clear; clc
% 设置要遍历的文件夹路径
folderPath = './data/mat';  % 请替换为你实际的文件夹路径
%%% get eye tracker files 
% filter so we only get data that matches the filter. uses regexp
[files,~] = FileFromFolder(folderPath,[],'mat');
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
%%% get exp record
filtstr = '^EXP_Sub(\d+)_Ses(\d+)[_a-zA-Z]+\d+T\d+\.mat$';
matched = regexpi({files.name}.',filtstr,'tokens');
expfile_idx = ~cellfun(@isempty,matched);
expfiles   = files(expfile_idx);
sub_ses_exp = cellfun(@(x) [str2double(x{1}{1}), str2double(x{1}{2})], matched(expfile_idx), 'UniformOutput', false);
sub_ses_exp = cell2mat(sub_ses_exp(~cellfun(@isempty, matched(expfile_idx))));
if ~isequal(sub_ses_exp,sub_ses)
    unmatched_sub_ses = setdiff(sub_ses, sub_ses_exp, 'rows');
    warning('Files not matched... Sub%.0f_Ses%.0f', unmatched_sub_ses(1), unmatched_sub_ses(2))
end

%%
for p=1:nsess
    fprintf('Subject: %.0f, Session-%.0f\n', sub_ses(p,1), sub_ses(p,2))
    datPath = fullfile(folderPath, eyefiles(p).name);
    resPath = fullfile(folderPath, resfiles(p).name);
    expPath = fullfile(folderPath, expfiles(p).name);

    % 加载 Dat*.mat 文件
    dat  = load(datPath);
    resT = load(resPath);
    exp  = load(expPath);
    
    % 检查 expt.winRect 是否存在
    if ~isfield(dat, 'expt') || ~isfield(dat.expt, 'winRect')
        % 如果不存在，添加 expt.winRect 字段并赋值
        disp(['Adding expt.winRect to ' eyefiles(p).name]);
        % 如果 expt 字段不存在，先创建 expt
        if ~isfield(dat, 'expt')
            dat.expt = struct();
        end
        % 添加 winRect 字段并赋值
        dat.expt.winRect = exp.winRect;
    else
        disp([eyefiles(p).name ' already contains expt.winRect']);
    end

    % 检查 expt.restTime 是否存在
    if  ~isfield(dat.expt, 'restTime')
        disp(['Adding expt.restTime to ' eyefiles(p).name]);
        dat.expt.restTime = exp.restTime;
    else
        disp([eyefiles(p).name ' already contains expt.restTime']);
    end
    % 检查 expt.blockSize 是否存在
    if  ~isfield(dat.expt, 'blockSize')
        disp(['Adding expt.blockSize to ' eyefiles(p).name]);
        dat.expt.blockSize = exp.blockSize;
    else
        disp([eyefiles(p).name ' already contains expt.blockSize']);
    end
    % 检查 expt.bgWidth 是否存在
    if  ~isfield(dat.expt, 'bgWidth')
        disp(['Adding expt.bgWidth to ' eyefiles(p).name]);
        dat.expt.bgWidth = exp.bgWidth;
    else
        disp([eyefiles(p).name ' already contains expt.bgWidth']);
    end
    % 检查 expt.GaborSF 是否存在
    if  ~isfield(dat.expt, 'GaborSF')
        disp(['Adding expt.GaborSF to ' eyefiles(p).name]);
        dat.expt.GaborSF = exp.GaborSF;
    else
        disp([eyefiles(p).name ' already contains expt.GaborSF']);
    end
    % 检查 expt.GaborCyc 是否存在
    if  ~isfield(dat.expt, 'GaborCyc')
        disp(['Adding expt.GaborCyc to ' eyefiles(p).name]);
        dat.expt.GaborCyc = exp.GaborCyc;
    else
        disp([eyefiles(p).name ' already contains expt.GaborCyc']);
    end
    % 检查 expt.GaborWidth 是否存在
    if  ~isfield(dat.expt, 'GaborWidth')
        disp(['Adding expt.GaborWidth to ' eyefiles(p).name]);
        dat.expt.GaborWidth = exp.GaborWidth;
    else
        disp([eyefiles(p).name ' already contains expt.GaborWidth']);
    end
    % 检查 expt.GaborOrient 是否存在
    if  ~isfield(dat.expt, 'GaborOrient')
        disp(['Adding expt.GaborOrient to ' eyefiles(p).name]);
        dat.expt.GaborOrient = exp.GaborOrient;
    else
        disp([eyefiles(p).name ' already contains expt.GaborOrient']);
    end
    % 检查 expt.difficulty 是否存在
    if ~isfield(exp, 'difficulty')
        exp.difficulty = exp.tgContrast/exp.threshold;
    end
    if  ~isfield(dat.expt, 'difficulty') 
        disp(['Adding expt.difficulty to ' eyefiles(p).name]);
        dat.expt.difficulty = exp.difficulty;
    else
        disp([eyefiles(p).name ' already contains expt.difficulty']);
    end
    % 保存更新后的数据到原文件
    save(datPath, '-struct', 'dat');
    
    % 补全results table
    if ~any(strcmp(resT.results.Properties.VariableNames, 'headDist'))
        % 如果没有 headDist 列，则将其填充为 exp.ut.distance 的值
        disp(['Adding headDist to ' resfiles(p).name]);
        % 获取 exp.ut.distance 的值
        distanceValue = exp.ut.distance;
        % 填充 headDist 列
        resT.results.headDist = repmat(distanceValue, height(resT.results), 1);
    else
        disp([resfiles(p).name ' already contains headDist.']);
    end
    save(resPath,'-struct', 'resT')
end

