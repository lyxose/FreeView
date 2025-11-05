function [resfiles, sub_ses_res, select_sess] = get_eye_data_files(dirs_fix)
    % 数据文件获取函数：从指定目录筛选有效的被试-会话文件列表
    % 获取所有.mat文件
    [files,~] = FileFromFolder(dirs_fix,[],'mat');
    % 正则筛选文件名格式：Dat_SubXX_SesYY.mat
    filtstr = '^Dat_Sub(\d+)_Ses(\d+).mat$';
    matched = regexpi({files.name}.',filtstr,'tokens');
    resfile_idx = ~cellfun(@isempty,matched);
    resfiles   = files(resfile_idx);
    % 提取被试ID和会话ID
    sub_ses_res = cellfun(@(x) [str2double(x{1}{1}), str2double(x{1}{2})], matched(resfile_idx), 'UniformOutput', false);
    sub_ses_res = cell2mat(sub_ses_res(~cellfun(@isempty, matched(resfile_idx))));
    % 按 subID 升序排序
    [~, sortIdx] = sort(sub_ses_res(:,1));
    sub_ses_res = sub_ses_res(sortIdx,:);
    resfiles = resfiles(sortIdx);
    % 选择会话1（可根据需要修改筛选条件）
    select_sess = find((sub_ses_res(:,2)==1)==1)';
end