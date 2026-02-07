%% 合并分段的实验数据文件
% 该脚本合并实验数据的两个session，完整的480trial实验
% 支持对Ses2进行截断处理（保留部分trials）
% 处理逻辑：Ses1 + Ses2(可选截断)
% 最后把原始两个分段文件移到split子文件夹

clear; clc;

%% 参数设置
baseDir = 'E:\Desktop\临时文件\AttenSamp\FreeView\FreeView_Paper_sup\Processed_data';
subID = 'Sub2';  % 修改被试号
ses1_file = fullfile(baseDir, ['Dat_' subID '_Ses1.mat']);
ses2_file = fullfile(baseDir, ['Dat_' subID '_Ses2.mat']);
output_file = fullfile(baseDir, ['Dat_' subID '_Ses1.mat']);  % 输出覆盖Ses1的位置
split_folder = fullfile(baseDir, 'split');

% *** 截断参数：设置是否对Ses2进行截断 ***
% 设置为 -1 或 [] 表示不截断（全部保留）
% 设置为正整数表示只保留前N行
ses2_truncate_rows = 420;  % Sub2的情况：只保留前420行（舍弃后60行）

%% 检查文件是否存在
if ~isfile(ses1_file)
    error('Ses1文件不存在: %s', ses1_file);
end
if ~isfile(ses2_file)
    error('Ses2文件不存在: %s', ses2_file);
end

%% 加载两个文件
fprintf('加载Ses1数据...\n');
load(ses1_file, 'data', 'expT', 'expt', 'geometry', 'settings', 'systemInfo');
data_ses1 = data;
expT_ses1 = expT;
expt_ses1 = expt;
geometry_ses1 = geometry;
settings_ses1 = settings;
systemInfo_ses1 = systemInfo;

fprintf('加载Ses2数据...\n');
load(ses2_file, 'data', 'expT', 'expt', 'geometry', 'settings', 'systemInfo');
data_ses2 = data;
expT_ses2 = expT;
expt_ses2 = expt;
geometry_ses2 = geometry;
settings_ses2 = settings;
systemInfo_ses2 = systemInfo;

%% 对Ses2进行截断处理（如果需要）
if ~isempty(ses2_truncate_rows) && ses2_truncate_rows > 0 && ses2_truncate_rows < height(expT_ses2)
    fprintf('对Ses2进行截断处理: %d行 -> %d行\n', height(expT_ses2), ses2_truncate_rows);
    expT_ses2 = expT_ses2(1:ses2_truncate_rows, :);
    
    % 截断gaze数据
    if isfield(data_ses2.gaze, 'left') && isstruct(data_ses2.gaze.left)
        fields_left = fieldnames(data_ses2.gaze.left);
        for i = 1:length(fields_left)
            field = fields_left{i};
            if isnumeric(data_ses2.gaze.left.(field)) && size(data_ses2.gaze.left.(field), 1) > 0
                data_ses2.gaze.left.(field) = data_ses2.gaze.left.(field)(1:ses2_truncate_rows, :);
            end
        end
    end
    
    if isfield(data_ses2.gaze, 'right') && isstruct(data_ses2.gaze.right)
        fields_right = fieldnames(data_ses2.gaze.right);
        for i = 1:length(fields_right)
            field = fields_right{i};
            if isnumeric(data_ses2.gaze.right.(field)) && size(data_ses2.gaze.right.(field), 1) > 0
                data_ses2.gaze.right.(field) = data_ses2.gaze.right.(field)(1:ses2_truncate_rows, :);
            end
        end
    end
    
    if isfield(data_ses2.gaze, 'systemTimeStamp')
        data_ses2.gaze.systemTimeStamp = data_ses2.gaze.systemTimeStamp(1:ses2_truncate_rows);
    end
end

%% 合并expT表（直接拼接）
fprintf('合并expT表...\n');
expT_merged = [expT_ses1; expT_ses2];
fprintf('  原Ses1: %d行，原Ses2: %d行，合并后: %d行\n', ...
    height(expT_ses1), height(expT_ses2), height(expT_merged));

%% 合并data.gaze
fprintf('合并gaze数据...\n');
data_merged = data_ses1;  % 使用Ses1的结构作为基础

% 合并left和right eye数据
if isfield(data_ses1.gaze, 'left') && isfield(data_ses2.gaze, 'left')
    % 合并left眼数据
    fields_left = fieldnames(data_ses1.gaze.left);
    for i = 1:length(fields_left)
        field = fields_left{i};
        data_left_1 = data_ses1.gaze.left.(field);
        data_left_2 = data_ses2.gaze.left.(field);
        if isnumeric(data_left_1) && isnumeric(data_left_2)
            data_merged.gaze.left.(field) = [data_left_1; data_left_2];
        end
    end
end

if isfield(data_ses1.gaze, 'right') && isfield(data_ses2.gaze, 'right')
    % 合并right眼数据
    fields_right = fieldnames(data_ses1.gaze.right);
    for i = 1:length(fields_right)
        field = fields_right{i};
        data_right_1 = data_ses1.gaze.right.(field);
        data_right_2 = data_ses2.gaze.right.(field);
        if isnumeric(data_right_1) && isnumeric(data_right_2)
            data_merged.gaze.right.(field) = [data_right_1; data_right_2];
        end
    end
end

% 合并systemTimeStamp
if isfield(data_ses1.gaze, 'systemTimeStamp') && isfield(data_ses2.gaze, 'systemTimeStamp')
    data_merged.gaze.systemTimeStamp = [data_ses1.gaze.systemTimeStamp; data_ses2.gaze.systemTimeStamp];
end

%% 使用第一个session的其他元数据
data = data_merged;
expt = expt_ses1;
geometry = geometry_ses1;
settings = settings_ses1;
systemInfo = systemInfo_ses1;

%% 创建split文件夹
if ~isfolder(split_folder)
    mkdir(split_folder);
    fprintf('创建split文件夹: %s\n', split_folder);
end

%% 移动原始文件到split文件夹
fprintf('移动原始文件到split文件夹...\n');
try
    movefile(ses1_file, fullfile(split_folder, ['Dat_' subID '_Ses1.mat']));
    fprintf('  已移动: Dat_%s_Ses1.mat\n', subID);
    
    movefile(ses2_file, fullfile(split_folder, ['Dat_' subID '_Ses2.mat']));
    fprintf('  已移动: Dat_%s_Ses2.mat\n', subID);
catch ME
    warning('移动文件时出错: %s', ME.message);
end

%% 保存合并后的文件
fprintf('保存合并后的文件: %s\n', output_file);
expT = expT_merged;  % 重命名为原始变量名
save(output_file, 'data', 'expT', 'expt', 'geometry', 'settings', 'systemInfo', '-v7.3');

fprintf('\n========== 合并完成 ==========\n');
fprintf('合并后的数据：\n');
fprintf('  - expT: %d行 x %d列\n', height(expT), width(expT));
fprintf('  - data.gaze.systemTimeStamp: %d个样本\n', length(data_merged.gaze.systemTimeStamp));
fprintf('  - 输出文件: %s\n', output_file);
fprintf('  - 备份文件夹: %s\n', split_folder);
