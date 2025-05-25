%% 表格合并与数据预处理

myDir = fileparts(mfilename('fullpath'));
cd(myDir);
                                dirs.home       = cd;
cd data;                        dirs.data       = cd;
        cd samples_ophak;       dirs.samples    = cd;
cd ..;  cd fixDet;              dirs.fix        = cd;
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

binSize = 25; % in pixel
dotsize = 10;

drawTraj = false;
drawCorrDist = true;
% only plot the selected fixations
Start_nFix = 2;
End_nFix = 2; 
% for Start_nFix = 2:3
% for End_nFix = Start_nFix:3


%%% get result table with eye data
[files,~] = FileFromFolder(dirs.fix,[],'mat');
filtstr = '^Dat_Sub(\d+)_Ses(\d+).mat$';
matched = regexpi({files.name}.',filtstr,'tokens');
resfile_idx = ~cellfun(@isempty,matched);
resfiles   = files(resfile_idx);
sub_ses_res = cellfun(@(x) [str2double(x{1}{1}), str2double(x{1}{2})], matched(resfile_idx), 'UniformOutput', false);
sub_ses_res = cell2mat(sub_ses_res(~cellfun(@isempty, matched(resfile_idx))));

select_sess = 1:length(resfiles);
% 初始化总表
combinedResT = table();


% 循环读取并合并表格
for p = select_sess
    % 加载当前会话数据
    currentResT = load(fullfile(dirs.fix, [resfiles(p).fname '.mat'])).expT;
    if size(currentResT, 2) == 16
        % 获取原列名
        originalVars = currentResT.Properties.VariableNames;
        
        % 插入空列（默认追加到末尾）
        currentResT.seed = zeros(size(currentResT, 1),1)-1; % 先添加空列到末尾
        
        % 调整列顺序：将新列插入原第15列之后
        currentResT = movevars(currentResT, 'seed', 'After', originalVars{15}); 
    end
    % 检查表格结构一致性（可选）
    if ~isempty(combinedResT) && ~isequal(currentResT.Properties.VariableNames, combinedResT.Properties.VariableNames)
        currentResT = movevars(currentResT, 'seed', 'After', originalVars{15}); 
    end
    if ~isempty(combinedResT) && ~isequal(currentResT.Properties.VariableNames, combinedResT.Properties.VariableNames)
        error('表格结构不一致，请检查变量名是否相同');
    end
    
    % 垂直拼接表格[4,7](@ref)
    combinedResT = vertcat(combinedResT, currentResT);
end

% 提取RT列数据
rtData = combinedResT.key1RT;
meanRT = mean(rtData);
stdRT = std(rtData);
rtData(abs(rtData - meanRT) > 3*stdRT) = NaN; % 3σ原则
rtData = rmmissing(rtData);

%% 分布统计与可视化
% 创建绘图窗口
figure('Position', [100 100 1200 500])

% 绘制直方图[13,14](@ref)
subplot(1,2,1)
h = histogram(rtData, 'Normalization', 'probability',...
    'BinWidth', 0.1, 'FaceColor', [0.2 0.4 0.8], 'EdgeColor', 'w');
title('RT分布直方图')
xlabel('反应时间（s）')
ylabel('概率密度')
grid on

% 添加统计标注
text(0.7, 0.8,...
    sprintf('均值 = %.1f ± %.1f s\n样本数 = %d',...
    mean(rtData), std(rtData), numel(rtData)),...
    'Units', 'normalized', 'FontSize', 10);

% 绘制概率密度曲线[12](@ref)
subplot(1,2,2)
pd = fitdist(rtData, 'Normal'); % 正态分布拟合[9](@ref)
x = linspace(min(rtData), max(rtData), 1000);
cdf_values = cdf(pd, x);
plot(x, cdf_values, 'LineWidth', 2, 'Color', [0.8 0.2 0.2])
title('累计概率密度曲线')
xlabel('反应时间（s）')
ylabel('概率密度')
grid on

% 保存图形
% saveas(gcf, fullfile(distpath, 'RT_distribution_combined.png'));