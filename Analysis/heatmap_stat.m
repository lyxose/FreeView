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
cd ..;                                  dirs.home               = cd;
    cd Analysis;                            dirs.ana            = cd;
            cd Processed_data;                  dirs.data       = cd;
                    cd samples_ophak;           dirs.samples    = cd;
            cd ..;  cd fixDet;                  dirs.fix        = cd;
            cd ..;  cd msgs_ophak;              dirs.msgsO      = cd;
            cd ..;  
    cd ..;  cd function_library;        dirs.funclib            = cd;  
    cd ..;  cd function_library_cus;    dirs.funclib_cus        = cd;
    cd ..;
cd ..;
    cd Data;                            dirs.mat                = cd;
 
cd(dirs.ana);
addpath(genpath(dirs.funclib));                 % add dirs to path
addpath(genpath(dirs.funclib_cus));

binSize = 25; % in pixel
dotsize = 10;

skip_corr=false;% ignore correct trials
drawTraj = false;
drawEach = false;
    tagNum = true;
drawCorrDist = false; % show expected distribution
fitModel = false;
% only plot the selected fixations
Start_nFix  = 4;
End_nFix    = 4; 
atLeast_nFix = 7;%End_nFix+5;
angbinSize = 15;
learn_stage_n = 100;
% for Start_nFix = 2:3
% for End_nFix = Start_nFix:3

R_max = 7.5;
R_min = 0;

%%% get result table with eye data
[files,~] = FileFromFolder(dirs.fix,[],'mat');
filtstr = '^Dat_Sub(\d+)_Ses(\d+).mat$';
matched = regexpi({files.name}.',filtstr,'tokens');
resfile_idx = ~cellfun(@isempty,matched);
resfiles   = files(resfile_idx);
sub_ses_res = cellfun(@(x) [str2double(x{1}{1}), str2double(x{1}{2})], matched(resfile_idx), 'UniformOutput', false);
sub_ses_res = cell2mat(sub_ses_res(~cellfun(@isempty, matched(resfile_idx))));

% select_sess = 1:length(resfiles);
% select_sess = 5:7;
exclude_sess = [];
% exclude_sess = [1,2,3,10,18,25];
% select_sess = 10:14;
% select_sess = [1:10 14:length(resfiles)]; % single subj
select_sess = find((sub_ses_res(:,2)==1)==1)';
% select_sess = find((sub_ses_res(:,1)>=12)==1)';
% select_sess = [2,3,4,5,7,8,10,11,16,17];

alpha = 0.05; % traj

figure(1);  % 在循环外创建 figure，固定窗口
% figure(2);  % 在循环外创建 figure，固定窗口

% heatpath = sprintf('./results/comb_heatmap-%.0f-%.0fto%.0f/',length(select_sess),Start_nFix, End_nFix);
heatpath = sprintf('./results/comb_heatmap/');
trajpath = sprintf('./results/comb_trajec-%.0f-%.0fto%.0f/',length(select_sess),Start_nFix, End_nFix);
corrpath = sprintf('./results/comb_corr-%.0f-%.0fto%.0f/',length(select_sess),Start_nFix, End_nFix);
distpath = sprintf('./results/comb_dist-%.0f-%.0fto%.0f/',length(select_sess),Start_nFix, End_nFix);

if ~exist(heatpath,'dir')
    mkdir(heatpath);
end
% if ~exist(trajpath,'dir')
%     mkdir(trajpath);
% end
% if ~exist(corrpath,'dir')
%     mkdir(corrpath);
% end
% if ~exist(distpath,'dir')
%     mkdir(distpath);
% end
fixPos = [];
stTs = [];   % color will be depended on the start time 
triSpl = []; % to split each trial by the end index
tgLoc = [];
% 
% for ecc = [2,4,6]
% for ori = 0:45:315
%% 提取和合并数据
for p=select_sess
    psp(1) = length(stTs)+1;
    if ismember(p,exclude_sess)
        continue
    end
    sess = load(fullfile(dirs.mat,sprintf("Dat_Sub%.0f_Ses%.0f.mat",sub_ses_res(p,1),sub_ses_res(p,2))),'expt','geometry');
    resT = load(fullfile(dirs.fix,[resfiles(p).fname '.mat'])).expT;
    try
        img_width = sess.expt.winRect(3);  % 二维空间的宽度
        img_height = sess.expt.winRect(4); % 二维空间的高度
    catch
        img_width=1920;
        img_height=1080;
    end
    ut = UT(sess.geometry.displayArea.width/10, img_width, mean(resT.headDist,1, 'omitnan'), false);
% for ori
%         rows = all(bsxfun(@eq, [resT.ECC, resT.Orient], [ecc,ori]), 2); 
%     for k = 1:max(ic)
%         fixPos = [];
%         rows = ~isnan(resT.judge);
%         iresT = resT(rows,:); % subset of resT, all in the EO ECC-Ori condition 
    iresT = resT;
    for i=1:height(iresT)
        if i < learn_stage_n
            continue
        end
        if skip_corr && ~isnan(iresT.key2RT(i)) || isempty(iresT.dat(i).fix)
            continue
        end
        tFixPos = transpose([iresT.dat(i).fix.xpos; iresT.dat(i).fix.ypos]); % this trial eye trajectory
        tTime = iresT.dat(i).fix.startT;  % start time stamp
        tFixPos = tFixPos(tTime>0,:);
        tTime = tTime(tTime>0);
        % 剔除一开始就在中心1°以内的trial，并跳过长度等于1的trial
        while length(tFixPos(:,1))>1 && norm(tFixPos(1,:)-[img_width, img_height]/2)<ut.deg2pix(1)
            if length(tFixPos(:,1))==1
                break
            end
            tFixPos = tFixPos(2:end,:); % remove the start fixation
            tTime = tTime(2:end,:);
        end
        
        % cut to early fixations
        % 按剔除后总fixation数量筛选trial，以免选取到最后几个已经位于目标附近的trial
        if length(tTime(:,1))>=atLeast_nFix  %
            % exclude some trial that 被试分心 or 眨眼 导致fixation信号开始晚于1s
            if tTime(1)>1000
%                 sprintf("Dat_Sub%.0f_Ses%.0f.mat",sub_ses_res(p,1),sub_ses_res(p,2))
%                 sprintf('i=%.0f',i)
%                 disp(iresT.dat(i).fix.startT);
                continue
            end
            if End_nFix == 0
                tFixPos = tFixPos(Start_nFix:end,:);
                tTime = tTime(Start_nFix:end,:);
            else
%                     tFixPos = tFixPos(Start_nFix+7:End_nFix+7,:)-tFixPos(Start_nFix:End_nFix,:)+ones(End_nFix-Start_nFix+1,2).*[img_width/2,img_height/2];                
                tFixPos = tFixPos(Start_nFix:End_nFix,:);
                tTime = tTime(Start_nFix:End_nFix,:);
            end
        else
            continue % drop short trials
        end

%             tFixPos = arrayfun(@(i) ut.Pol2Rect([1,-1].*(ut.Rect2Pol(tFixPos(i, :)-[img_width, img_height]/2))-[0,ori]).*[1,-1]+[img_width, img_height]/2, ...
%                         1:size(tFixPos, 1), 'UniformOutput', false);
%             tFixPos = cell2mat(tFixPos');
        % concate all the selected results
        tsp(1)  = length(stTs)+1;
        stTs = [stTs;tTime];
        fixPos = [fixPos; tFixPos];
        tsp(2) = length(stTs);  
        triSpl = [triSpl;tsp];
        ttgLoc = ut.deg2pix([iresT.Xtarg(i),iresT.Ytarg(i)]).*[1,-1]+[img_width,img_height]./2;
        tgLoc = [tgLoc;ttgLoc];
    end
    psp(2) = length(stTs);
    pSpl(p,:) = psp;
end
% end
% end
%% 统计角度范围(四分)
% 全局
fixPol = ut.Rect2Pol([1,-1].*(fixPos-[img_width,img_height]./2));
% fixPol(:,2) = mod(fixPos(:,2),360);
rangeCount = zeros(1,36000);
for i = 1:36000
    rangeCount(i) = sum(mod(fixPol(:,2)-i/100,360) < angbinSize ,"all");
end
figure()
plot((0.01:0.01:360)+angbinSize/2,rangeCount)
yLimits = ylim;
for o = 0:45:315
    line([o,o], yLimits, 'Color', 'r', 'LineStyle', '--');     % 红色虚线
    line([o,o]-22.5, yLimits, 'Color', 'b', 'LineStyle', '--');     % 红色虚线
end
title(sprintf('Fixation count %.0f to %.0f of %.0f+ in %.0f° bin  n = %.0f',Start_nFix,End_nFix,atLeast_nFix,angbinSize,length(select_sess)-length(exclude_sess)))

% 合并
rangeCount2 = zeros(1,450);
for i = 1:450
    rangeCount2(i) = sum(mod(fixPol(:,2)-i/10,45) < angbinSize ,"all");
end
figure()
plot((0.1:0.1:45)+angbinSize/2,rangeCount2)
yLimits = ylim;
line([22.5 22.5], yLimits, 'Color', 'b', 'LineStyle', '--'); % 蓝色虚线 
line([45 45], yLimits, 'Color', 'r', 'LineStyle', '--');     % 红色虚线
title(sprintf('Fixation count %.0f to %.0f of %.0f+ in %.0f° bin  n = %.0f',Start_nFix,End_nFix,atLeast_nFix,angbinSize,length(select_sess)-length(exclude_sess)))

%%
Fs = 100; % 采样频率（Hz）
N=length(rangeCount)*4;
f = (0:N-1)*(Fs/N); % 频率轴
data = rangeCount-mean(rangeCount);
X = fft([data,data,data,data]); % 快速傅里叶变换
X_mag = abs(X); % 计算幅度谱

% 由于对称性，只考虑前半部分
N_half = N/2;
f_half = f(2:N_half);
X_mag_half = X_mag(2:N_half);
figure()
lambda = 1./f_half;
window = lambda<360;
plot(lambda(window), X_mag_half(window));

%% 柱状图
n_bin = 16; % 自定义区间数量
shift = 360/n_bin/2; % 0度的方位正好在bin的中间
edges = linspace(0, 360, n_bin+1)-shift; % 生成0到360的n+1个等分边界
[counts, ~] = histcounts(mod(fixPol(:,2)+shift,360)-shift, edges); % counts为各区间元素数量
counts = counts - ([counts(end),counts(1:end-1)] + [counts(2:end),counts(1)])/2; % detrend

red = [1 0 0];              %   0               90              180             270
blue= [0 0 1];              %       22.5    67.5    112.5   157.5   202.5   247.5   292.5   337.5
pink = [1 0.7529 0.7961];   %           45              135             225             315    
color_map = repmat([red;blue;pink;blue], ceil(n_bin/4), 1); % 红、蓝、红、浅蓝交替重复
color_map = color_map(1:n_bin, :); % 截取与数据长度匹配的部分

figure()
bar((edges(1:end-1)+edges(2:end))/2,counts,'FaceColor', 'flat', 'CData', color_map)
% bar(((edges(1:end-1)+edges(2:end))/2-shift,counts)

%% group level
for st_type = 1:3
allCounts=[];
for p = select_sess
    pfixPol = fixPol(pSpl(p,1):pSpl(p,2),2);
    [pCounts, ~] = histcounts(mod(pfixPol+shift,360)-shift, edges); % counts为各区间元素数量
    pCounts = pCounts - ([pCounts(end),pCounts(1:end-1)] + [pCounts(2:end),pCounts(1)])/2; % detrend
    allCounts = [allCounts;pCounts];
end
xtags = (edges(1:end-1)+edges(2:end))/2;
color_map = repmat([red;blue;pink;blue], ceil(n_bin/4), 1); % 红、蓝、红、浅蓝交替重复

% % axis vs. gap
if st_type == 1
    allCounts = [sum(allCounts(:,1:2:16),2),sum(allCounts(:,2:2:16),2)];
    xtags = {'Axis','Gap'};
    diff_data = allCounts(:,1) - allCounts(:,2);
    [~, p_diff] = ttest(diff_data);
    color_map = color_map(1:2,:);
% card vs. ordi
elseif st_type == 2
    allCounts = [sum(allCounts(:,1:4:16),2), sum(allCounts(:,3:4:16),2)];
    xtags = {'Card','Ordi'};
    diff_data = allCounts(:,1) - allCounts(:,2);
    [~, p_diff] = ttest(diff_data);
    color_map = color_map([1,3],:);
end

% 数据准备
mean_data = mean(allCounts, 1); % 计算均值
std_data = std(allCounts, 0, 1); % 计算标准差
n_samples = size(allCounts, 1); % 样本数量

% 绘制柱状图
figure()
hBar = bar(mean_data,'FaceColor', 'flat', 'CData', color_map);
if st_type==3
    xticks(1:length(mean_data));
end
xticklabels(xtags);      
hold on;

% 获取柱子位置并添加误差线
xPos = hBar.XEndPoints; % 动态获取柱子中心坐标
errorbar(xPos, mean_data, std_data,...
'LineStyle', 'none',...
'Color', [0.2 0.2 0.2],...
'LineWidth', 1.2); % 绘制误差线

% 单样本t检验
significance = zeros(1, length(mean_data));
for i = 1:length(mean_data)
[~, p] = ttest(allCounts(:,i), 0); % 单样本t检验
if p < 0.001
sigSymbol = '​​***';
elseif p < 0.01
sigSymbol = '**​​';
elseif p < 0.05
sigSymbol = '*';
else
sigSymbol = 'n.s.';
end
significance(i) = p;

% 标注显著性符号
yPos = mean_data(i) + 1.2*std_data(i);  % 动态调整标注高度
text(xPos(i), yPos, sigSymbol,...
    'HorizontalAlignment', 'center',...
    'VerticalAlignment', 'bottom',...
    'FontSize', 12,...
    'FontWeight', 'bold');  % 显著性标注[1,2](@ref)
end

% 添加参考线（0值线）
yline(0, '--', 'LineWidth', 1.2, 'Color', [0.5 0.5 0.5]);

% 配对检验
if length(xtags)==2
    text_str = sprintf('p = %.3f', p_diff);

    % 在归一化坐标系下添加文本
    text(0.05, 0.95, text_str,... 
        'Units', 'normalized',...      % 使用归一化坐标系[7,8](@ref)
        'HorizontalAlignment', 'left',...  % 左对齐[4,5](@ref)
        'VerticalAlignment', 'top',...     % 顶部对齐[4,5](@ref)
        'FontSize', 12,...                % 字号控制[2,5](@ref)
        'FontWeight', 'bold',...          % 字体加粗[2](@ref)
        'Color', 'k');                   % 黑色字体[2](@ref)
end

% 图形美化
set(gca, 'FontSize', 12, 'LineWidth', 1.2)
xlabel('Orientation', 'FontSize', 14)
ylabel('Detrended Count', 'FontSize', 14)
title(sprintf('Fixation heatmap %.0f to %.0f of %.0f+ n = %.0f'             ,Start_nFix,End_nFix,atLeast_nFix,length(select_sess)-length(exclude_sess)), 'FontSize', 16)
box off
end
%% 轴vs非轴
% allCounts=[];
% for p = select_sess
%     pfixPol = fixPol(pSpl(p,1):pSpl(p,2),2);
%     [pCounts, ~] = histcounts(mod(pfixPol+shift,360)-shift, edges); % counts为各区间元素数量
%     pCounts = pCounts - ([pCounts(end),pCounts(1:end-1)] + [pCounts(2:end),pCounts(1)])/2; % detrend
%     allCounts = [allCounts;pCounts];
% end
% 
% allCounts = [sum(allCounts(:,[1,5]),2), sum(allCounts(:,[3,7]),2)];
% xtags = {'Card','Ordi'};
% 
% % 配色方案（红蓝对比）[8](@ref)
% color_map = [0.8 0.2 0.2; 0.2 0.4 0.8]; 
% 
% % 数据统计
% mean_data = mean(allCounts, 1);
% std_data = std(allCounts, 0, 1);
% n_samples = size(allCounts, 1);
% 
% % 创建图形容器
% figure('Position', [100 100 800 600])
% hBar = bar(mean_data, 'FaceColor', 'flat', 'CData', color_map,...
%     'BarWidth', 0.7);
% xticklabels(xtags);
% hold on;
% 
% % 误差线绘制（95%置信区间）[3](@ref)
% xPos = hBar.XEndPoints;
% errorbar(xPos, mean_data, std_data/sqrt(n_samples)*1.96,... 
%     'LineStyle', 'none', 'Color', [0.3 0.3 0.3],...
%     'LineWidth', 1.5, 'CapSize', 15);
% 
% % 执行配对t检验
% diff_data = allCounts(:,1) - allCounts(:,2);
% [~, p] = ttest(diff_data);
% 
% % 动态计算标注参数
% y_max = max(mean_data) + max(std_data)*1.5;
% x_limits = [min(xPos), max(xPos)];
% 
% if p < 0.001
%     symbol = '***';
% elseif p < 0.01
%     symbol = '**';
% elseif p < 0.05
%     symbol = '*';
% else
%     symbol = '';
% end
% 
% % 显著性标注系统
% if p < 0.05
%     line(x_limits, [y_max y_max], 'Color', 'k', 'LineWidth', 1.5);
%     text(mean(x_limits), y_max*1.05, symbol,...
%         'FontSize', 16, 'HorizontalAlignment', 'center');
% end
% 
% % 图形美化
% set(gca, 'FontSize', 14, 'LineWidth', 1.5,...
%     'YTick', linspace(0, max(mean_data)*1.2, 6),...
%     'Box', 'off');
% ylabel('Normalized Response', 'FontWeight', 'bold', 'FontSize', 16)
% title('Paired Comparison with Significance', 'FontSize', 18)
% 
% 辅助函数
%% 统计角度范围(六分)
% % 全局
% fixPol = ut.Rect2Pol([1,-1].*(fixPos-[img_width,img_height]./2));
% % fixPol(:,2) = mod(fixPos(:,2),360);
% rangeCount = zeros(1,36000);
% for i = 1:36000
%     rangeCount(i) = sum(mod(fixPol(:,2)-i/100,360) < angbinSize ,"all");
% end
% figure()
% plot((0.01:0.01:360)+angbinSize/2,rangeCount)
% yLimits = ylim;
% for o = 30:60:360
%     line([o,o], yLimits, 'Color', 'r', 'LineStyle', '--');     % 红色虚线
%     line([o,o]-30, yLimits, 'Color', 'b', 'LineStyle', '--');     % 红色虚线
% end
% title(sprintf('Fixation count %.0f to %.0f of %.0f+ in %.0f° bin  n = %.0f',Start_nFix,End_nFix,atLeast_nFix,angbinSize,length(select_sess)-length(exclude_sess)))
% % 合并
% rangeCount2 = zeros(1,600);
% for i = 1:600
%     rangeCount2(i) = sum(mod(fixPol(:,2)-i/10,60) < angbinSize ,"all");
% end
% figure()
% plot((0.1:0.1:60)+angbinSize/2,rangeCount2)
% yLimits = ylim;
% line([60 60], yLimits, 'Color', 'b', 'LineStyle', '--'); % 蓝色虚线
% line([30 30], yLimits, 'Color', 'r', 'LineStyle', '--');     % 红色虚线
% title(sprintf('Fixation count %.0f to %.0f of %.0f+ in %.0f° bin  n = %.0f',Start_nFix,End_nFix,atLeast_nFix,angbinSize,length(select_sess)-length(exclude_sess)))
% 
% %%
% Fs = 10; % 采样频率（Hz）
% N=length(rangeCount);
% f = (0:N-1)*(Fs/N); % 频率轴
% 
% X = fft(rangeCount); % 快速傅里叶变换
% X_mag = abs(X); % 计算幅度谱
% 
% % 由于对称性，只考虑前半部分
% N_half = N/2;
% f_half = f(1:N_half);
% X_mag_half = X_mag(1:N_half);
% 
% plot(f_half, X_mag_half);

%% 绘制无序热图
        figure(1);
        clf;  
        edges_x = linspace(1, img_width,  round(img_width/binSize));  % 网格的X边界
        edges_y = linspace(1, img_height, round(img_height/binSize)); % 网格的Y边界
        pCounts = hist3(fixPos, 'Edges', {edges_x(1:end-1), edges_y(1:end-1)});

%         % statistics for the saccade distance
%         disPos = fixPos(2:end,:)-fixPos(1:end-1,:);
%         disPos(triSpl(:,2),:)=NaN;
%         disPos=disPos(~isnan(disPos(:,1)),:)+[img_width/2,img_height/2];
%         counts = hist3(disPos, 'Edges', {edges_x(1:end-1), edges_y(1:end-1)});
        density_map = pCounts' / sum(pCounts(:));
%         figure;
        imagesc(edges_x, edges_y, density_map);  % 使用imagesc绘制热图
        pr = prctile(density_map(:), 99.99);
        clim([min(density_map(:)), pr]);
%         set(gcf, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
        axis image;
        colorbar; % 显示颜色条
        % 标记中央位置：用 rectangle 函数绘制一个圆形
        hold on;  % 保持当前图形
        rectangle('Position', [round(img_width/2)-5, round(img_height/2)-5, 10, 10], ...
                  'EdgeColor', '#E3170D', 'LineWidth', 0.7, 'Curvature', [1, 1]);  % 
%         rectangle('Position', [target_loc(1)-tgWidth/2, target_loc(2)-tgWidth/2, tgWidth, tgWidth], ...
%                   'EdgeColor', 'r', 'LineWidth', 0.7, 'Curvature', [1, 1]);  %        
        % 标记目标位置
        mark_colors = [red;pink];
        for ecc = [2,4,6]
        for ori = 0:45:315
            tgWidth = ut.deg2pix(sess.expt.GaborWidth);
            target_loc = ut.Pol2Rect([ecc,ori]);
            target_loc = target_loc.*[1,-1]+[img_width, img_height]/2;
            
            rectangle('Position', [target_loc(1)-tgWidth/2, target_loc(2)-tgWidth/2, tgWidth, tgWidth], ...
                      'EdgeColor', mark_colors(mod(ori,90)/45+1,:), 'LineWidth', 0.7, 'Curvature', [1, 1]);
        end
        end
        % 标记fixation范围
%         tgECC_max = ut.deg2pix(R_max);
% %         tgECC_min = ut.deg2pix(R_min);
% %         rectangle('Position', [img_width/2-tgECC_min, img_height/2-tgECC_min, 2*tgECC_min, 2*tgECC_min], ...
% %                   'EdgeColor', 'w', 'LineWidth', 0.7, 'Curvature', [1, 1]);
%         rectangle('Position', [img_width/2-tgECC_max, img_height/2-tgECC_max, 2*tgECC_max, 2*tgECC_max], ...
%                   'EdgeColor', 'w', 'LineWidth', 0.7, 'Curvature', [1, 1]);

%         % 标记target范围：
%         % 绘制红色散点（tgLoc是1600x2的坐标矩阵，格式为[x, y]）
%         scatter(tgLoc(:,1), tgLoc(:,2), 5, 'r', 'filled', ...
%             'MarkerEdgeColor', 'r','MarkerEdgeAlpha', 0.3, 'MarkerFaceAlpha', 0.3); 
%         % 调整坐标轴范围与热图一致（避免散点超出热图显示区域）
%         xlim([min(edges_x), max(edges_x)]); 
%         ylim([min(edges_y), max(edges_y)]);
%         % 设置散点透明度（可选，增强热图与散点的对比）
%         set(gca, 'Layer', 'top');  % 确保散点显示在热图上方
        

        % 标记背景范围
        bgWidth = ut.deg2pix(sess.expt.bgWidth);
        rectangle('Position', [img_width/2-bgWidth/2, img_height/2-bgWidth/2, bgWidth, bgWidth], ...
                  'EdgeColor', '#E3170D', 'LineWidth', 0.7, 'Curvature', [1, 1]);
        title(sprintf('Fixation heatmap %.0f to %.0f of %.0f+ n = %.0f'             ,Start_nFix,End_nFix,atLeast_nFix,length(select_sess)-length(exclude_sess)));
        exportgraphics(gcf, sprintf([heatpath,'%.0fto%.0f_%.0fsess_heat.png'],Start_nFix,End_nFix,length(select_sess)-length(exclude_sess)), 'Resolution', 600);
        savefig(sprintf([heatpath,'%.0fto%.0f_%.0fsess_heat.fig']            ,Start_nFix,End_nFix,length(select_sess)-length(exclude_sess)))
        hold off; 
%%
        if drawTraj
            figure(2); clf; hold on; axis equal; 
            colormap(cool);  % 选择traj的颜色映射，例如 'parula', 'jet', 'hot'
            clim([0,1]);
            axis([0 img_width 0 img_height]);  % 横纵坐标范围为 img_width 和 img_height
            set(gca, 'YDir', 'reverse');  % y轴零点在左上角
            
            % 归一化时间值到 [0,1]，用于颜色映射
    %             normTime = time / maxTime;
    %         normTime = tTime / max(tTime);
            rankTime = tiedrank(stTs);
            normTime = rankTime/max(rankTime);
            
            % 绘制渐变线段/点
    %         for i = 1:size(triSpl,1)
    %             tsp = triSpl(i,:);
    %             xpos = fixPos(tsp(1):tsp(2),1);
    %             ypos = fixPos(tsp(1):tsp(2),2);
    %             tTime = normTime(tsp(1):tsp(2));
    %             for j = 1:length(xpos)%-1
    %                 ct = tTime(j);  % 当前时间的归一化值
    %                 color = cool(256);  % 获取 colormap 颜色
    %                 idx = max(1, round(ct * 255));  % 找到颜色索引
    %                 scatter(xpos(j), ypos(j), dotsize, color(idx, :), 'filled', 'MarkerEdgeAlpha', alpha, 'MarkerFaceAlpha', alpha);
    % %                 plot(xpos(j:j+1), ypos(j:j+1), '-', 'Color', [color(idx, :), alpha], 'LineWidth', 1.5);
    %             end
    %         end
            % 预计算所有点的颜色
            color = cool(256);  % 获取 colormap 颜色
            all_xpos = [];
            all_ypos = [];
            all_colors = [];
            
            for i = 1:size(triSpl, 1)
                tsp = triSpl(i, :);
                xpos = fixPos(tsp(1):tsp(2), 1);
                ypos = fixPos(tsp(1):tsp(2), 2);
                tTime = normTime(tsp(1):tsp(2));
                
                % 预计算每个点的颜色
                idx = max(1, round(tTime * 255));  % 找到颜色索引
                colors = color(idx, :);
                % 逐个检查轨迹
                if drawEach
                    if tagNum
                        for ii = 1:length(xpos)
                            text(xpos(ii), ypos(ii), num2str(ii), ...
                                 'FontSize', 5, 'Color', colors(ii,:));
                        end
                    else
                        scatter(xpos, ypos, dotsize, colors, 'filled');
                    end
                    rectangle('Position', [round(img_width/2)-5, round(img_height/2)-5, 10, 10], ...
                              'EdgeColor', 'k', 'LineWidth', 0.7, 'Curvature', [1, 1]);  % 红色圆形
                    rectangle('Position', [tgLoc(i,:)-5, 10, 10], ...
                              'EdgeColor', 'r', 'LineWidth', 0.7, 'Curvature', [1, 1]);  % 红色圆形
                    
                    input(sprintf('Next:%d',i+1))
                    cla;
                end
                
                % 收集所有点的坐标和颜色
                all_xpos = [all_xpos; xpos];
                all_ypos = [all_ypos; ypos];
                all_colors = [all_colors; colors];
            end
            
            % 一次性绘制所有点
            scatter(all_xpos, all_ypos, dotsize, all_colors, 'filled', 'MarkerEdgeAlpha', alpha, 'MarkerFaceAlpha', alpha);
    
            % 设置图像的坐标范围和 y 轴的零点在左上角
            % 标记中心位置
            rectangle('Position', [round(img_width/2)-5, round(img_height/2)-5, 10, 10], ...
                      'EdgeColor', 'k', 'LineWidth', 0.7, 'Curvature', [1, 1]);  % 红色圆形
            % 标记目标位置
            for ecc = [2,4,6]
            for ori = 0:45:315
                tgWidth = ut.deg2pix(sess.expt.GaborWidth);
                target_loc = ut.Pol2Rect([ecc,ori]);
                target_loc = target_loc.*[1,-1]+[img_width, img_height]/2;
        
                rectangle('Position', [target_loc(1)-tgWidth/2, target_loc(2)-tgWidth/2, tgWidth, tgWidth], ...
                          'EdgeColor', 'k', 'LineWidth', 0.7, 'Curvature', [1, 1]);
            end
            end
            % 标记背景范围
            bgWidth = ut.deg2pix(sess.expt.bgWidth);
            rectangle('Position', [img_width/2-bgWidth/2, img_height/2-bgWidth/2, bgWidth, bgWidth], ...
                      'EdgeColor', 'k', 'LineWidth', 0.7, 'Curvature', [1, 1]);
    
            c = colorbar;  % 显示颜色条
            qtls = prctile(stTs,round(c.Ticks*100));
            c.TickLabels = arrayfun(@(x) sprintf('%.2f', x), qtls, 'UniformOutput', false);
            c.Label.String = 'Start Time (ms)';
            title(sprintf('Fixation Trajectory %.0f to %.0f of %.0f+  %.0f subj',Start_nFix,End_nFix,atLeast_nFix,length(select_sess)));
    %         xlabel('X Position');
    %         ylabel('Y Position');
            exportgraphics(gcf, sprintf([trajpath,'%.0fto%.0f_%.0fsess_tra.png'],Start_nFix,End_nFix,length(select_sess)), 'Resolution', 600);
            savefig(sprintf([trajpath,'%.0fto%.0f_%.0fsess_tra.fig'],Start_nFix,End_nFix,length(select_sess)))
            hold off; 
        end
        
        if drawCorrDist
            % bin的中心坐标
            [ex,ey]=meshgrid((edges_x(1:end-1)+edges_x(2:end))/2,...
                             (edges_y(1:end-1)+edges_y(2:end))/2);
            ex = ex-img_width/2;
            ey = ey-img_height/2;
            ec = sqrt(ex.^2+ey.^2);
            cir_mask = ec<bgWidth/2 & ec>ut.deg2pix(1);
            obs_data=density_map(cir_mask);
            obs_data=obs_data./sum(obs_data(:));
            exp_data_uni = ones(length(obs_data),1)./length(obs_data);
            sigma_space = 1.5.^(-10:25);
            corr_4gau = zeros(length(sigma_space),1);
            dist_4gau = zeros(length(sigma_space),1);
            corr_2dgau = zeros(length(sigma_space),1);
            dist_2dgau = zeros(length(sigma_space),1);
            
            if fitModel          
            h = waitbar(0, 'Processing...'); % 初始进度为 0            
            for k = 1:length(sigma_space)
            
            sigma = sigma_space(k);
            hold on
%             exp_distr_4gau = gaumap(img_width, img_height, sigma);
            sigma_r = 200;
            exp_distr_4gau = OctagonalPolmap(img_width, img_height, sigma_r, sigma);
            bin_values_4gau = distr2bin(exp_distr_4gau, img_width, img_height, edges_x, edges_y);
            bin_values_4gau(~cir_mask) = 0;
            bin_values_4gau = bin_values_4gau./sum(bin_values_4gau(:));
            exp_data_4gau = bin_values_4gau(cir_mask);
            corr_4gau(k) = corr(exp_data_4gau,obs_data);
            dist_4gau(k) = norm(exp_data_4gau-obs_data);
            
            exp_distr_2dgau=  fspecial('gaussian', [img_height, img_width], sigma);
            bin_values_2dgau = distr2bin(exp_distr_2dgau, img_width, img_height, edges_x, edges_y);
            bin_values_2dgau(~cir_mask) = 0;
            bin_values_2dgau = bin_values_2dgau./sum(bin_values_2dgau(:));
            exp_data_2dgau = bin_values_2dgau(cir_mask);
            corr_2dgau(k) = corr(exp_data_2dgau,obs_data);
            dist_2dgau(k) = norm(exp_data_2dgau-obs_data);
            
            % imagesc(edges_x, edges_y, bin_values_4gau);
            % rectangle('Position', [img_width/2-bgWidth/2, img_height/2-bgWidth/2, bgWidth, bgWidth], ...
            %           'EdgeColor', '#E3170D', 'LineWidth', 0.7, 'Curvature', [1, 1]);
            % rectangle('Position', [img_width/2-5, img_height/2-5, 10, 10], ...
            %           'EdgeColor', '#E3170D', 'LineWidth', 0.7, 'Curvature', [1, 1]);
            % axis image;
            % colorbar; % 显示颜色条
            % hold off
            waitbar(k / length(sigma_space), h, sprintf('Processing... %d%%', round(k / length(sigma_space) * 100)));
            end
            
            figure(4); clf; hold on; 
            plot(sigma_space,corr_4gau , 'LineWidth', 1.5)
            plot(sigma_space,corr_2dgau, 'LineWidth', 1.5)
            plot(sigma_space,ones(k,1).*corr(exp_data_uni,obs_data), 'LineWidth', 1.5)
            legend('4gau','2dgau','uniform')
            % 设置横坐标为对数刻度
            set(gca, 'XScale', 'log');
            xlabel('sigma')
            ylabel('Correlation')
            title(sprintf('Heatmap correlation %.0f to %.0f of %.0f+  %.0f subj',Start_nFix,End_nFix,atLeast_nFix,length(select_sess)));
            exportgraphics(gcf, sprintf([corrpath,'%.0fto%.0f_%.0fsess_corr.png'],Start_nFix,End_nFix,length(select_sess)), 'Resolution', 600);
            hold off
            
            
            figure(5); clf; hold on; 
            plot(sigma_space,dist_4gau , 'LineWidth', 1.5)
            plot(sigma_space,dist_2dgau, 'LineWidth', 1.5)
            plot(sigma_space,ones(k,1).*norm(exp_data_uni-obs_data), 'LineWidth', 1.5)
            legend('4gau','2dgau','uniform')
            % 设置横坐标为对数刻度
            set(gca, 'XScale', 'log');
            xlabel('sigma')
            ylabel('Euclidean distance')
            title(sprintf('Heatmap distance %.0f to %.0f of %.0f+  %.0f subj',Start_nFix,End_nFix,atLeast_nFix,length(select_sess)));
            exportgraphics(gcf, sprintf([distpath,'%.0fto%.0f_%.0fsess_dist.png'],Start_nFix,End_nFix,length(select_sess)), 'Resolution', 600);
            hold off
            end

        end
 
% end
% end
% end
% end
if fitModel
    sigma=200;
    exp_distr_4gau = OctagonalPolmap(img_width, img_height, 100, 0.265, 100);
    bin_values_4gau = distr2bin(exp_distr_4gau, img_width, img_height, edges_x, edges_y);
    bin_values_4gau(~cir_mask) = 0;
    bin_values_4gau = bin_values_4gau./sum(bin_values_4gau(:));
    exp_data_4gau = bin_values_4gau(cir_mask);
    corr_4gau_best = corr(exp_data_4gau,obs_data)
    dist_4gau_best = norm(exp_data_4gau-obs_data)
    figure(11);
    imagesc(edges_x, edges_y, bin_values_4gau);
    rectangle('Position', [img_width/2-bgWidth/2, img_height/2-bgWidth/2, bgWidth, bgWidth], ...
              'EdgeColor', '#E3170D', 'LineWidth', 0.7, 'Curvature', [1, 1]);
    rectangle('Position', [img_width/2-5, img_height/2-5, 10, 10], ...
              'EdgeColor', '#E3170D', 'LineWidth', 0.7, 'Curvature', [1, 1]);
    % 标记目标位置
    for ecc = [2,4,6]
    for ori = 0:45:315
        tgWidth = ut.deg2pix(sess.expt.GaborWidth);
        target_loc = ut.Pol2Rect([ecc,ori]);
        target_loc = target_loc.*[1,-1]+[img_width, img_height]/2;
        rectangle('Position', [target_loc(1)-tgWidth/2, target_loc(2)-tgWidth/2, tgWidth, tgWidth], ...
                  'EdgeColor', '#E3170D', 'LineWidth', 0.7, 'Curvature', [1, 1]);
    end
    end
    axis image;
    title('Σ4 Gau')
    colorbar; % 显示颜色条
    
    corr_4gau(k) = corr(exp_data_4gau,obs_data);
    
    exp_distr_2dgau=  fspecial('gaussian', [img_height, img_width], sigma);
    bin_values_2dgau = distr2bin(exp_distr_2dgau, img_width, img_height, edges_x, edges_y);
    bin_values_2dgau(~cir_mask) = 0;
    bin_values_2dgau = bin_values_2dgau./sum(bin_values_2dgau(:));
    exp_data_2dgau = bin_values_2dgau(cir_mask);
    corr_2dgau(k) = corr(exp_data_2dgau,obs_data);
    dist_2dgau(k) = norm(exp_data_2dgau-obs_data);
    figure(7);
    imagesc(edges_x, edges_y, bin_values_2dgau);
    rectangle('Position', [img_width/2-bgWidth/2, img_height/2-bgWidth/2, bgWidth, bgWidth], ...
              'EdgeColor', '#E3170D', 'LineWidth', 0.7, 'Curvature', [1, 1]);
    rectangle('Position', [img_width/2-5, img_height/2-5, 10, 10], ...
              'EdgeColor', '#E3170D', 'LineWidth', 0.7, 'Curvature', [1, 1]);
    axis image;
    title('2d Gau')
    colorbar; % 显示颜色条
end

%     c == [ECC, Ori]
%     
%     % 假设dat.fix.xpos和dat.fix.ypos是1x51的double类型数据
%     % 设定图像的尺寸和分辨率
%     img_width = expt.winRect(3);  % 二维空间的宽度
%     img_height = expt.winRect(4); % 二维空间的高度
%     
%     % 假设dat.fix.xpos和dat.fix.ypos是1x51的double类型数据
%     xpos = dat.fix.xpos;
%     ypos = dat.fix.ypos;
%     % 1. 对二维平面中的点进行核密度估计（KDE），以估算每个位置的点出现概率
%     % 创建一个网格来表示二维空间
%     edges_x = linspace(1, img_width, 50);  % 网格的X边界
%     edges_y = linspace(1, img_height, 50); % 网格的Y边界
%     
%     % 2. 使用二维直方图，计算数据在网格上的分布
%     counts = hist3([xpos', ypos'], 'Edges', {edges_x, edges_y});
%     
%     % 3. 将结果转换为二维密度
%     density_map = counts' / sum(counts(:));  % 将计数转换为密度
%     
%     % 4. 绘制热图
%     figure;
%     imagesc(edges_x, edges_y, density_map);  % 使用imagesc绘制热图
%     axis xy; % 使y轴方向与图像坐标系一致
%     colorbar; % 显示颜色条
%     title('Density Heatmap');
%     xlabel('X Position');
%     ylabel('Y Position');
% 
% end
% 
% 
% 
% 
% 
% 
% 
% sessionFileName = sprintf('%s.mat',files(1).subj);
% expt = load(fullfile(dirs.mat,sessionFileName),'expt').expt;
% geom = load(fullfile(dirs.mat,sessionFileName),'geometry').geometry.displayArea;
% if nfiles~=height(res_Table)
%     error('some trials are losted')
% end
% 
% for p=1:nfiles
%     % load fix data
%     dat  = load(fullfile(dirs.fix,[files(p).fname '.mat'])); dat = dat.dat;
%     if isempty(dat.time)
%         warning('no data for %s, empty file',files(p).fname);
%         continue;
%     end
%     
%     % get msgs
%     msgs    = loadMsgs(fullfile(dirs.msgsO,[files(p).fname '.txt']));
%     [times,what,msgs] = parseMsgs_Form(msgs);
% 
% %     if ~strcmp(lastRead,sessionFileName)
% %         lastRead = sessionFileName;
% %         fInfo = [sess.expt.stim.fInfo];
% %     end
% %     qWhich= strcmp({fInfo.name},what{1});
%     
%     % load img, if only one
% %     if ~~exist(fullfile(dirs.AOImasks,what{1}),'file')
% %         img.data = imread(fullfile(dirs.AOImasks,what{1}));
% %     elseif ~~exist(fullfile(sess.expt.stim(qWhich).fInfo.folder,what{1}),'file')
% %         img.data = imread(fullfile(sess.expt.stim(qWhich).fInfo.folder,what{1}));
% %     else
% %     grayImg = uint8(sess.expt.stim*255);
%     ut = UT(geom.width/10, expt.winRect(3), res_Table.headDist(p));
%     bgCenter = expt.winRect(3:4)/2;
%     tgCenter_ = [res_Table.ECC(p),res_Table.Orient(p)];
%     tgCenter = ut.Pol2Rect(tgCenter_);
%     
%     grayImg = genStim(expt.winRect, ut, res_Table.bgContrast(p), ...
%         1, tgCenter, expt.GaborSF, expt.GaborWidth, ...
%         expt.GaborOrient, expt.bgWidth, res_Table.seed(p));
%     img.data     = uint8 (cat(3,grayImg,grayImg,grayImg)*255);
% 
% 
% 
% %     if ~isempty(img)
% %         % get position on screen
% %         stimRect = [0,0,size(img.data,2),size(img.data,1)];
% %         img.x    = linspace(stimRect(1),stimRect(3),size(img.data,2));
% %         img.y    = linspace(stimRect(2),stimRect(4),size(img.data,1));
% %     end
% %     
% %     % plot
% %     if ~ishghandle(fhndl)
% %         fhndl = figure('Units','normalized','Position',[0 0 1 1]);  % make fullscreen figure
% %     else
% %         figure(fhndl);
% %         clf;
% %     end
% %     set(fhndl,'Visible','on');  % assert visibility to bring window to front again after keypress
% %     drawFix(dat,dat.fix,[dat.I2MCopt.xres dat.I2MCopt.yres],img,[dat.I2MCopt.missingx dat.I2MCopt.missingy],sprintf('subj %s, trial %03d, stim: %s',files(p).subj,files(p).runnr,what{1}));
% %     pause
% %     if ~ishghandle(fhndl)
% %         return;
% %     end
% end
% if ishghandle(fhndl)
%     close(fhndl);
% end
% 
% rmpath(genpath(dirs.funclib));                  % cleanup path

