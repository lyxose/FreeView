
% 2025-09-21 by LuYuxin
% 调试时勿重复运行以免丢失变量
clear variables; clear global; clear mex; fclose('all'); clc
rootDir = fileparts(mfilename('fullpath'));
addpath(fullfile(rootDir,'/Methods/function_library_cus/ANA'))    

%%
dbstop if error % for debugging: trigger a debug point when an error occurs
% ========================= 板块开关设置 =========================
% 通过设置下方 logical 变量控制各功能板块是否运行（数据准备及序列分析只需运行一次，即可改为false）
RUN_DATA_PREP          = true;    % 数据准备与筛选（只需运行一次，即可改为false）
PLOT_HEATMAP           = true;
PLOT_BASIC_STATS       = true;    % 基本分布可视化（fixTable/trial分布）
PLOT_SECTORS           = false;    % 注视点在原空间（不同颜色扇区）上的分布图
PLOT_PIE               = false;    % 选定时窗饼图
PLOT_SCAN_DIAGRAM      = false;    % 注视点扇区扫描示意图
PLOT_ANG_SCAN          = false;    % 角度扫描曲线
PLOT_ANG_PROP_SCAN     = false;    % 角度占比扫描曲线
PLOT_TIME_SERIES       = true;    % 时间/注视序列分析与绘图
PLOT_SECTOR_CORR       = false;    % 扇区相关性/Fisher z分析与热图
PLOT_16BIN             = false;    % 柱状图统计（16-bin/Axis-Gap/Card-Obli）
PLOT_AXIS_EFFECT       = false;    % 轴主效应柱状图
PLOT_OBLI_EFFECT       = false;    % 斜主效应柱状图
PLOT_TRIAL_SLIDING     = false;    % trial-level方窗滑动窗口分析与绘图
PLOT_TRIAL_SLIDING_GAU = true;    % trial-level高斯滑动窗口分析与绘图
PLOT_TRIALWIZE         = false;    % trialwise比例动态
PLOT_SPECTRUM          = false;    % 频谱分析（FFT）
SAVE_FIGURES           = false;    % 保存所有图为PNG
SAVE_FIXATION_VIDEO    = false;    % 生成注视点视频（较慢，默认关闭）

BAR_BY_COUNT           = true;    % 使用fix计数（而非时间曲线求均值）统计柱状图效应
EFFECT_BY_PROPORTION   = true;    % 使用所占比例（而非zscore）绘制柱状图效应

% PLOT_SUBJ_TRIAL_CURVES = false;    % 逐被试trial-level曲线
% ========================= 共用参数 =========================
heat_binSize = 25; % for heatmap, in pixel

CardColor = [251,  4, 255]/255;   % [245, 229, 38]/255;           %   0               90              180             270
ObliColor = [160, 95, 255]/255;  % [174, 198, 55]/255; %           45              135             225             315    
GapColor  = [57, 198, 255]/255;   % [68, 51, 205]/255;            %       22.5    67.5    112.5   157.5   202.5   247.5   292.5   337.5
AxisColor = mix_RGB_by_HSL(CardColor, ObliColor); % 混合色，作为 Axis 颜色
n_bin_FT = 16;
cmap16_FT = repmat([CardColor; GapColor; ObliColor; GapColor], ceil(n_bin_FT/4), 1);
cmap16_FT = cmap16_FT(1:n_bin_FT,:);

skip_corr = false;  % ignore correct trials
keep_nFix = 11;     % ignore 极端数据
keep_Time = 4000;   % ignore 极端数据
angbinSize = 11.25; % for angle curve
% default R_min R_max
R_min = 1;
R_max = 7.5;
% specific ring band
% R_min = 1;
% R_max = 3;
% R_min = 3;
% R_max = 5;
% R_min = 5;
% R_max = 7;



statX = 'Time';  % 横轴可选: 'fixation', 'time'
% statY = 'zScore'; % 纵轴可选: 'zScore', 'Proportion'
% normMode_FT = 'zScore'; % 'zscore','minmax','sum1','demean','mean1' 
doSmooth = false; % 是否对 time-count 和 time-density 结果进行高斯平滑

shift_FT = 360 / n_bin_FT / 2;                       % shift 使 0° 落在 bin 中心
edges_FT = linspace(0,360,n_bin_FT+1) - shift_FT;    % edges，供 histcounts 使用


% ========================= 主脚本 =========================
% ver = 'v1.5'; % 
vers = {'v1','v1.5','v2'};
figHandles = findall(0, 'Type', 'figure');
if ~isempty(figHandles)
    disp('即将关闭所有图窗，按空格键继续...');
    pause;
    close all;
end
for verc = vers
ver = strrep(verc{1}, '.', '_');
disp(['######## Processing version: ', verc{1},' ########']);
% 版本特异性参数设置
verDir = sprintf('E:\\Desktop\\临时文件\\AttenSamp\\FreeView_%s', verc{1}); % 数据主目录
matDir = '\Data';   % Dat_Subxx_Sesxx.mat 所属文件夹到主目录的相对路径
getThrData = true; % 是否获取v1.5的阈值数据（仅在RUN_DATA_PREP时有效）
if getThrData
    learn_stage_n = 0;  % trials in learning stage, to make the results comparable to FreeView_v1.5
else
    learn_stage_n = 72;  % trials in learning stage, to make the results comparable to FreeView_v1.5
end
last_trial =  480;  % total trials in formal experiment
switch lower(verc{1})
    case 'v1'
        exclude_sub = [22]; % 1号被试在~4200ms后没有数据; 22号被试只完成了一半不到（225trial）; 15,17,21,27号被试提前结束
    case 'v1.5'
        exclude_sub = [8,15]; % 8 & 15 眼动控制精度问题，数据不可靠；4号被试整体fixation number过少，没有大于10的fixnumber；16和17有两次阈限阶段，可能影响trial level结果。
        if getThrData
            thrDir = '\Analysis\Threshold\Processed_data\fixDet'; % 阈值数据目录，相对于 verDir
            learn_stage_n = 0;      % learning stage is in threshold stage
            last_trial =  480;  % 截取480-72，从而与v1和v2匹配
        else
            learn_stage_n = 0;      % learning stage is in formal experiment
            last_trial =  480-72;  % total trials in formal experiment
        end
    case 'v2'
        exclude_sub = [1,2,3,10];% 1,2,3在调试阶段，程序不同；10号被试只完成了357 trials (键盘故障)，少于-3sd
        matDir = '\Data\Formal'; 
    otherwise
        error('ver must be "v1" or "v1.5" or "v2".');
end

% 设置路径
[dirs, resultDir] = setup_paths(verDir, matDir, rootDir, verc{1});
if strcmpi(ver, 'v1_5') && getThrData
    dirs.thrDir = fullfile(verDir, thrDir); % 只要设置了路径，就会读取并合并阈限数据到正式数据之前
end

timeRes_FT = 10; % ms 
win_left = 1000; win_right = 4000; % 时间窗(ms)，闭区间
xlab = 'Time (ms)';
if doSmooth; smooth_sigma_ms = 100; end  % 高斯平滑窗口（ms）


% ---- 数据准备与筛选 ----
if RUN_DATA_PREP
    [resfiles, sub_ses_res, select_sess] = get_eye_data_files(dirs.fix);
    pairs_FT.(ver) = sub_ses_res(select_sess, :);
    pairs_FT.(ver) = pairs_FT.(ver)(~ismember(pairs_FT.(ver)(:,1), exclude_sub), :); % exclude subjects
    % 提取和合并数据
    [fixTable, Nsubj, img_width, img_height, ut, center, digPlace.(ver)] = build_fixTable(select_sess, exclude_sub, sub_ses_res, dirs, resfiles, learn_stage_n, last_trial, skip_corr); 

    % 数据筛选与变量准备
    [fixTable, start_FT, dur_FT, angles_FT, xpos_FT, ypos_FT, sub_FT, ses_FT, tri_FT, dnfix_FT,subj_stats] = filter_fixTable_for_analysis(fixTable, R_max);
    % 外部筛选（可变）：
    % 筛选时间窗口：
    win_select_Fixs = (win_left<=start_FT & start_FT<=win_right) | (win_left<=start_FT+dur_FT & start_FT+dur_FT<=win_right);
    % 再筛选特定的半径区间：
    win_select_Fixs = win_select_Fixs & (fixTable.r(~fixTable.dropFix)>=R_min & fixTable.r(~fixTable.dropFix) <= R_max);
    % 保存表格到csv
    writetable(fixTable, sprintf('ALL_fixTable_%dSubj.csv',Nsubj));
    CleanedTable.(ver) = fixTable(~fixTable.dropFix,:);
end

% ---- 1-4s热图 ----
if PLOT_HEATMAP
    plot_fixTable_heatmap(xpos_FT(win_select_Fixs), ypos_FT(win_select_Fixs), 25, 1920, 1080);
    set(gcf, 'Name', [verc{1}, ' (', map_labels(ver), ')', '--全时窗热图'], 'NumberTitle', 'off');
    title('Fixation Position Heatmap');
end


% 注视点在原空间（不同颜色扇区）上的分布图
if PLOT_SECTORS
    plot_fixation_sectors(xpos_FT(win_select_Fixs), ypos_FT(win_select_Fixs), ...
                         img_width, img_height, R_max, repmat([AxisColor;GapColor],n_bin_FT/2,1), edges_FT, ut);
    set(gcf, 'Name', [verc{1}, ' (', map_labels(ver), ')', '--注视点扇区分布(AG)'], 'NumberTitle', 'off');
    plot_fixation_sectors(xpos_FT(win_select_Fixs), ypos_FT(win_select_Fixs), ...
                         img_width, img_height, R_max, cmap16_FT, edges_FT, ut);
    set(gcf, 'Name', [verc{1}, ' (', map_labels(ver), ')', '--注视点扇区分布(COG)'], 'NumberTitle', 'off');
end

% ---- 角度扫描曲线 ----
if PLOT_ANG_SCAN && RUN_DATA_PREP % 角度扫描不支持跳过预处理！！
    startAngle = -22.5/2;
    angles_FT_ = angles_FT(win_select_Fixs);
    sub_FT_   = sub_FT(win_select_Fixs);
    ses_FT_   = ses_FT(win_select_Fixs);
    % 角度滑动统计未做时间加权（与柱状图/时程曲线不等价）！且边界处理是宽松的，即只要注视点的开始或结束时间在时间窗内即被考虑
    if PLOT_SCAN_DIAGRAM
    plot_fixation_scan_diagram(xpos_FT(win_select_Fixs), ypos_FT(win_select_Fixs), img_width, img_height, cmap16_FT, edges_FT, angbinSize, startAngle,360,0);
    set(gcf, 'Name', [verc{1}, ' (', map_labels(ver), ')', '--全角度扫描示意图'], 'NumberTitle', 'off');
    end
    [centers360_FT, m360_FT, se360_FT, subCounts360n_FT] = analyze_angle_curve(angles_FT_, sub_FT_, ses_FT_, pairs_FT.(ver), Nsubj, angbinSize, 360, 'zScore',startAngle);
    plot_angle_curve(subCounts360n_FT, centers360_FT, 360, angbinSize, CardColor, GapColor, ObliColor, 'zScore',[-1, 1]);
    set(gcf, 'Name', [verc{1}, ' (', map_labels(ver), ')', '--全角度扫描'], 'NumberTitle', 'off');

    % 45° fold: 将360°的曲线分为8段，每段宽度为45°，对每段内的数据取平均
    % 计算45°fold：将centers360_FT mod 45后相同的点分组，行内平均
    binSize = angbinSize;
    foldPeriod = 45;
    if PLOT_SCAN_DIAGRAM
    plot_fixation_scan_diagram(xpos_FT(win_select_Fixs), ypos_FT(win_select_Fixs), img_width, img_height, cmap16_FT, edges_FT, binSize, startAngle,foldPeriod,0);
    set(gcf, 'Name', [verc{1}, ' (', map_labels(ver), ')', '--45°fold扫描示意图'], 'NumberTitle', 'off');
    end
    mod_angles = mod(centers360_FT - startAngle, foldPeriod);
    [uniq_mod, ~, ic] = unique(round(mod_angles, 8)); % 防止浮点误差
    nbins45 = numel(uniq_mod);
    subCounts45n_FT = zeros(Nsubj, nbins45);
    for k = 1:nbins45
        idx = (ic == k);
        subCounts45n_FT(:,k) = mean(subCounts360n_FT(:,idx), 2, 'omitnan');
    end
    centers45 = uniq_mod + startAngle;
    plot_angle_curve(subCounts45n_FT, centers45, 45, binSize, CardColor, GapColor, ObliColor, 'zScore',[-.4, .8]);
    set(gcf, 'Name', [verc{1}, ' (', map_labels(ver), ')', '--45°fold扫描'], 'NumberTitle', 'off');

    % 90° fold: 将360°的曲线分为4段，每段宽度为90°，对每段内的数据取平均
    foldPeriod = 90;
    binSize = angbinSize;
    if PLOT_SCAN_DIAGRAM
    plot_fixation_scan_diagram(xpos_FT(win_select_Fixs), ypos_FT(win_select_Fixs), img_width, img_height, cmap16_FT, edges_FT, binSize, startAngle,90);
    set(gcf, 'Name', [verc{1}, ' (', map_labels(ver), ')', '--90°fold扫描示意图'], 'NumberTitle', 'off');
    end
    mod_angles90 = mod(centers360_FT - startAngle, foldPeriod);
    [uniq_mod90, ~, ic90] = unique(round(mod_angles90, 8)); % 防止浮点误差
    nbins90 = numel(uniq_mod90);
    subCounts90n_FT = zeros(Nsubj, nbins90);
    for k = 1:nbins90
        idx = (ic90 == k);
        subCounts90n_FT(:,k) = mean(subCounts360n_FT(:,idx), 2, 'omitnan');
    end
    centers90 = uniq_mod90 + startAngle;
    plot_angle_curve(subCounts90n_FT, centers90, 90, binSize, CardColor, GapColor, ObliColor, 'zScore',[-.4, .8]);
    set(gcf, 'Name', [verc{1}, ' (', map_labels(ver), ')', '--90°fold扫描'], 'NumberTitle', 'off');
end


% ---- 角度占比扫描曲线 ----
if PLOT_ANG_PROP_SCAN && RUN_DATA_PREP % 角度扫描不支持跳过预处理！！
    angles_FT_ = angles_FT(win_select_Fixs);
    sub_FT_   = sub_FT(win_select_Fixs);
    ses_FT_   = ses_FT(win_select_Fixs);
    startAngle = -22.5/2;
    if PLOT_SCAN_DIAGRAM
    plot_fixation_scan_diagram(xpos_FT(win_select_Fixs), ypos_FT(win_select_Fixs), img_width, img_height, cmap16_FT, edges_FT, angbinSize, startAngle,360,0);
    set(gcf, 'Name', [verc{1}, ' (', map_labels(ver), ')', '--全角度扫描示意图'], 'NumberTitle', 'off');
    end
    % 角度滑动统计未做时间加权（与柱状图/时程曲线不等价）！且边界处理是宽松的，即只要注视点的开始或结束时间在时间窗内即被考虑
    [centers360_FT, m360_FT, se360_FT, subCounts360n_FT] = analyze_angle_curve(angles_FT_, sub_FT_, ses_FT_, pairs_FT.(ver), Nsubj, angbinSize, 360, 'sum1', startAngle);
    plot_angle_curve(subCounts360n_FT, centers360_FT, 360, angbinSize, CardColor, GapColor, ObliColor, 'Proportion',[0.02, 0.04]);
    set(gcf, 'Name', [verc{1}, ' (', map_labels(ver), ')', '--全角度占比扫描'], 'NumberTitle', 'off');
    % 45° fold: 将360°的曲线分为8段，每段宽度为45°，直接求比例和
    binSize = angbinSize;
    foldPeriod = 45;
    if PLOT_SCAN_DIAGRAM
    plot_fixation_scan_diagram(xpos_FT(win_select_Fixs), ypos_FT(win_select_Fixs), img_width, img_height, cmap16_FT, edges_FT, binSize, startAngle,foldPeriod,0);
    set(gcf, 'Name', [verc{1}, ' (', map_labels(ver), ')', '--45°fold扫描示意图'], 'NumberTitle', 'off');
    end
    mod_angles = mod(centers360_FT - startAngle, foldPeriod);
    [uniq_mod, ~, ic] = unique(round(mod_angles, 8)); % 防止浮点误差
    nbins45 = numel(uniq_mod);
    subCounts45n_FT = zeros(Nsubj, nbins45);
    for k = 1:nbins45
        idx = (ic == k);
        subCounts45n_FT(:,k) = sum(subCounts360n_FT(:,idx), 2, 'omitnan');
    end
    centers45 = uniq_mod + startAngle;
    plot_angle_curve(subCounts45n_FT, centers45, 45, binSize, CardColor, GapColor, ObliColor, 'Proportion',[0.22, 0.30], false);
    set(gcf, 'Name', [verc{1}, ' (', map_labels(ver), ')', '--45°fold占比扫描'], 'NumberTitle', 'off');
    % ---------- Fit N°-period cosine with fixed constant = chance ----------
    % 针对v1_5版本，调整统计角度和r刻度标签位置
    if strcmpi(ver, 'v1_5')
        statAngle = 22.5;
        rtickTheta = 5;
        rtickHAlign = 'left';
        rtickVAlign = 'top';
    else
        statAngle = 0;
        rtickTheta = 15;
        rtickHAlign = 'right';
        rtickVAlign = 'top';
    end
    [fit45.(ver), fit45_stats.(ver)] = fit_periodic_cosine_phase(subCounts45n_FT, centers45, 45, 0.25, ...
        'plot', true, 'histColor', [0.5 0.5 0.5], 'scatterColor', [0 0 0], 'ExpColor', AxisColor, ...
        'statAngle', statAngle, 'rtickTheta', rtickTheta, 'rtickHAlign', rtickHAlign, 'rtickVAlign', rtickVAlign);  
    set(gcf, 'Name', [verc{1}, ' (', map_labels(ver), ')', '--45°fold相位vtest'], 'NumberTitle', 'off');
    AngelCurve45.(ver).data = subCounts45n_FT;
    AngelCurve45.(ver).centers = centers45;
    % 90° fold: 将360°的曲线分为4段，每段宽度为90°，直接求比例和
    foldPeriod = 90;
    binSize = angbinSize;
    if PLOT_SCAN_DIAGRAM
    plot_fixation_scan_diagram(xpos_FT(win_select_Fixs), ypos_FT(win_select_Fixs), img_width, img_height, cmap16_FT, edges_FT, binSize, startAngle,foldPeriod,0);
    set(gcf, 'Name', [verc{1}, ' (', map_labels(ver), ')', '--90°fold扫描示意图'], 'NumberTitle', 'off');
    end
    mod_angles90 = mod(centers360_FT - startAngle, foldPeriod);
    [uniq_mod90, ~, ic90] = unique(round(mod_angles90, 8)); % 防止浮点误差
    nbins90 = numel(uniq_mod90);
    subCounts90n_FT = zeros(Nsubj, nbins90);
    for k = 1:nbins90
        idx = (ic90 == k);
        subCounts90n_FT(:,k) = sum(subCounts360n_FT(:,idx), 2, 'omitnan');
    end
    centers90 = uniq_mod90 + startAngle;
    plot_angle_curve(subCounts90n_FT, centers90, 90, binSize, CardColor, GapColor, ObliColor, 'Proportion',[0.105, 0.162], false);
    set(gcf, 'Name', [verc{1}, ' (', map_labels(ver), ')', '--90°fold占比扫描'], 'NumberTitle', 'off');
    AngelCurve90.(ver).data = subCounts90n_FT;
    AngelCurve90.(ver).centers = centers90;
    % 180° fold: 将360°的曲线分为2段，每段宽度为180°，直接求比例和
    % foldPeriod = 180;
    % binSize = angbinSize;
    % if PLOT_SCAN_DIAGRAM
    % plot_fixation_scan_diagram(xpos_FT(win_select_Fixs), ypos_FT(win_select_Fixs), img_width, img_height, cmap16_FT, edges_FT, binSize, startAngle,foldPeriod,0);
    % set(gcf, 'Name', [verc{1}, ' (', map_labels(ver), ')', '--180°fold扫描示意图'], 'NumberTitle', 'off');
    % end
    % mod_angles180 = mod(centers360_FT - startAngle, foldPeriod);
    % [uniq_mod180, ~, ic180] = unique(round(mod_angles180, 8)); % 防止浮点误差
    % nbins180 = numel(uniq_mod180);
    % subCounts180n_FT = zeros(Nsubj, nbins180);
    % for k = 1:nbins180
    %     idx = (ic180 == k);
    %     subCounts180n_FT(:,k) = sum(subCounts360n_FT(:,idx), 2, 'omitnan');
    % end
    % centers180 = uniq_mod180 + startAngle;
    % plot_angle_curve(subCounts180n_FT, centers180, 180, binSize, CardColor, GapColor, ObliColor, 'Proportion',[0.05, 0.085], true);
    % set(gcf, 'Name', [verc{1}, ' (', map_labels(ver), ')', '--180°fold占比扫描'], 'NumberTitle', 'off');
end

% ---- 柱状图统计 ----
if PLOT_16BIN || PLOT_AXIS_EFFECT || PLOT_OBLI_EFFECT
    % bar plot统一口径
    % if BAR_BY_COUNT
        % 统计每个被试在选定时间窗内的注视点角度分布（16-bin计数）
        bin_counts = zeros(Nsubj, n_bin_FT);
        for si = 1:Nsubj
            bin_counts(si,:) = histcounts(mod(angles_FT(win_select_Fixs & sub_FT == pairs_FT.(ver)(si,1)) + shift_FT, 360) - shift_FT, edges_FT);
        end
        
        % 16bin柱状图数据
        mean_sectors.(ver) = normalize_by_dim(bin_counts - 0.5 * (bin_counts(:,[end,1:end-1]) + bin_counts(:,[2:end,1])), 'zScore'); % 环形邻居去趋势 then normalize
        mean_sectors_prop.(ver) = normalize_by_dim(bin_counts, 'sum1'); % 转为比例
        % Axis/Gaps/Card/Obli效应
        Axis_Effect.(ver) = sum(bin_counts(:,1:2:16),2) ./ sum(bin_counts,2);   % 占比数据
        Obli_Effect.(ver) = [sum(bin_counts(:,1:4:16),2), sum(bin_counts(:,3:4:16),2)]./ sum(bin_counts(:,1:2:16),2); % 占比数据./ sum(bin_counts,2); %
        Obli_Effect_inall.(ver) = [sum(bin_counts(:,1:4:16),2), sum(bin_counts(:,3:4:16),2)]./ sum(bin_counts,2); %        
        Gap_Mean.(ver) = sum(bin_counts(:,2:2:16),2) ./ sum(bin_counts,2);
end

% ---- 饼图绘制 ----
if PLOT_PIE
    % axis vs gap
    all_effects = [mean(Obli_Effect_inall.(ver)(:,1)); mean(Obli_Effect_inall.(ver)(:,2)); mean(Gap_Mean.(ver))];
    plot_pie(all_effects, [CardColor; ObliColor; GapColor], 'ax', [], 'show_text', true, 'legend', {'Card.','Obli.','Gap'});
    set(gcf, 'Name', [verc{1}, ' (', map_labels(ver), ')', '--选定时窗饼图'], 'NumberTitle', 'off');
    % axis vs gap
    all_effects = [mean(Axis_Effect.(ver)); mean(Gap_Mean.(ver))];
    plot_pie(all_effects, [AxisColor; GapColor], 'ax', [], 'show_text', true,  'legend', {'Axis','Gap'});
    set(gcf, 'Name', [verc{1}, ' (', map_labels(ver), ')', '--选定时窗饼图（Axis）'], 'NumberTitle', 'off');
    % card vs obli
    all_effects = [mean(Obli_Effect.(ver)(:,1)); mean(Obli_Effect.(ver)(:,2))];
    plot_pie(all_effects, [CardColor; ObliColor], 'ax', [], 'show_text', true, 'legend', {'Card.','Obli.'});
    set(gcf, 'Name', [verc{1}, ' (', map_labels(ver), ')', '--选定时窗饼图（Card vs Obli）'], 'NumberTitle', 'off');
end

% ---- 16bin分别统计 ----
if PLOT_16BIN
    labels_16bin = arrayfun(@(x) sprintf('%.0f°',x), (edges_FT(1:end-1)+edges_FT(2:end))/2, 'uni', 0);
    plot_bar_multi(mean_sectors.(ver), cmap16_FT, labels_16bin, ...
        'ylabel', 'Detrend Z Score', ...
        'xlabel', '', ...
        'title', sprintf('Mean Sector Density (%d–%d ms)', win_left, win_right), ...
        'xtickMode', 'deg', ...
        'showInd', true, ...
        'showIndNum', false);
    set(gcf, 'Name', [verc{1}, ' (', map_labels(ver), ')', '--16扇区分别统计柱状图'], 'NumberTitle', 'off');

    plot_bar_multi(mean_sectors_prop.(ver), cmap16_FT, labels_16bin, ...
        'ylabel', 'Proportion', ...
        'xlabel', '', ...
        'title', sprintf('Mean Sector Density (%d–%d ms)', win_left, win_right), ...
        'xtickMode', 'deg', ...
        'showInd', true, ...
        'showIndNum', false, ...
        'baseline',  1/16, ...
        'baselineLabel', 'Chance');
    set(gcf, 'Name', [verc{1}, ' (', map_labels(ver), ')', '--16扇区分别统计柱状图（比例）'], 'NumberTitle', 'off');
end

% ---- 轴主效应统计 ----
if PLOT_AXIS_EFFECT
    plot_violin_prop(Axis_Effect.(ver), AxisColor, {'Axis'}, ...
        'ylabel', 'Proportion', ...
        'xlabel', '', ...
        'title', sprintf('Axis Effect (%d–%d ms)', win_left, win_right), ...
        'showInd', true, ...
        'showIndNum', false, ...
        'showIndLink', false, ...
        'chanceLevel', 0.5);
    set(gcf, 'Name', [verc{1}, ' (', map_labels(ver), ')', '--AG总效应'], 'NumberTitle', 'off');
end

% ---- 斜主效应统计 ----
if PLOT_OBLI_EFFECT
    plot_violin_prop(Obli_Effect.(ver)(:,1), cmap16_FT(1,:), {'Card'}, ...
        'ylabel', 'Proportion', ...
        'xlabel', '', ...
        'title', sprintf('Cardinal Effect (%d–%d ms)', win_left, win_right), ...
        'showInd', true, ...
        'showIndNum', false, ...
        'showIndLink', false, ...
        'chanceLevel', 0.5, ...
        'chanceLabel', 'Chance');
    set(gcf, 'Name', [verc{1}, ' (', map_labels(ver), ')', '--CO总效应'], 'NumberTitle', 'off');

    plot_violin_prop(Obli_Effect_inall.(ver), cmap16_FT([1,3],:), {'Card','Obli'}, ...
        'ylabel', 'Proportion', ...
        'xlabel', '', ...
        'title', sprintf('Cardinal Effect (%d–%d ms)', win_left, win_right), ...
        'showInd', true, ...
        'showIndNum', false, ...
        'showIndLink', true, ...
        'chanceLevel', mean(Gap_Mean.(ver))/2, ...
        'chanceLabel', 'Gap/2');
    set(gcf, 'Name', [verc{1}, ' (', map_labels(ver), ')', '--CO总效应SI'], 'NumberTitle', 'off');
end



% ---- 时序数据绘图 ----
if PLOT_TIME_SERIES
    % 时序数据计算
    [time_series_Mat.(ver), timeCenters_FT] = process_TimeCourse_count( ...
        angles_FT, start_FT, dur_FT, sub_FT, ses_FT, pairs_FT.(ver), Nsubj, ...
        timeRes_FT, keep_Time, edges_FT, shift_FT, n_bin_FT);
    if doSmooth
        sigma_frames = smooth_sigma_ms / timeRes_FT;
        time_series_Mat.(ver) = gaussian_smooth_along_dim(time_series_Mat.(ver), sigma_frames, 2);
    end    
    xWin = timeCenters_FT;
    % this_tseri_Mat = normalize_by_dim(time_series_Mat.(ver), 'zScore');
    % this_tseri_Mat = this_tseri_Mat - 0.5 * (this_tseri_Mat(:,:,[end,1:end-1])+this_tseri_Mat(:,:,[2:end,1]));
    % ylab = 'Detrend Z Score';
    % this_tseri_Mat = normalize_by_dim(time_series_Mat.(ver), 'sum1');
    ylab = 'Proportion';

    % 时序16扇区热图
    % plot_sector_time_heatmap(timeCenters_FT, n_bin_FT, this_tseri_Mat, edges_FT, ylab, [0.04,0.085]);
    % set(gcf, 'Name', [verc{1}, ' (', map_labels(ver), ')', '--16扇区热图（比例）'], 'NumberTitle', 'off');
    % xline(win_left, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.2);

    % plot_sector_stats_time_series(this_tseri_Mat, timeCenters_FT, Nsubj);
    % xline(win_left, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.2);
    % set(gcf, 'Name', [verc{1}, ' (', map_labels(ver), ')', '--16扇区变异性时程'], 'NumberTitle', 'off');
    
    this_tseri_Mat = normalize_by_dim(time_series_Mat.(ver), 'sum1'); % 不能用normalize后的zscore数据
    % this_tseri_Mat = this_tseri_Mat ./ sum(this_tseri_Mat,3);

    seriesAxis.(ver) = squeeze(sum(this_tseri_Mat(:,:,1:2:16),3)); % nsubj x nfix/time
    seriesGap.(ver)  = squeeze(sum(this_tseri_Mat(:,:,2:2:16),3));
        
    % Axis 占比时程   
    cfg.doStats = true;
    cfg.ylabel  = 'Proportion';
    cfg.xlabel  = xlab;
    cfg.chanceLabel = 'Chance';
    cfg.chanceLevel = 0.5;
    cfg.statTail= 'both';
    plot_single_prop(seriesAxis.(ver), xWin, AxisColor, 'Axis', cfg);
    xline(win_left, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.2, 'HandleVisibility', 'off');
    set(gcf, 'Name', [verc{1}, ' (', map_labels(ver), ')', '--AG时程'], 'NumberTitle', 'off');
    title('Time Course of Axis Effect');
    
    % Card 时程--单线占比版本
    seriesCard.(ver) = squeeze(sum(this_tseri_Mat(:,:,1:4:16),3)./seriesAxis.(ver));
    seriesObli.(ver) = squeeze(sum(this_tseri_Mat(:,:,3:4:16),3)./seriesAxis.(ver));
    cfg.statTail= 'right';
    % plot_comparison(seriesCard.(ver), seriesObli.(ver), xWin, cmap16_FT([1,3],:), {'Card','Oblique'}, cfg);
    plot_single_prop(seriesCard.(ver), xWin, CardColor, 'Cardinal', cfg);
    xline(win_left, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.2);
    set(gcf, 'Name', [verc{1}, ' (', map_labels(ver), ')', '--CO时程'], 'NumberTitle', 'off');
    title('Time Course of Cardinal Effect');

    % Card vs Oblique 时程--双线对比版本
    seriesCard.(ver) = squeeze(sum(this_tseri_Mat(:,:,1:4:16),3));  % 原始比例可与gap对比
    seriesObli.(ver) = squeeze(sum(this_tseri_Mat(:,:,3:4:16),3));  % 
    cfg.statTail= 'both';
    % if strcmpi(ver, 'v2')
    %     cfg.chanceLabel = 'Gap/2';
    %     cfg.chanceLevel = mean(Gap_Mean.(ver))/2;
    % else
    %     cfg.chanceLabel = 'Chance';
    %     cfg.chanceLevel = 0.25;
    % end
    cfg.chanceLabel = 'Gap/2';
    cfg.chanceLevel = seriesGap.(ver)/2;
    % cfg.chanceLabel = 'Chance';
    % cfg.chanceLevel = 0.5/2;
    plot_comparison(seriesCard.(ver), seriesObli.(ver), xWin, cmap16_FT([1,3],:), {'Card','Oblique'}, cfg);
    xline(win_left, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.2, 'HandleVisibility', 'off');
    set(gcf, 'Name', [verc{1}, ' (', map_labels(ver), ')', '--CO时程SI'], 'NumberTitle', 'off');
    title('Time Course of Cardinal Effect');
end


 %%

% ---- trial-level滑动窗口 ----
total_trial = last_trial - learn_stage_n;
cfg_win = struct();
cfg_win.edges_FT     = edges_FT;
cfg_win.shift_FT     = shift_FT;
cfg_win.cmap16_FT    = cmap16_FT;
cfg_win.win_left     = win_left;
cfg_win.win_right    = win_right;
cfg_win.step_trials  = 1;
cfg_win.total_trial  = total_trial;
cfg_win.digPlace     = digPlace.(ver);
cfg_win.mode         = 'continuous';
cfg_win.doStats      = true;
cfg_win.learn_stage_n= learn_stage_n;
cfg_win.normMode_FT  = 'zScore';
cfg_win.xlabel       = 'Trial number';
cfg_win.ylabel       = 'Proportion';
cfg_win.doSmooth     = false;
cfg_win.ver          = ver;
cfg_win.AxisColor    = AxisColor;
cfg_win.Card_in_Axis = true; % 计算Card占Axis的比例
cfg_win.cutoff       = 300;  % 计算学习趋势时只选取前300trial
% v1部分被试提前结束，需要跳过，生成tailMask: nSubj x total_trial，超出每个被试maxNtrial的trial位置标记为true
tailMask = false(height(subj_stats), total_trial);
for i = 1:height(subj_stats)
    if ~isnan(subj_stats.maxNtrial(i))
        tailMask(i, subj_stats.maxNtrial(i)+1:end) = true;
    end
end
cfg_win.tailMask = tailMask;

GauWin = 40;
BoxWin = 120;
%  trial-level矩形滑动窗口分析与绘图 
if PLOT_TRIAL_SLIDING
    [TseriesAxis1, TseriesGap1, TseriesCard1, TseriesObli1, xWin1] = compute_sliding_window_series(CleanedTable.(ver), cfg_win, total_trial);
    plot_sliding_window_analysis(TseriesAxis1, TseriesGap1, TseriesCard1, TseriesObli1, xWin1, cfg_win);
    % plot_sliding_window_analysis_corr(TseriesAxis1, TseriesGap1, TseriesCard1, TseriesObli1, xWin1, cfg_win);
    set(gcf, 'Name', [verc{1}, ' (', map_labels(ver), ')', '--trial方窗滑窗'], 'NumberTitle', 'off');
    % set(gcf, 'Name', [verc{1}, ' (', map_labels(ver), ')', '--trial方窗滑窗_带corr'], 'NumberTitle', 'off');

    cfg2 = cfg_win;
    cfg2.mode         = 'block';
    cfg2.blockSize    = 96;
    cfg2.win_trials   = BoxWin;
    if cfg2.blockSize>cfg_win.win_trials
        [TseriesAxis2, TseriesGap2, TseriesCard2, TseriesObli2, xWin2] = compute_sliding_window_series(CleanedTable.(ver), cfg2, total_trial);
        plot_sliding_window_analysis(TseriesAxis2, TseriesGap2, TseriesCard2, TseriesObli2, xWin2, cfg2);
        % plot_sliding_window_analysis_corr(TseriesAxis2, TseriesGap2, TseriesCard2, TseriesObli2, xWin2, cfg_win);
        set(gcf, 'Name', [verc{1}, ' (', map_labels(ver), ')', '--trial方窗滑窗分Block'], 'NumberTitle', 'off');
        % set(gcf, 'Name', [verc{1}, ' (', map_labels(ver), ')', '--trial方窗滑窗分Block_带corr'], 'NumberTitle', 'off');
    end
end
%  trial-level高斯滑动窗口分析与绘图 
if PLOT_TRIAL_SLIDING_GAU %&& 0
    cfg_gau = cfg_win;
    cfg_gau.win_trials = GauWin; % 高斯核标准差（trials），可调整

    [TseriesAxisG, TseriesGapG, TseriesCardG, TseriesObliG, xWinG] = compute_gaussian_window_series(CleanedTable.(ver), cfg_gau, total_trial);
    plot_sliding_window_analysis(TseriesAxisG, TseriesGapG, TseriesCardG, TseriesObliG, xWinG, cfg_gau);
    set(gcf, 'Name', [verc{1}, ' (', map_labels(ver), ')', '--trial高斯滑窗'], 'NumberTitle', 'off');
    cfg_gau.Card_in_Axis = false; % 计算Card和Obli的原始比例
    [TseriesAxisG, TseriesGapG, TseriesCardG, TseriesObliG, xWinG] = compute_gaussian_window_series(CleanedTable.(ver), cfg_gau, total_trial);    
    plot_sliding_window_analysis_obli(TseriesObliG, xWinG, cfg_gau);
    set(gcf, 'Name', [verc{1}, ' (', map_labels(ver), ')', '--trial高斯滑窗（Obli）_SI'], 'NumberTitle', 'off');
end
if PLOT_TRIALWIZE
    [TseriesAxis, TseriesGap, TseriesCard, TseriesObli, xWin_trialwise] = compute_trialwise_series(CleanedTable.(ver), cfg_win, total_trial);
    plot_sliding_window_analysis(TseriesAxis, TseriesGap, TseriesCard, TseriesObli, xWin_trialwise, cfg_win);
    set(gcf, 'Name', [verc{1}, ' (', map_labels(ver), ')', '--trialwise比例动态'], 'NumberTitle', 'off');
end
%%
% cfg_win.Card_in_Axis = false; % 计算Card和Obli的原始比例
% %  trial-level矩形滑动窗口分析与绘图, with GAP
% if PLOT_TRIAL_SLIDING
%     [TseriesAxis1, TseriesGap1, TseriesCard1, TseriesObli1, xWin1] = compute_sliding_window_series(CleanedTable.(ver), cfg_win, total_trial);
%     plot_sliding_window_analysis(TseriesAxis1, TseriesGap1, TseriesCard1, TseriesObli1, xWin1, cfg_win);
%     set(gcf, 'Name', [verc{1}, ' (', map_labels(ver), ')', '--trial方窗滑窗SI'], 'NumberTitle', 'off');

%     cfg2 = cfg_win;
%     cfg2.mode         = 'block';
%     cfg2.blockSize    = 96;
%     cfg2.win_trials   = BoxWin;
%     if cfg2.blockSize>cfg_win.win_trials
%         [TseriesAxis2, TseriesGap2, TseriesCard2, TseriesObli2, xWin2] = compute_sliding_window_series(CleanedTable.(ver), cfg2, total_trial);
%         plot_sliding_window_analysis(TseriesAxis2, TseriesGap2, TseriesCard2, TseriesObli2, xWin2, cfg2);
%         set(gcf, 'Name', [verc{1}, ' (', map_labels(ver), ')', '--trial滑窗分Block_SI'], 'NumberTitle', 'off');
%     end
% end
 
% %  trial-level高斯滑动窗口分析与绘图 
% if PLOT_TRIAL_SLIDING_GAU %&& 0
%     cfg_gau = cfg_win;
%     cfg_gau.win_trials = GauWin; % 高斯核标准差（trials），可调整
%     [TseriesAxisG, TseriesGapG, TseriesCardG, TseriesObliG, xWinG] = compute_gaussian_window_series(CleanedTable.(ver), cfg_gau, total_trial);
%     plot_sliding_window_analysis(TseriesAxisG, TseriesGapG, TseriesCardG, TseriesObliG, xWinG, cfg_gau);
%     % plot_sliding_window_analysis_corr(TseriesAxisG, TseriesGapG, TseriesCardG, TseriesObliG, xWinG, cfg_win);
%     set(gcf, 'Name', [verc{1}, ' (', map_labels(ver), ')', '--trial高斯滑窗SI'], 'NumberTitle', 'off');
%     % set(gcf, 'Name', [verc{1}, ' (', map_labels(ver), ')', '--trial高斯滑窗_corr'], 'NumberTitle', 'off');
% end

% ---- 逐被试trial-level曲线 ----
% if PLOT_SUBJ_TRIAL_CURVES 
%     plot_subject_1curves(TseriesAxis1 - TseriesGap1, xWin1, CardColor, 'Axis - Gap');
%     ylabel(ylab);
%     set(gcf, 'Name', [verc{1}, ' (', map_labels(ver), ')', '--AG逐被试trial时程'], 'NumberTitle', 'off');

%     plot_subject_2curves(TseriesCard1, TseriesObli1, xWin1, cmap16_FT([1,3],:), 'Card', 'Obli');
%     ylabel(ylab);
%     set(gcf, 'Name', [verc{1}, ' (', map_labels(ver), ')', '--CO逐被试trial时程'], 'NumberTitle', 'off');
% end

% ---- 基本分布可视化 ---- 
if PLOT_BASIC_STATS
    plot_fixTable_distributions(CleanedTable.(ver), img_width, img_height, heat_binSize, ut, center, verc{1});
    set(gcf, 'Name', [verc{1}, ' (', map_labels(ver), ')', '--fixation统计'], 'NumberTitle', 'off');

    plot_trial_distributions(CleanedTable.(ver));
    set(gcf, 'Name', [verc{1}, ' (', map_labels(ver), ')', '--trial统计'], 'NumberTitle', 'off');

    mask_win = (start_FT >= win_left) & (start_FT <= win_right);
    plot_fixTable_heatmap(xpos_FT(mask_win), ypos_FT(mask_win), 25, img_width, img_height);
    set(gcf, 'Name', [verc{1}, ' (', map_labels(ver), ')', sprintf('--%d-%dms热图', win_left, win_right)], 'NumberTitle', 'off');
    title(sprintf('Fixation Position Heatmap (%d-%dms)', win_left, win_right));

end


% ---- 扇区相关性/Fisher z分析 ----
if PLOT_SECTOR_CORR
    [stat_axis, stat_card, time_ms, fz_axis_obs, fz_card_obs, perm_z8, perm_z4] = analyze_sector_fisherz(permute(time_series_Mat.(ver), [2, 3, 1]), timeRes_FT);
    plot_sector_fisherz_clusters(stat_axis, stat_card, time_ms, fz_axis_obs, fz_card_obs, perm_z8, perm_z4, CardColor, GapColor, ObliColor);
    xline(win_left, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.2);
    set(gcf, 'Name', [verc{1}, ' (', map_labels(ver), ')', '--16扇区Fisher z时程'], 'NumberTitle', 'off');
end

% ---- 频谱分析（FFT） ----
if PLOT_SPECTRUM         
    spec_result = analyze_subject_spectrum(TseriesAxis1);
    plot_subject_spectrum(spec_result, CardColor);
    set(gcf, 'Name', [verc{1}, ' (', map_labels(ver), ')', '--AG trial频谱'], 'NumberTitle', 'off');
end

% ---- 保存所有图到PNG ----
if SAVE_FIGURES
    figHandles = findall(0, 'Type', 'figure');
    for iFig = 1:numel(figHandles)
        fig = figHandles(iFig);
        figName = get(fig, 'Name');
        if isempty(figName)
            figName = sprintf('Figure_%d', fig.Number);
        end
        figName = regexprep(figName, '[\/:*?"<>|]', '_');
        saveas(fig, fullfile(resultDir, [figName, '.png']));
        % exportgraphics(fig, fullfile(resultDir, [figName, '_hq.png']), 'Resolution', 300);
    end
end

% ---- 生成注视点视频（较慢，默认关闭） ----
if SAVE_FIXATION_VIDEO
    win_left = 0; % ms
    win_right= keep_Time; % ms
    frame_rate = 50; % fps
    frame_time_res = 1000/frame_rate; % ms
    n_frames = round((win_right - win_left) / frame_time_res);
    video_width = 1920; % pixels
    video_height= 1080; % pixels
    make_fixation_sector_video(start_FT, dur_FT, xpos_FT, ypos_FT, win_left, win_right, frame_rate, video_width, video_height, R_max, cmap16_FT, edges_FT, ut);
end

if SAVE_FIGURES
    close all; 
end

end

%% 交互效应
cd (rootDir); cd Results;
if numel(vers) >= 2
    compare_groups = strrep(vers, 'v1.5', 'v1_5');
    % Axis Effect
    [effect_all, group_all, cmp_labels] = CrossV_collect_data(Axis_Effect, compare_groups);
    cmp_labels = map_labels(cmp_labels);
    CrossV_plot_violin_anova(effect_all, group_all, cmp_labels, 'Axis Proportion');
    print(gcf, 'Axis_Effect.png', '-dpng', '-r300');
    % Obli Effect (Card proportion in Axis)
    [effect_cmp, group_cmp, cmp_labels] = CrossV_collect_data(Obli_Effect, compare_groups);
    cmp_labels = map_labels(cmp_labels);
    CrossV_plot_violin_anova(effect_cmp, group_cmp, cmp_labels, 'Card. Proportion', 0.5);
    print(gcf, 'Obli_Effect.png', '-dpng', '-r300');
end

%% 角度扫描对比
cd (rootDir); cd Results;


%% ACC angle scan
% ---- ACC angle scan: 按角度滑动统计 trial-level accuracy 并绘图 ----

[centersACC, meanACC, seACC, subjACC] = analyze_angle_accuracy(CleanedTable.(ver), pairs_FT.(ver), Nsubj, angbinSize, 360, startAngle);
% 绘图（复用已有绘图接口，与 angle proportion 保持一致）
plot_angle_curve(subjACC, centersACC, 360, angbinSize, CardColor, GapColor, ObliColor, 'Proportion', [0, 1]);
set(gcf, 'Name', [verc{1}, ' (', map_labels(ver), ')', '--角度-Accuracy扫描'], 'NumberTitle', 'off');


%% 交互辅助函数

% function [effect_all, group_all, ver_names] = CrossV_collect_data(EffectStruct, compare_groups, do_diff)
    % 从EffectStruct中提取指定组的效应数据，准备进行交互效应分析

% function CrossV_plot_bar_anova(data, group, labels, ylabel_str, chanceLevel)
    % 绘制分组柱状图（占比数据），以chance level为基线，显著性与chance level比较
%% 角度扫描对比辅助函数
% helper: 找到包含峰值的连通区并返回宽度
% 比较 v1 与 v1_5 的 45°-fold 相位曲线（将 v1_5 平移 22.5° 后叠加，并计算 FWHM）
if exist('AngelCurve45','var') && isfield(AngelCurve45,'v1') && isfield(AngelCurve45,'v1_5')
    foldP = 45;
    shift_ang = 22.5;

    % extract
    c1 = AngelCurve45.v1.centers(:)';
    data1 = AngelCurve45.v1.data; % nsub x nbin
    c2 = AngelCurve45.v1_5.centers(:)';
    data2 = AngelCurve45.v1_5.data; % nsub x nbin

    % center to [-foldP/2, foldP/2)
    cc1 = mod(c1,foldP); cc1(cc1>=foldP/2) = cc1(cc1>=foldP/2) - foldP;
    cc2 = mod(c2 + shift_ang,foldP); cc2(cc2>=foldP/2) = cc2(cc2>=foldP/2) - foldP;

    % sort centers and reorder data accordingly
    [cc1s, idx1] = sort(cc1); data1s = data1(:, idx1);
    [cc2s, idx2] = sort(cc2); data2s = data2(:, idx2);
    % --- 高斯平滑参数 ---
    sigma = 0; % 单位: step
    smooth1 = data1s;
    smooth2 = data2s;
    if sigma>0
        for si = 1:size(data1s,1)
            smooth1(si,:) = gaussian_smooth_along_dim(data1s(si,:), sigma, 2, 'circular');
        end
        for si = 1:size(data2s,1)
            smooth2(si,:) = gaussian_smooth_along_dim(data2s(si,:), sigma, 2, 'circular');
        end
    end

    % compute per-subject FWHM using smoothed curves
    n1 = size(smooth1,1); n2 = size(smooth2,1);
    fwhm1 = nan(n1,1); fwhm2 = nan(n2,1);
    for si = 1:n1
        fwhm1(si) = compute_subject_FWHM(cc1s, smooth1(si,:), 0, foldP);
    end
    for si = 1:n2
        fwhm2(si) = compute_subject_FWHM(cc2s, smooth2(si,:), 0, foldP);
    end

    % remove NaN
    fwhm1 = fwhm1(~isnan(fwhm1));
    fwhm2 = fwhm2(~isnan(fwhm2));

    % unpaired statistics
    if numel(fwhm1) >= 2 && numel(fwhm2) >= 2
        [h_t,p_t,ci,stats_t] = ttest2(fwhm1, fwhm2);
        p_rank = ranksum(fwhm1, fwhm2);
    else
        h_t = NaN; p_t = NaN; stats_t = []; p_rank = NaN;
    end

    % print summary
    fprintf('Subjects (v1): %d, (v1_5 shifted): %d\n', numel(fwhm1), numel(fwhm2));
    fprintf('Mean FWHM (v1, smoothed) = %.3f deg, (v1_5 shifted, smoothed) = %.3f deg\n', nanmean(fwhm1), nanmean(fwhm2));
    if ~isnan(p_t)
        fprintf('Unpaired t-test: t(%d)=%.3f p=%.4f\n', stats_t.df, stats_t.tstat, p_t);
        fprintf('Wilcoxon ranksum p = %.4f\n', p_rank);
    else
        warning('Not enough subjects for statistics.');
    end

    % --- 颜色设置（与交互效应一致） ---
    color_v1   = [0.2 0.6 1];    % 蓝色
    color_v1_5 = [1 0.5 0.2];    % 橙色

    % plot mean curves with SE (smoothed)
    m1 = mean(smooth1,1,'omitnan'); se1 = std(smooth1,0,1,'omitnan')./sqrt(max(1,size(smooth1,1)));
    m2 = mean(smooth2,1,'omitnan'); se2 = std(smooth2,0,1,'omitnan')./sqrt(max(1,size(smooth2,1)));

    fig = figure('Name','AngelCurve45_compare_smoothed','Color','w');
    hold on;
    % patch SE
    patch([cc1s, fliplr(cc1s)], [m1+se1, fliplr(m1-se1)], color_v1, 'FaceAlpha', 0.25, 'EdgeColor','none');
    patch([cc2s, fliplr(cc2s)], [m2+se2, fliplr(m2-se2)], color_v1_5, 'FaceAlpha', 0.18, 'EdgeColor','none');
    h1 = plot(cc1s, m1, '-','Color', color_v1,  'LineWidth', 2);
    h2 = plot(cc2s, m2, '-','Color', color_v1_5, 'LineWidth', 2);
    xlim([-foldP/2, foldP/2]);
    xlabel(sprintf('Angle (centered, -%.1f to +%.1f°)', foldP/2, foldP/2));
    ylabel('Proportion (mean, smoothed)');
    legend([h1,h2], {'v1','v1\_5 (shifted 22.5°)'}, 'Location','best');
    box off; hold off;
else
    warning('AngelCurve45.v1 or AngelCurve45.v1_5 missing — cannot compute per-subject FWHM.');
end

% ---------------- helper: compute subject FWHM (no interpolation)
function w = compute_subject_FWHM(centers, y, peakAngle, period)
    % centers: 1 x nbins (sorted, in [-period/2, period/2))
    % y: 1 x nbins
    % peakAngle: angle to take as peak reference (e.g., 0)
    % baseline uses values at +/- period/2 (i.e., +/-22.5 for 45deg fold)
    w = NaN;
    if isempty(centers) || isempty(y)
        return;
    end
    % find nearest index to peakAngle
    [~, idxPeak] = min(abs(centers - peakAngle));
    peakVal = y(idxPeak);

    % baseline: mean of bins nearest to +period/2 and -period/2
    target = period/2;
    [~, idx_b_pos] = min(abs(centers - target));
    [~, idx_b_neg] = min(abs(centers + target));
    baselineVal = mean(y(unique([idx_b_pos, idx_b_neg])),'omitnan');

    half = baselineVal + (peakVal - baselineVal)/2;

    mask = y >= half;
    if ~any(mask)
        return;
    end
    % find runs
    d = diff([0, mask(:)', 0]);
    starts = find(d==1);
    ends = find(d==-1)-1;

    % handle wrap-around (mask at both ends)
    if mask(1) && mask(end)
        % merge last and first runs
        % last run starts at starts(end), first run ends at ends(1)
        starts = [starts(1) + (starts(end)-starts(end)), starts(2:end)]; %#ok<NASGU>
        % simpler: treat combined run indices explicitly
        % combined start = starts(end), combined end = ends(1)
        combinedStart = starts(end);
        combinedEnd = ends(1);
        % create arrays with combined run as first element
        starts2 = [combinedStart, starts(2:end-0)];
        ends2 = [combinedEnd, ends(2:end)];
        starts = starts2; ends = ends2;
    end

    % find the segment containing peak
    segIdx = find(starts <= idxPeak & ends >= idxPeak, 1);
    if isempty(segIdx)
        % try to find nearest segment by distance
        dstart = min(abs(starts - idxPeak));
        dend = min(abs(ends - idxPeak));
        if ~isempty(dstart) && ~isempty(dend)
            [~, segIdx] = min([dstart, dend]);
            segIdx = segIdx(1);
        else
            return;
        end
    end

    sIdx = starts(segIdx);
    eIdx = ends(segIdx);
    xstart = centers(sIdx);
    xend = centers(eIdx);
    w = xend - xstart;
    if w < 0
        w = w + period;
    end
end

%% 正确率角度扫描

function [centers, groupMean, groupSE, subjCurves] = analyze_angle_accuracy(visTable, pairs, nsbj, binSize, foldPeriod, startAngle)
% 对每个被试 (按 pairs 列表) 使用滑动角度窗统计 trial-level accuracy（比例）
% visTable: table of trials/fixations（必须包含 trial-level target angle 与 correct 标志）
% pairs: nsbj x 2 matrix (sub, ses)
% nsbj: number of subjects
% binSize, foldPeriod, startAngle: 同 analyze_angle_curve 用法
%
% 返回：
% centers: 角度中心向量
% groupMean, groupSE: 均值与SE（按被试先求被试内比例，再组均值）
% subjCurves: nsbj x ncenters，每被试每角度窗的 accuracy（0-1 或 NaN）

% --- 务实型列名识别 ---
varnames = lower(visTable.Properties.VariableNames);
findcol = @(pat) find(contains(varnames,pat),1);
subcol = findcol('sub');
if isempty(subcol), subcol = findcol('subject'); end
sescol = findcol('ses'); if isempty(sescol), sescol = findcol('session'); end
trialcol = findcol('trial'); 
% target angle 列（多选尝试）
targcol = findcol('target') ; 
if isempty(targcol), targcol = findcol('targ'); end
if isempty(targcol), targcol = findcol('angle'); end
% correct 列
corrcol = findcol('corr'); if isempty(corrcol), corrcol = findcol('correct'); end
if isempty(corrcol), corrcol = findcol('iscorrect'); end
if isempty(corrcol), corrcol = findcol('acc'); end

if isempty(subcol) || isempty(sescol) || isempty(trialcol) || isempty(targcol) || isempty(corrcol)
    error('analyze_angle_accuracy: cannot find required columns in visTable (need sub, ses, trial, target-angle, correct).');
end

subcol = visTable.Properties.VariableNames{subcol};
sescol = visTable.Properties.VariableNames{sescol};
trialcol = visTable.Properties.VariableNames{trialcol};
targcol = visTable.Properties.VariableNames{targcol};
corrcol = visTable.Properties.VariableNames{corrcol};

% 确定扫描分辨率（与 analyze_angle_curve 保持一致）
if foldPeriod == 360
    step = 0.01;
else
    step = 0.1;
end
centers = startAngle + (0:round(foldPeriod/step)-1)*step;
nbins = numel(centers);

subjCurves = nan(nsbj, nbins);

% 每个被试/会话抽取 trial-level target angles 与 correctness（按 trial 去重）
for si = 1:nsbj
    subjID = pairs(si,1); sesID = pairs(si,2);
    mask = visTable.(subcol) == subjID & visTable.(sescol) == sesID;
    if ~any(mask)
        continue;
    end
    tbl_sub = visTable(mask,:);
    % 以 trial 为单位取第一个记录的 target angle 与 correctness（若 visTable 已是 trial-level 则无影响）
    [uTrials, ia] = unique(tbl_sub.(trialcol),'stable');
    targAngles = tbl_sub.(targcol)(ia);
    corrVals = tbl_sub.(corrcol)(ia);
    % 将 correctness 二值化（尽量兼容 logical/0-1/true/false/string）
    if iscell(corrVals)
        corrBins = cellfun(@(x) any(strcmpi(x,{'1','true','t','correct','yes'})), corrVals);
    else
        corrBins = double(corrVals);
    end
    % 若 correctness 非 0/1，则尝试归一化到 0/1
    if ~all(ismember(unique(corrBins(~isnan(corrBins))),[0 1]))
        corrBins = double(corrBins ~= 0);
    end

    % 对每个滑动中心统计该试次的 accuracy（角度窗口内的 trial 的正确率）
    for ii = 1:nbins
        center = centers(ii);
        leftEdge = center - binSize/2;
        angs_shifted = mod(targAngles - leftEdge, foldPeriod);
        idx = angs_shifted < binSize;
        if any(idx)
            subjCurves(si,ii) = mean(corrBins(idx));
        else
            subjCurves(si,ii) = NaN;
        end
    end
end

% 组内平均与 SE
groupMean = nanmean(subjCurves,1);
validN = sum(~isnan(subjCurves),1);
groupSE = nanstd(subjCurves,0,1) ./ sqrt(max(1, validN));

end
%% 实用小函数（脚本内定义）

% 路径设置函数：根据版本自动设置数据和结果目录，并切换工作路径
function [dirs, resultDir] = setup_paths(verDir, matDir, wkRootDir, ver)
    % 设置数据目录
    dirs.mat = fullfile(verDir, matDir);
    dirs.fix = fullfile(verDir, 'Analysis', 'Processed_data', 'fixDet');
    % 结果目录（按版本分文件夹）
    resultDir = fullfile(wkRootDir, 'Results', ver);
    if ~exist(resultDir, 'dir')
        mkdir(resultDir);
    end
    cd(resultDir); % 切换工作路径到结果文件夹
end

% function out_labels = map_labels(in_labels)
  % 版本标签映射函数：将版本代码映射为更易读的标签

% function [resfiles, sub_ses_res, select_sess] = get_eye_data_files(dirs_fix)
    % 数据文件获取函数：从指定目录筛选有效的被试-会话文件列表

% function s = sig_symbol(p)
% 将 p 值转为显著性符号

% function out = ternary(cond, a, b)
% 简单三元运算


% function [stat, details] = cluster_based_permutation_test(seriesA, seriesB, varargin)
    % seriesA, seriesB: nSubj x nTime 矩阵
    % 可选参数：
    % 'n_perm' - 置换次数，默认5000
    % 'alpha'  - 显著性水平，默认0.05
    % 'tail'   - 检验尾数，'both'（默认），'left'，'right'
    % 输出：
    % stat - 结构体，包含以下字段：
    %   .prob         - 每个时间点的p值
    %   .t_obs       - 观察到的t值
    %   .cluster_inds - 观察到的显著簇的起止索引
    %   .cluster_p   - 每个簇的p值
    % details - 结构体，包含以下字段：
    %   .t_threshold  - t值阈值
    %   .perm_max     - 每次置换的最大簇质量
    %   .obs_clusters  - 观察到的簇信息（起始索引，结束索引，簇质量）
    %   .tail         - 检验尾数

% function smoothed = gaussian_smooth_along_dim(mat, sigma, dim)
    % 高斯平滑函数：对输入矩阵沿指定维度进行高斯平滑

% function [sIdx, eIdx] = find_runs(mask)
%  辅助函数：在逻辑掩码中查找连续为真的区间起止索引

%% ---- 辅助函数：绘制时程对比图 ----
% function plot_comparison(seriesA, seriesB, xWin, colors, labels, cfg)
    % 绘制两组时程数据对比，进行组间显著性比较


% function plot_single(series, xWin, color, label, cfg)
    % 绘制单一时程数据，进行与0的显著性比较

% function plot_single_prop(series, xWin, color, label, cfg)
    % 绘制占比数据（0~1），与chance level（默认0.5）比较显著性

% function plot_bar_multi(data, colors, labels, varargin)
    % 多列输入的绘制（如16-bin），与原代码一致

% function plot_violin_prop(data, colors, labels, varargin)
    % 数据为比例，展示分布特征

%%

% ---- 预处理函数定义 ----
% function [fixTable, Nsubj, img_width, img_height, ut, center, digPlace] = build_fixTable(select_sess, exclude_sub, sub_ses_res, dirs, resfiles, learn_stage_n, last_trial, skip_corr)
    % 构建注视点数据表格

% ---- 筛选函数定义 ----
% function [fixTable, start_FT, dur_FT, angles_FT, xpos_FT, ypos_FT, sub_FT, ses_FT, tri_FT, dnfix_FT] = filter_fixTable_for_analysis(fixTable, R_max)
    % 1. 剔除起始fixation在中心1°以内的fixation
    % 2. 剔除最末连续几个位于目标1°以内的fixation
    % 3. 剔除r>R_max（背景范围）的fixation
    % 4. 剔除最末几个指向目标方向（误差连续减小）的fixation
    % 输出数据表格并标记离群数据


%% --------------------- 分析函数实现 ---------------------

% function analyze_axis_obli_correlation(Axis_Effect, Obli_Effect, varargin)
    % 分析 Axis Effect 与 Obli Effect 之间的相关性

% function [sub_fix_bin_count, nfix_FT] = process_fixation_count_series( ...
%     angles_FT, dnfix_FT, sub_FT, ses_FT, pairs_FT, nsbj_FT, keep_nFix, ...
%     n_bin_FT, shift_FT, edges_FT)
% 按注视序号聚合角度数据为每被试 × 注视序号 × 16-bin 的计数矩阵

% function [sub_time_bin_count, timeCenters_FT] = process_TimeCourse_count( ...
%         angles_FT, start_FT, dur_FT, sub_FT, ses_FT, pairs_FT, nsbj_FT, ...
%         timeRes_ms, keep_Time_ms, edges_FT, shift_FT, n_bin_FT)
% process_TimeCourse_count
%   将每个注视（angles_FT, start_FT, dur_FT）按时间轴汇总为“被试 × 时间 × 16-bin”的矩阵。

% function [sub_time_bin_density, video_time_axis] = process_time_density_series( ...
%     start_FT, dur_FT, xpos_FT, ypos_FT, sub_FT, ses_FT, pairs_FT, ...
%     video_width, video_height, ut, cfg)
    % process_time_DENSITY_series 计算每个被试在时间序列上的空间密度分布


% 辅助函数：处理单个被试
% function [sector_density, dur] = single_subject_time_density(...
%     si, pairs_FT, sub_FT, ses_FT, start_FT, dur_FT, xpos_FT, ypos_FT, ...
%     n_frames, frame_time_res, video_width, video_height, ...
%     radius, gauss_template, sector_map, n_sectors, ...
%     write_video)

% 辅助函数：保存视频
% function save_subject_density_video(video_frames, pair, frame_time_res, n_frames)

% 
% function [sub_fix_bin_density, nfix_FT] = process_fix_density_series( ...
    % start_FT, dur_FT, xpos_FT, ypos_FT, dnfix_FT, sub_FT, ses_FT, pairs_FT, ...
    % video_width, video_height, ut, cfg)
% PROCESS_FIX_DENSITY_SERIES 计算每个被试按注视序号的空间密度分布

% 辅助函数：处理单个被试（按注视序号）
% function [sector_density, dur] = single_subject_fix_density(...

% function [TseriesAxis, TseriesGap, TseriesCard, TseriesObli, xWin] = compute_sliding_window_series(visTable, cfg, total_trial)
% COMPUTE_SLIDING_WINDOW_SERIES 对 visTable 执行基于 trial 的滑动窗口统计

% function [TseriesAxis, TseriesGap, TseriesCard, TseriesObli, xWin] = compute_gaussian_window_series(visTable, cfg, total_trial)
% COMPUTE_GAUSSIAN_WINDOW_SERIES 对 visTable 执行基于 trial 的高斯滑动窗口统计

% function spec_result = analyze_subject_spectrum(seriesAxis)
% 对输入的 nsbj x n_trials 时序数据做频谱分析和置换检验
% 返回结构体 spec_result，包含所有分析结果和可直接用于绘图的字段

% function [centers, groupMean, groupSE, subjCurves] = analyze_angle_curve(angles_FT, sub_FT, ses_FT, pairs_FT, nsbj_FT, binSize, foldPeriod, normMode)
    % 对每个被试(会话)的所有 angles_FT 进行沿角度轴的滑动计数

% function normCurves = normalize_by_dim(subjCurves, normMode, normDim)
%NORMALIZE_BY_dim Normalize subject curves along specified dimension.

% function [stat_axis, stat_card, time_ms, fz_axis_obs, fz_card_obs, perm_z8, perm_z4] = analyze_sector_fisherz(sector_counts_sub, frame_time_res)
    % 对 sector_counts_sub 做扇区相关性分析和簇置换检验
    % 返回 Axis/Card 的统计结构体和相关序列

% 封装为函数 fit_periodic_cosine_phase
% 用法示例：
% [fitResult, stats] = fit_periodic_cosine_phase(subCountsFold, centersFold, period, chance, ...
%     'plot', true, 'histColor', [0.5 0.5 0.5], 'scatterColor', [1 0 0], 'ExpColor', [1 0 0], ...
%     'ax', gca, 'statAngle', 0, 'rtickTheta', 120, 'rtickHAlign', 'right', 'rtickVAlign', 'top');
%
% 参数说明：
%   subCountsFold: Nsubj x nbins，被试数据（已fold）
%   centersFold:   1 x nbins，中心角度（deg）
%   period:        fold周期（deg），如45
%   chance:        常数基线（如0.25）
%   可选参数:
%       'plot'         : 是否绘图（默认true）
%       'histColor'    : 直方图颜色（默认[0.5 0.5 0.5]）
%       'scatterColor' : 散点颜色（默认[1 0 0]）
%       'ExpColor'     : 突出统计角度的径向线颜色（默认[1 0 0]，红色）
%       'ax'           : 指定绘图的axes句柄（默认新建figure）
%       'statAngle'    : 统计角度（deg，默认0，极坐标右侧）
%       'rtickTheta'   : r刻度标签的极角（deg，默认120，周期内单位）
%       'rtickHAlign'  : r刻度标签水平对齐（默认'right'）
%       'rtickVAlign'  : r刻度标签垂直对齐（默认'top'）

function [fitResult, stats] = fit_periodic_cosine_phase(subCountsFold, centersFold, period, chance, varargin)
    % Parse optional arguments
    p = inputParser;
    addParameter(p, 'plot', true, @(x)islogical(x)||isnumeric(x));
    addParameter(p, 'histColor', [0.5 0.5 0.5]);
    addParameter(p, 'scatterColor', [1 0 0]);
    addParameter(p, 'ExpColor', [1 0 0]);
    addParameter(p, 'ax', []);
    addParameter(p, 'statAngle', 0); % 新增：统计角度（deg）
    addParameter(p, 'rtickTheta', 120); % 新增：r刻度标签极角（周期内单位，默认120）
    addParameter(p, 'rtickHAlign', 'right'); % 新增：r刻度标签水平对齐
    addParameter(p, 'rtickVAlign', 'top');   % 新增：r刻度标签垂直对齐
    parse(p, varargin{:});
    doPlot = logical(p.Results.plot);
    histColor = p.Results.histColor;
    scatterColor = p.Results.scatterColor;
    ExpColor = p.Results.ExpColor;
    ax = p.Results.ax;
    statAngle = p.Results.statAngle;
    rtickTheta = p.Results.rtickTheta;
    rtickHAlign = p.Results.rtickHAlign;
    rtickVAlign = p.Results.rtickVAlign;

    [Nsubj, nbins] = size(subCountsFold);
    omega = 2*pi / period; % rad per degree
    theta_deg = centersFold(:)'; % 1 x nbins
    cos_term = cos(omega * theta_deg);
    sin_term = sin(omega * theta_deg);

    phases = nan(Nsubj,1);
    amps   = nan(Nsubj,1);
    beta12 = nan(Nsubj,2);

    for si = 1:Nsubj
        y = subCountsFold(si,:); % 1 x nbins (比例)
        valid = ~isnan(y);
        X = [cos_term(valid)', sin_term(valid)']; % nbins x 2
        b = X \ (y(valid)' - chance);
        beta1 = b(1); beta2 = b(2);
        R = hypot(beta1, beta2);              % 振幅
        phi = atan2(beta2, beta1);            % 相位（rad）
        phi = mod(phi + pi, 2*pi) - pi;       % [-pi, pi]
        phases(si) = phi;
        amps(si)   = R;
        beta12(si,:) = [beta1, beta2];
    end

    % Preferred angle within degree domain for the period:
    theta0_deg = phases / omega;
    theta0_deg_mod = mod(theta0_deg, period);

    % ---- Circular test: v-test of phi toward statAngle ----
    valid_idx = ~isnan(phases);
    angles_test = phases(valid_idx);
    n_angles = numel(angles_test);

    % 统计角度折算到cos相位空间
    phi0 = omega * mod(statAngle, period); % 统计角度对应的cos相位(rad)

    if n_angles >= 2
        if exist('circ_vtest','file') == 2
            [vtest_p, vtest_stat] = circ_vtest(angles_test, phi0);
            vtest_method = 'circ_vtest';
        elseif exist('circ_rtest','file') == 2
            [vtest_p, vtest_stat] = circ_rtest(angles_test - phi0);
            vtest_method = 'circ_rtest (fallback)';
        else
            R = abs(sum(exp(1i*(angles_test - phi0))));
            Rbar_tmp = R / n_angles;
            Z = n_angles * Rbar_tmp^2;
            p_approx = exp(-Z) * (1 + (2*Z - Z^2)/(4*n_angles) - (24*Z - 132*Z^2 + 76*Z^3 - 9*Z^4)/(288*n_angles^2));
            vtest_p = min(max(p_approx,0),1);
            vtest_stat = Z;
            vtest_method = 'rayleigh_approx (no CircStat)';
        end
        mean_dir = angle(mean(exp(1i*angles_test))); % mean phi in radians
        Rbar = abs(mean(exp(1i*angles_test)));
        mean_theta0_deg = mod((mean_dir / omega), period);
    else
        vtest_p = NaN; vtest_stat = NaN; vtest_method = 'insufficient';
        mean_dir = NaN; Rbar = NaN; mean_theta0_deg = NaN;
    end

    % 输出结构体
    fitResult.phase_rad = phases;
    fitResult.phase_deg = rad2deg(phases);
    fitResult.amp = amps;
    fitResult.beta = beta12;
    fitResult.theta0_deg = theta0_deg;
    fitResult.theta0_deg_mod = theta0_deg_mod;
    fitResult.vtest_method = vtest_method;
    fitResult.vtest_p = vtest_p;
    fitResult.vtest_stat = vtest_stat;
    fitResult.mean_phi = mean_dir;
    fitResult.mean_theta0_deg = mean_theta0_deg;
    fitResult.Rbar = Rbar;
    fitResult.valid_idx = valid_idx;

    stats.n_valid = sum(valid_idx);
    stats.mean_amp = nanmean(amps);
    stats.std_amp = nanstd(amps);
    stats.mean_phi = mean_dir;
    stats.mean_theta0_deg = mean_theta0_deg;
    stats.Rbar = Rbar;
    stats.vtest_p = vtest_p;
    stats.vtest_stat = vtest_stat;
    stats.vtest_method = vtest_method;

    % 打印
    fprintf('--- %d°-period cosine fit (C = %.3f, statAngle = %.1f°) ---\n', period, chance, statAngle);
    fprintf('有效被试数: %d / %d\n', sum(valid_idx), Nsubj);
    fprintf('平均振幅 (mean R): %.4f, std: %.4f\n', nanmean(amps), nanstd(amps));
    fprintf('phi mean (rad)=%.4f => preferred theta0 mean (deg in 0-%d)=%.2f deg, resultant length Rbar = %.4f\n', mean_dir, period, mean_theta0_deg, Rbar);
    fprintf('环形检验方法: %s, 统计量: %.4f, p = %.4g\n', vtest_method, vtest_stat, vtest_p);

    % ---- Visualization: polar histogram (folded to 0-period) + scatter points ----
    if doPlot
        % 允许传入ax，否则新建figure
        if isempty(ax) || ~isvalid(ax)
            f = figure('Color','none','Position',[120,80,270,240]);
            ax = polaraxes('Parent',f);
        else
            axes(ax); %#ok<LAXES>
        end
        theta_plot = theta0_deg_mod / period * 2*pi; % 0..2pi
        nBins = 16;
        binWidth = 2 * pi / nBins;
        binEdges = linspace(0, 2*pi, nBins+1) - binWidth/2;
        h1 = polarhistogram(ax, theta_plot(~isnan(theta_plot)), binEdges, ...
            'Normalization', 'probability', 'FaceAlpha', 0.5, 'FaceColor', histColor, 'EdgeColor', histColor);
        hold(ax, 'on');
        r_plot = amps;
        rmin = 1e-16;
        r_plot_log = log10(max(r_plot, rmin));
        r_plot_log_scaled = (r_plot_log + 3) / 10;
        h2 = polarplot(ax, theta_plot, r_plot_log_scaled, '.', 'MarkerSize', 12, 'Color', scatterColor);

        % 突出标记统计角度对应的径向线
        % 统计角度statAngle（deg）对应极坐标theta
        theta_stat = mod(statAngle, period) / period * 2*pi;
        rmax = 0.32;
        polarplot(ax, [theta_stat theta_stat], [0 rmax], '--', 'LineWidth', 2, 'Color', ExpColor);

        pax = ax;
        pax.RLim = [0 0.3];
        log_ticks = log10([0.001 0.01 0.1 1]);
        log_ticks_scaled = (log_ticks +3) /10;
        log_vals= {'0.001(0)', '0.01(0.1)', '0.1(0.2)', '1(0.3)'};
        pax.RTick = log_ticks_scaled;
        pax.RTickLabel = {};
        % 变换rtickTheta为极坐标角度（deg），周期内单位转为极坐标角度
        theta = deg2rad(rtickTheta * (360/period));
        for i = 1:length(log_vals)
            r = log_ticks_scaled(i);
            text(ax, theta, r, log_vals{i}, 'HorizontalAlignment', rtickHAlign, 'VerticalAlignment', rtickVAlign, 'FontSize', 10);
        end
        title(ax, sprintf('Phases folded to 0-%d°, n=%d', period, sum(valid_idx)));
        % 设置theta ticks为加粗
        thetaticks = 0:40:320;
        set(ax, 'ThetaTick', thetaticks);
        set(ax, 'ThetaTickLabel', arrayfun(@(x) sprintf('%d°', x), (thetaticks)/(360/period), 'UniformOutput', false));
        ax.ThetaAxis.FontWeight = 'bold';
        ax.Color = 'none'; % 背景透明
        if isempty(ax.Parent) || ~strcmp(get(ax.Parent,'Type'),'figure')
            set(gcf, 'Color', 'none');
        end
        hold(ax, 'off');

        % 标记p值在指定极角最外圈外侧一点点
        r_p = pax.RLim(2) * 1.03;
        p_str = sprintf('V test\np = %.3f', vtest_p);
        if vtest_p < 0.05
            p_col = [1 0 0];
        else
            p_col = [0 0 0];
        end
        text(ax, 59/360*2*pi, r_p, p_str, 'Color', p_col, 'FontWeight', 'bold', 'FontSize', 8, ...
            'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');
    end
end


%% ----------------------- 可视化函数实现 ---------------------
% function plot_fixTable_distributions(visTable, img_width, img_height, heat_binSize, ut, center, ver)
% plot_fixTable_distributions 可视化注视点表格的主要分布特征

% function plot_trial_distributions(visTable)
% plot_trial_distributions 绘制每trial的fixation数量分布（左）与trial时间分布（右）

% function plot_fixTable_heatmap(xpos_FT, ypos_FT, heat_binSize, img_width, img_height)
% plot_fixTable_heatmap 绘制基于注视点坐标的二维热图（密度图），y轴翻转使0在上

% function plot_sector_time_heatmap(timeCenters_FT, nSector, sub_time_bin_FT, edges_FT, ylab)
% 绘制按时间的16扇区组均值热图

% function plot_sector_stats_time_series(sub_time_bin_FT, timeCenters_FT, nsbj)
% 绘制扇区组统计（range/std）随时间的均值±SE曲线

% function plot_sliding_window_analysis_corr(TseriesAxis, TseriesGap, TseriesCard, TseriesObli, xWin, cfg)
% 上图：Axis-Gap（v1_5时为Gap-Axis，需翻转），下图：Card-Obli

% function plot_sliding_window_analysis(TseriesAxis, TseriesGap, TseriesCard, TseriesObli, xWin, cfg)
% PLOT_SLIDING_WINDOW_ANALYSIS 绘制滑动窗口分析结果

% function plot_subject_spectrum(spec_result, CardColor)
% 绘制频谱分析结果，支持自定义主色 CardColor

% function plot_sector_fisherz_clusters(stat_axis, stat_card, time_ms, fz_axis_obs, fz_card_obs, perm_z8, perm_z4, CardColor, GapColor, ObliColor)
% PLOT_SECTOR_FISHERZ_CLUSTERS 绘制扇区相关性 Fisher z 时程及显著性簇

% ---- 扇区滑动绘图函数 ----
% function plot_angle_curve(centers, groupMean, groupSE, foldPeriod, binSize, nsbj_FT, CardColor, GapColor, ObliColor, normMode)
    % 绘制按角度滑动计数曲线，带组均值和SE阴影

% function make_fixation_sector_video(start_FT, dur_FT, xpos_FT, ypos_FT, win_left, win_right, frame_rate, video_width, video_height, R_max, cmap16_FT, sector_edges, ut)
% 生成注视点视频，背景为16扇区，每帧显示当前注视点及扇区计数
function plot_fixation_scan_diagram(xpos, ypos, img_width, img_height, cmap16, sector_edges, angbinSize, startAngle, foldPeriod, centerAngle)
% 绘制注视点扫描示意图，支持多个等距小扇形（foldPeriod控制）
% xpos, ypos: 注视点坐标
% img_width, img_height: 视频尺寸
% cmap16: 16色扇区配色
% sector_edges: 16个扇区边界
% angbinSize: 小扇形圆心角（度）
% startAngle: 圆弧起始角度（度）
% centerAngle: 扇形中心位置（度）
% foldPeriod: 小扇形周期（度），如360/4=90时，绘制4个扇形

if nargin < 10 || isempty(centerAngle)
    centerAngle = 0;
end 
if nargin < 9 || isempty(foldPeriod)
    foldPeriod = 360;
end
if nargin < 8 || isempty(startAngle)
    startAngle = 0;
end
if nargin < 7 || isempty(angbinSize)
    angbinSize = 11.25;
end

center_pix = [img_width, img_height] / 2;
distFromCenter = sqrt((xpos - center_pix(1)).^2 + (ypos - center_pix(2)).^2);
Rcircle = max(distFromCenter) * 1.01;

fig = figure('Name', 'Fixation Scan Diagram', 'Color', 'w');
ax = axes('Parent', fig);
hold(ax, 'on');
set(ax, 'YDir', 'reverse', 'DataAspectRatio', [1 1 1], ...
    'XTick', [], 'YTick', [], 'Box', 'off', 'Visible', 'off');

% 计算注视点相对于中心的角度 
angles = atan2d(-(ypos - center_pix(2)), xpos - center_pix(1)); % y轴向下，负号
angles = mod(angles, 360);

% 获取所有小扇形的中心角度
n = round(360 / foldPeriod);
sector_centers = mod(centerAngle + (0:n-1)*foldPeriod, 360);

% 判断哪些点落在任意一个小扇形内
in_sector = false(size(angles));
for k = 1:numel(sector_centers)
    theta0 = sector_centers(k);
    theta1 = theta0 - angbinSize/2;
    theta2 = theta0 + angbinSize/2;
    theta1_mod = mod(theta1, 360);
    theta2_mod = mod(theta2, 360);
    if theta1_mod <= theta2_mod
        in_sector = in_sector | (angles >= theta1_mod & angles <= theta2_mod);
    else
        in_sector = in_sector | (angles >= theta1_mod | angles <= theta2_mod);
    end
end

% 绘制落在灰色小扇形内的点（深灰色）
scatter(ax, xpos(in_sector), ypos(in_sector), 1, [0.3 0.3 0.3], 'filled', ...
    'MarkerFaceAlpha', 0.8);

% 绘制其他fixation点（浅灰色）
scatter(ax, xpos(~in_sector), ypos(~in_sector), 1, [0.7 0.7 0.7], 'filled', ...
    'MarkerFaceAlpha', 0.5);

% 绘制35%浅灰色实线圆圈（表示背景范围）
theta_circ = linspace(0, 2*pi, 360);
plot(ax, center_pix(1) + Rcircle*cos(theta_circ), ...
         center_pix(2) - Rcircle*sin(theta_circ), ...
         'Color', [0.8 0.8 0.8], 'LineWidth', 1.5);

% 绘制所有灰色半透明小扇形
for k = 1:numel(sector_centers)
    theta0 = sector_centers(k);
    theta_patch = linspace(theta0 - angbinSize/2, theta0 + angbinSize/2, 60);
    xp = center_pix(1) + Rcircle * cosd(theta_patch);
    yp = center_pix(2) - Rcircle * sind(theta_patch);
    patch(ax, [center_pix(1), xp, center_pix(1)], [center_pix(2), yp, center_pix(2)], ...
        [0.5 0.5 0.5], 'FaceAlpha', 0.4, 'EdgeColor', 'none');
end

% 绘制8条扇区中心虚线（直径，透明度0.8，颜色由cmap16控制）
n_sector = 16;
for si = 1:n_sector/2
    midA = mod((sector_edges(si) + sector_edges(si+1))/2, 360);
    x1 = center_pix(1) + Rcircle * cosd(midA);
    y1 = center_pix(2) - Rcircle * sind(midA);
    x2 = center_pix(1) + Rcircle * cosd(midA+180);
    y2 = center_pix(2) - Rcircle * sind(midA+180);
    color = cmap16(si,:);
    plot(ax, [x1 x2], [y1 y2], '--', 'Color', [color 0.8], 'LineWidth', 1.2);
end
% 绘制16个角度数字（整数显示整数，非整数保留一位小数）
labelPos = zeros(n_sector, 2);
labR = Rcircle * 1.12;
for si = 1:n_sector
    a1 = sector_edges(si);
    a2 = sector_edges(si+1);
    midA = mod((a1 + a2)/2, 360);
    labelPos(si,:) = center_pix + [labR * cosd(midA), -labR * sind(midA)];
    if abs(midA - round(midA)) > 1e-6 % 非整数
        lbl = sprintf('%.1f°', midA);
    else
        lbl = sprintf('%.0f°', midA);
    end
    text(labelPos(si,1), labelPos(si,2), lbl, ...
        'Parent', ax, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
        'FontWeight', 'bold', 'FontSize', 9, 'Color', [0.06 0.06 0.06]);
end

% 绘制带箭头的同心圆弧（只在第一个扇形中心）
arc_span = 30; % 圆弧跨度（度）
arc_theta = linspace(startAngle, startAngle+arc_span, 60);
arcR = Rcircle * 1.23;
arc_x = center_pix(1) + arcR * cosd(arc_theta);
arc_y = center_pix(2) - arcR * sind(arc_theta);
plot(ax, arc_x, arc_y, 'Color', [0.2 0.2 0.2 0.85], 'LineWidth', 2);
% 箭头
arrow_len = 18; % 单位: 像素
arrow_angle = startAngle + arc_span;
arr_base = [center_pix(1) + arcR * cosd(arrow_angle), center_pix(2) - arcR * sind(arrow_angle)];
theta_rad = deg2rad(arrow_angle);
tangent_dir = [-sin(theta_rad), -cos(theta_rad)];
wing_angle = 25;
Rmat = [cosd(wing_angle), -sind(wing_angle); sind(wing_angle), cosd(wing_angle)];
wing1_dir = (Rmat * (-tangent_dir)')';
wing2_dir = (Rmat' * (-tangent_dir)')';
arr1 = arr_base + arrow_len * wing1_dir / norm(wing1_dir);
arr2 = arr_base + arrow_len * wing2_dir / norm(wing2_dir);
plot(ax, [arr1(1) arr_base(1) arr2(1)], [arr1(2) arr_base(2) arr2(2)], ...
    'Color', [0.2 0.2 0.2 0.85], 'LineWidth', 2);

% 设置绘图范围
plot_radius = Rcircle * 1.25;
xlim(ax, center_pix(1) + [-plot_radius, plot_radius]);
ylim(ax, center_pix(2) + [-plot_radius, plot_radius]);
hold(ax, 'off');

end
