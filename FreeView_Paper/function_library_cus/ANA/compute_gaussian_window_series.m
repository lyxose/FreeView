function [TseriesAxis, TseriesGap, TseriesCard, TseriesObli, xWin] = compute_gaussian_window_series(visTable, cfg, total_trial)
% COMPUTE_GAUSSIAN_WINDOW_SERIES 对 visTable 执行基于 trial 的高斯滑动窗口统计
%
% 输入：
%   visTable - fixation table
%   cfg      - 配置结构体，需包含:
%       .edges_FT, .shift_FT, .win_left, .win_right, .digPlace
%       .win_trials   : 高斯核的标准差（单位：trials）
%       .learn_stage_n  : 跳过的学习阶段 trial 数
%   total_trial - 总 trial 数
%
% 输出：
%   Tseries... : nsbj × total_trial 矩阵，高斯平滑后的效应时序
%   xWin       : 1 × total_trial 窗口中心位置向量 (1:total_trial)

    % ---- 1. 参数验证与准备 ----
    assert(istable(visTable) && ~isempty(visTable), 'visTable 必须是非空表格');
    requiredFields = {'edges_FT','shift_FT','win_left','win_right','digPlace','win_trials'};
    for i = 1:numel(requiredFields)
        assert(isfield(cfg, requiredFields{i}), ['cfg 缺少字段: ', requiredFields{i}]);
    end
    
    if ~isfield(cfg, 'learn_stage_n'), cfg.learn_stage_n = 0; end
    
    trialDigits = 10^(cfg.digPlace+1);
    sesTrial = mod(visTable.TriID, trialDigits);
    sesTrial(sesTrial==0) = trialDigits;

    timeOverlap = (visTable.startT < cfg.win_right) & ((visTable.startT + visTable.dur) > cfg.win_left);
    
    pairs_SW = unique([visTable.subID, visTable.sessID], 'rows', 'stable');
    nsbj_SW  = size(pairs_SW, 1);
    
    skip_n = max(0, cfg.learn_stage_n);
    sesTrial_adj = sesTrial - skip_n;
    valid_mask = sesTrial_adj > 0 & sesTrial_adj <= total_trial;

    % ---- 2. 逐 Trial 计算瞬时 16-bin 计数 ----
    per_trial_counts = zeros(nsbj_SW, total_trial, 16);

    for ti = 1:total_trial
        % 筛选当前 trial 的有效注视点
        maskTrial = sesTrial_adj == ti & timeOverlap & valid_mask;
        if ~any(maskTrial), continue; end
        
        for si = 1:nsbj_SW
            m = maskTrial & visTable.subID == pairs_SW(si,1) & visTable.sessID == pairs_SW(si,2);
            if ~any(m), continue; end
            
            % 计算16-bin计数（不做去趋势）
            per_trial_counts(si, ti, :) = histcounts(mod(visTable.theta(m) + cfg.shift_FT, 360) - cfg.shift_FT, cfg.edges_FT);
        end
    end
    
    % ---- 3. 应用高斯滑动窗口（带边缘归一化） ----
    sigma = cfg.win_trials;
    radius = ceil(3 * sigma);
    x_kernel = -radius:radius;
    gauss_kernel = exp(-0.5 * (x_kernel / sigma).^2);
    
    smoothed_counts = zeros(size(per_trial_counts));
    
    for si = 1:nsbj_SW
        for bi = 1:16
            series = squeeze(per_trial_counts(si, :, bi));
            is_valid = series > 0 | any(per_trial_counts(si, :, :) > 0, 3); % 至少有一个bin有计数
            
            series_padded = series;
            series_padded(~is_valid) = 0;
            
            numerator = conv(series_padded, gauss_kernel, 'same');
            denominator = conv(double(is_valid), gauss_kernel, 'same');
            denominator(denominator < 1e-6) = 1; % 避免除以零，保持为0
            
            smoothed_series = numerator ./ denominator;
            smoothed_counts(si, :, bi) = smoothed_series;
        end
    end

    % ---- 4. 计算最终效应时序（占比） ----
    total_counts = sum(smoothed_counts, 3); % nsbj × total_trial
    total_counts(total_counts < eps) = NaN; % 避免除以零
    
    TseriesAxis = sum(smoothed_counts(:, :, 1:2:16), 3) ./ total_counts;
    TseriesGap  = sum(smoothed_counts(:, :, 2:2:16), 3) ./ total_counts;
    
    axis_counts = sum(smoothed_counts(:, :, 1:2:16), 3);
    axis_counts(axis_counts < eps) = NaN;
    
    TseriesCard = sum(smoothed_counts(:, :, 1:4:16), 3) ./ axis_counts;
    TseriesObli = sum(smoothed_counts(:, :, 3:4:16), 3) ./ axis_counts;
    
    xWin = 1:total_trial;
end
