function [TseriesAxis, TseriesGap, TseriesCard, TseriesObli, xWin] = compute_sliding_window_series(visTable, cfg, total_trial)
% COMPUTE_SLIDING_WINDOW_SERIES 对 visTable 执行基于 trial 的滑动窗口统计
%
% 输入：
%   visTable - fixation table（包含 subID, sessID, theta, startT, dur, TriID/Ntrial 等字段）
%   cfg      - 配置结构体，包含：
%       .edges_FT       : 16-bin 的角度边界向量
%       .shift_FT       : 角度偏移量（使 0° 落在 bin 中心）
%       .win_left       : 时间窗左边界（ms）
%       .win_right      : 时间窗右边界（ms）
%       .win_trials     : 滑动窗口大小（trials）
%       .step_trials    : 滑动步长（trials）
%       .total_trial    : 每会话总 trial 数（用于计算窗口数）
%       .digPlace       : TriID 中 trial 位数（用于提取会话内 trial 序号）
%       .mode           : 'continuous'（连续滑窗）或 'block'（按 block 对齐）
%       .blockSize      : 若 mode='block'，block 大小（如 96）
%       .learn_stage_n  : 跳过开头被舍弃的 trials 数
%       .doSmooth       : 是否对时间序列做平滑（默认 false）
%   total_trial - 传入总 trial 数（用于跳过最后被舍弃的 trials）
%
% 输出：
%   TseriesAxis, TseriesGap, TseriesCard, TseriesObli : nsbj × nWin 矩阵
%   xWin : 1 × nWin 窗口中心位置向量

    % ---- 参数验证与默认值 ----
    assert(istable(visTable) && ~isempty(visTable), 'visTable 必须是非空表格');
    requiredFields = {'edges_FT','shift_FT','win_left','win_right', ...
                      'win_trials','step_trials','total_trial','digPlace'};
    for i = 1:numel(requiredFields)
        assert(isfield(cfg, requiredFields{i}), ['cfg 缺少字段: ', requiredFields{i}]);
    end

    if ~isfield(cfg, 'mode'), cfg.mode = 'continuous'; end
    if ~isfield(cfg, 'blockSize'), cfg.blockSize = 96; end
    if ~isfield(cfg, 'learn_stage_n'), cfg.learn_stage_n = 0; end

    % ---- 提取会话内 trial 序号 ----
    trialDigits = 10^(cfg.digPlace+1);
    sesTrial = mod(visTable.TriID, trialDigits);
    sesTrial(sesTrial==0) = trialDigits;

    % ---- 时间窗筛选 ----
    timeOverlap = (visTable.startT < cfg.win_right) & ...
                  ((visTable.startT + visTable.dur) > cfg.win_left);

    % ---- 被试×会话列表 ----
    pairs_SW = unique([visTable.subID, visTable.sessID], 'rows', 'stable');
    nsbj_SW  = size(pairs_SW, 1);

    % ---- 跳过开头被舍弃的 trials ----
    skip_n = max(0, cfg.learn_stage_n);
    sesTrial_adj = sesTrial - skip_n; % 调整 trial 序号
    valid_mask = sesTrial_adj > 0 & sesTrial_adj <= total_trial; % 跳过开头和结尾被舍弃的 trials

    % ---- 根据模式设置窗口参数 ----
    switch lower(cfg.mode)
        case 'continuous'
            % 连续滑窗：1..total_trial
            nWin = total_trial - cfg.win_trials + 1;
            t0s  = 1:cfg.step_trials:nWin;
            t1s  = t0s + cfg.win_trials - 1;
            xWin = t0s + (cfg.win_trials - 1) / 2;
            nIterations = 1; % 只迭代一次

        case 'block'
            % 按 block 对齐滑窗
            nBlocks   = floor(total_trial / cfg.blockSize);
            nWinBlock = cfg.blockSize - cfg.win_trials + 1;
            if nWinBlock<=0
                warning('窗口大小 (%d)必须小于Block 大小 (%d) ', cfg.win_trials, cfg.blockSize);
                return
            end
            xWin      = (1:nWinBlock) + (cfg.win_trials - 1) / 2;
            nIterations = nBlocks; % 每个 block 一次迭代

        otherwise
            error('cfg.mode 必须是 "continuous" 或 "block"');
    end

    % ---- 结果容器 ----
    if strcmp(cfg.mode, 'continuous')
        TseriesAxis = zeros(nsbj_SW, numel(t0s));
        TseriesGap  = zeros(nsbj_SW, numel(t0s));
        TseriesCard = zeros(nsbj_SW, numel(t0s));
        TseriesObli = zeros(nsbj_SW, numel(t0s));
    else % 'block'
        blockAxis = nan(nsbj_SW, nWinBlock, nBlocks);
        blockGap  = nan(nsbj_SW, nWinBlock, nBlocks);
        blockCard = nan(nsbj_SW, nWinBlock, nBlocks);
        blockObli = nan(nsbj_SW, nWinBlock, nBlocks);
    end

    % ---- 主循环：逐窗口、逐被试统计 ----
    if strcmp(cfg.mode, 'continuous')
        bin_counts_all = zeros(nsbj_SW, numel(t0s), 16);
        for wi = 1:numel(t0s)
            t0 = t0s(wi); t1 = t1s(wi);
            maskWin = (sesTrial_adj >= t0) & (sesTrial_adj <= t1) & timeOverlap & valid_mask;
            if ~any(maskWin), continue; end
            for si = 1:nsbj_SW
                m = maskWin & visTable.subID==pairs_SW(si,1) & visTable.sessID==pairs_SW(si,2);
                if ~any(m), continue; end
                bin_counts_all(si, wi, :) = histcounts(mod(visTable.theta(m) + cfg.shift_FT, 360) - cfg.shift_FT, cfg.edges_FT);
            end
        end
        if cfg.doSmooth
            sigma = 3;
            bin_counts_all = gaussian_smooth_along_dim(bin_counts_all, sigma, 3);
        end
        TseriesAxis = sum(bin_counts_all(:,:,1:2:16), 3) ./ sum(bin_counts_all, 3);
        TseriesGap  = sum(bin_counts_all(:,:,2:2:16), 3) ./ sum(bin_counts_all, 3);
        TseriesCard = sum(bin_counts_all(:,:,1:4:16), 3) ./ sum(bin_counts_all(:,:,1:2:16), 3);
        TseriesObli = sum(bin_counts_all(:,:,3:4:16), 3) ./ sum(bin_counts_all(:,:,1:2:16), 3);
    else % 'block'
        for bi = 1:nBlocks
            base = (bi - 1) * cfg.blockSize;
            for t0 = 1:nWinBlock
                t1 = t0 + cfg.win_trials - 1;
                maskWin = (sesTrial_adj >= base + t0) & (sesTrial_adj <= base + t1) & timeOverlap & valid_mask;
                if ~any(maskWin), continue; end

                % 统计每个被试在当前窗口内的16-bin计数
                bin_counts = zeros(nsbj_SW, 16);
                for si = 1:nsbj_SW
                    m = maskWin & visTable.subID==pairs_SW(si,1) & visTable.sessID==pairs_SW(si,2);
                    if ~any(m), continue; end
                    
                    bin_counts(si,:) = histcounts(mod(visTable.theta(m) + cfg.shift_FT, 360) - cfg.shift_FT, cfg.edges_FT);
                end

                % 按主脚本方法计算各效应的占比
                blockAxis(si, t0, bi) = sum(bin_counts(si,1:2:16)) ./ sum(bin_counts(si,:));
                blockGap(si, t0, bi)  = sum(bin_counts(si,2:2:16)) ./ sum(bin_counts(si,:));
                blockCard(si, t0, bi) = sum(bin_counts(si,1:4:16)) ./ sum(bin_counts(si,1:2:16));
                blockObli(si, t0, bi) = sum(bin_counts(si,3:4:16)) ./ sum(bin_counts(si,1:2:16));
            end
        end

        % Block 模式：跨 block 平均
        TseriesAxis = squeeze(mean(blockAxis, 3, 'omitnan'));
        TseriesGap  = squeeze(mean(blockGap,  3, 'omitnan'));
        TseriesCard = squeeze(mean(blockCard, 3, 'omitnan'));
        TseriesObli = squeeze(mean(blockObli, 3, 'omitnan'));
    end
end
