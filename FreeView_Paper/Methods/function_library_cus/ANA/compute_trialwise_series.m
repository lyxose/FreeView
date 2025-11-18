function [TseriesAxis, TseriesGap, TseriesCard, TseriesObli, xWin] = compute_trialwise_series(visTable, cfg, total_trial)
% COMPUTE_TRIALWISE_SERIES 逐 trial 统计 Axis、Gap、Card、Obli 4组占比（不加窗）
%
% 输入：
%   visTable - fixation table
%   cfg      - 配置结构体，需包含:
%       .edges_FT, .shift_FT, .win_left, .win_right, .digPlace
%       .learn_stage_n  : 跳过的学习阶段 trial 数
%       .Card_in_Axis : Card 占比分母用 Axis（true）还是总数（false）
%   total_trial - 总 trial 数
%
% 输出：
%   Tseries... : nsbj × total_trial 矩阵，trialwise 占比
%   xWin       : 1 × total_trial 窗口中心位置向量 (1:total_trial)

    assert(istable(visTable) && ~isempty(visTable), 'visTable 必须是非空表格');
    requiredFields = {'edges_FT','shift_FT','win_left','win_right','digPlace'};
    for i = 1:numel(requiredFields)
        assert(isfield(cfg, requiredFields{i}), ['cfg 缺少字段: ', requiredFields{i}]);
    end
    if ~isfield(cfg, 'learn_stage_n'), cfg.learn_stage_n = 0; end
    if ~isfield(cfg, 'Card_in_Axis'), cfg.Card_in_Axis = false; end

    trialDigits = 10^(cfg.digPlace+1);
    sesTrial = mod(visTable.TriID, trialDigits);
    sesTrial(sesTrial==0) = trialDigits;

    timeOverlap = (visTable.startT < cfg.win_right) & ((visTable.startT + visTable.dur) > cfg.win_left);
    pairs_SW = unique([visTable.subID, visTable.sessID], 'rows', 'stable');
    nsbj_SW  = size(pairs_SW, 1);
    skip_n = max(0, cfg.learn_stage_n);
    sesTrial_adj = sesTrial - skip_n;
    valid_mask = sesTrial_adj > 0 & sesTrial_adj <= total_trial;

    per_trial_counts = zeros(nsbj_SW, total_trial, 16);
    for ti = 1:total_trial
        maskTrial = sesTrial_adj == ti & timeOverlap & valid_mask;
        if ~any(maskTrial), continue; end
        for si = 1:nsbj_SW
            m = maskTrial & visTable.subID == pairs_SW(si,1) & visTable.sessID == pairs_SW(si,2);
            if ~any(m), continue; end
            per_trial_counts(si, ti, :) = histcounts(mod(visTable.theta(m) + cfg.shift_FT, 360) - cfg.shift_FT, cfg.edges_FT);
        end
    end

    total_counts = sum(per_trial_counts, 3); % nsbj × total_trial
    axis_counts = sum(per_trial_counts(:, :, 1:2:16), 3);

    TseriesAxis = axis_counts ./ total_counts;
    TseriesGap  = sum(per_trial_counts(:, :, 2:2:16), 3) ./ total_counts;
    TseriesCard = sum(per_trial_counts(:, :, 1:4:16), 3) ./ ternary(cfg.Card_in_Axis, axis_counts, total_counts);
    TseriesObli = sum(per_trial_counts(:, :, 3:4:16), 3) ./ ternary(cfg.Card_in_Axis, axis_counts, total_counts);

    % 分母为0的地方填nan
    TseriesAxis(total_counts==0) = nan;
    TseriesGap(total_counts==0) = nan;
    if cfg.Card_in_Axis
        TseriesCard(axis_counts==0) = nan;
        TseriesObli(axis_counts==0) = nan;
    else
        TseriesCard(total_counts==0) = nan;
        TseriesObli(total_counts==0) = nan;
    end

    xWin = 1:total_trial;
end
