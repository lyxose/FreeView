function [sub_time_bin_count, timeCenters_FT] = process_TimeCourse_count( ...
        angles_FT, start_FT, dur_FT, sub_FT, ses_FT, pairs_FT, nsbj_FT, ...
        timeRes_ms, keep_Time_ms, edges_FT, shift_FT, n_bin_FT)
% PROCESS_TIME_count_SERIES
%   将每个注视（angles_FT, start_FT, dur_FT）按时间轴汇总为“被试 × 时间 × 16-bin”的矩阵。
%
% 输入：
%   angles_FT     - 向量，所有注视的角度（度，0-360），已按 fixTable 筛选
%   start_FT      - 向量，对应每个注视的开始时间（ms）
%   dur_FT        - 向量，对应每个注视的持续时长（ms）
%   sub_FT, ses_FT- 向量，对应每个注视的被试ID与会话ID（用于分组）
%   pairs_FT      - nx2 矩阵，n = nsbj_FT，列为 [subID, sessID]（唯一被试×会话对）
%   nsbj_FT       - 标量，被试（会话）数量（rows(pairs_FT)）
%   timeRes_ms    - 标量，时间分辨率（ms），例如 10
%   keep_Time_ms  - 标量，时间范围上限（ms），例如 5000
%   edges_FT      - 16-bin 的角度边界（用于 histcounts）
%   shift_FT      - shift（度），用于把 0° 放在 bin 中心（与 edges_FT 配套）
%   n_bin_FT      - bin 数量（通常 16）
%
% 输出：
%   sub_time_bin_count - nsbj × nTime × n_bin (double) ：每被试每时间点每扇区的 detrended 计数
%   timeCenters_FT   - 1 x nTime 向量，返回各时间 bin 的中心（ms）
%
% 细节：
%   - 时间 bin 为左闭右开 [t, t+timeRes)，对应 histcounts 的时间分箱。
%   - 对于某一时间 bin，若某注视的 start<=t 且 start+dur>t，则认为该注视在该 bin 内“活动”。
%   - 每个时间 bin 对角度向量做 histcounts
%     这样可突出局部相对于周围扇区的偏差（与脚本其余部分一致）。
%
% 注意：该函数仅汇总角度维度，不做时间平滑。后续可在绘图端做时域平滑（卷积）。

    if nargin < 12
        error('process_time_count_series: missing required arguments');
    end

    % prepare time bins
    nTime = ceil(keep_Time_ms / timeRes_ms);
    timeEdges = (0:nTime) * timeRes_ms;          % edges: [0, 10, 20, ... , keep_Time]
    timeCenters_FT = timeEdges(1:end-1) + timeRes_ms/2;  % bin center e.g. 5,15,...

    % allocate output: nsbj x nTime x n_bin
    sub_time_bin_count = zeros(nsbj_FT, nTime, n_bin_FT);

    % loop per subject-session pair (keeps memory access localized)
    for si = 1:nsbj_FT
        % mask the indices for this subject-session
        m_sub = (sub_FT == pairs_FT(si,1)) & (ses_FT == pairs_FT(si,2));
        if ~any(m_sub)
            warning('process_time_count_series: no data for subject %d session %d', pairs_FT(si,1), pairs_FT(si,2));
            continue;
        end

        % extract vectors for this subject
        angs_all = angles_FT(m_sub);
        starts   = start_FT(m_sub);
        durs     = dur_FT(m_sub);

        % for each time bin, find active fixations and compute detrended 16-bin counts
        for ti = 1:nTime
            t_left = timeEdges(ti);   % bin left edge (inclusive)
            % active fixations: start <= t_left < start+dur
            active_mask = (starts <= t_left) & ((starts + durs) > t_left);
            if ~any(active_mask)
                continue;
            end
            angs_active = angs_all(active_mask);

            % histogram into the 16 bins (with circular shift applied)
            cnt = histcounts(mod(angs_active + shift_FT, 360) - shift_FT, edges_FT);

            % store (row = subject, col = time, page = bin)
            sub_time_bin_count(si, ti, :) = cnt;
        end
    end
end

