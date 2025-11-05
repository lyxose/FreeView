function [sub_fix_bin_count, nfix_FT] = process_fixation_count_series( ...
    angles_FT, dnfix_FT, sub_FT, ses_FT, pairs_FT, nsbj_FT, keep_nFix, ...
    n_bin_FT, shift_FT, edges_FT)
% 按注视序号聚合角度数据为每被试 × 注视序号 × 16-bin 的计数矩阵
%
% 输入参数：
%   angles_FT   - 所有保留注视的角度（度，0-360）
%   dnfix_FT    - 每注视在 trial 中的剔除后序号（1..n）
%   sub_FT      - 被试ID（与 angles_FT 对齐）
%   ses_FT      - 会话ID（与 angles_FT 对齐）
%   pairs_FT    - nsbj_FT × 2，被试-会话对列表
%   nsbj_FT     - 有效被试数
%   keep_nFix   - 希望显示的最大注视序号
%   cmap16_FT   - 16×3颜色映射
%   n_bin_FT    - bin数量（通常16）
%   shift_FT    - bin中心偏移（deg）
%   edges_FT    - bin边界（17元素）
%
% 输出参数：
%   sub_fix_bin_count - nsbj_FT × nfix_FT × 16，每被试每注视序号的16-bin计数
%   nfix_FT        - 实际用于计算的注视序号上限
%
% 说明：
%   - 对每个被试-会话，按注视序号统计角度分布
%   - 结果用于时序分析和可视化

    % 上限注视数（至少 1）
    nfix_FT = max(1, round(keep_nFix));

    % 预分配容器：nsbj × nfix × 16
    sub_fix_bin_count = zeros(nsbj_FT, nfix_FT, n_bin_FT);

    % 对每个被试-会话遍历并统计
    for si = 1:nsbj_FT
        % 掩码：选出该被试会话的所有注视
        m_sub = (sub_FT == pairs_FT(si,1)) & (ses_FT == pairs_FT(si,2));
        if ~any(m_sub)
            warning('process_fixation_count_series: no data for subject %d session %d', pairs_FT(si,1), pairs_FT(si,2));
            continue;
        end

        % 被试的注视序号（剔除后）与角度
        dn = dnfix_FT(m_sub);      % 1..n
        angs = angles_FT(m_sub);   % deg 0..360

        % 每个注视序号 k 单独统计 16-bin 计数并做去趋势
        for k = 1:nfix_FT
            mk = (dn == k); % 逻辑掩码：该被试会话中所有 oblique==k 的注视点
            if ~any(mk)
                % 保持为 0
                continue;
            end
            cnt = histcounts( mod(angs(mk) + shift_FT, 360) - shift_FT, edges_FT );
            sub_fix_bin_count(si, k, :) = cnt;
        end
    end

    % 返回
    % sub_fix_bin_count: nsbj x nfix x 16
end

