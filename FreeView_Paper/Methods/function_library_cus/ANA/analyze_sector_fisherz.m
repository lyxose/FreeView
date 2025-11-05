function [stat_axis, stat_card, time_ms, fz_axis_obs, fz_card_obs, perm_z8, perm_z4] = analyze_sector_fisherz(sector_counts_sub, frame_time_res)
    % 对 sector_counts_sub 做扇区相关性分析和簇置换检验
    % 返回 Axis/Card 的统计结构体和相关序列
    nT   = size(sector_counts_sub,1);
    nSec = size(sector_counts_sub,2);
    nSub = size(sector_counts_sub,3);
    time_ms = (0:nT-1) * frame_time_res;
    idx_axis = 1:2:16;
    idx_card = 1:4:16;

    % 高斯平滑
    smooth_sigma_ms = 100;
    sigma_frames = max(0.5, smooth_sigma_ms / frame_time_res);
    radius = ceil(3 * sigma_frames);
    tvec = -radius:radius;
    gk = exp(-0.5 * (tvec ./ sigma_frames).^2);
    gk = gk / sum(gk);

    data2d = reshape(single(sector_counts_sub), nT, nSec * nSub);
    sm_data2d = zeros(size(data2d), 'single');
    for col = 1:size(data2d,2)
        sm_data2d(:,col) = conv(data2d(:,col), gk, 'same');
    end
    sector_counts_sub_sm = reshape(sm_data2d, nT, nSec, nSub);

    % 协方差矩阵
    S_all = zeros(nSec, nSec, nT, 'double');
    for t = 1:nT
        Xt = double(squeeze(sector_counts_sub_sm(t,:,:))');
        S_all(:,:,t) = cov(Xt, 1);
    end
    compute_fisherz_series = @(idx_set) arrayfun(@(tt) fisherz_from_cov(S_all(idx_set,idx_set,tt)), 1:nT);

    fz_axis_obs = compute_fisherz_series(idx_axis);
    fz_card_obs = compute_fisherz_series(idx_card);

    % 置换组合
    all8  = nchoosek(1:16, 8);
    row_axis = ismember(all8, idx_axis, 'rows');
    all8_no_axis = all8(~row_axis,:);
    n_all8 = size(all8_no_axis,1);
    n_perm8 = min(5000, n_all8);
    pick8 = all8_no_axis(randperm(n_all8, n_perm8), :);

    all4  = nchoosek(1:16, 4);
    row_card = ismember(all4, idx_card, 'rows');
    all4_no_card = all4(~row_card,:);
    n_perm4 = size(all4_no_card,1);
    pick4 = all4_no_card;

    perm_z8 = nan(n_perm8, nT);
    for p = 1:n_perm8
        idx = pick8(p,:);
        for tt = 1:nT
            perm_z8(p,tt) = fisherz_from_cov(S_all(idx,idx,tt));
        end
    end
    perm_z4 = nan(n_perm4, nT);
    for p = 1:n_perm4
        idx = pick4(p,:);
        for tt = 1:nT
            perm_z4(p,tt) = fisherz_from_cov(S_all(idx,idx,tt));
        end
    end

    thr95_axis = prctile(perm_z8, 95, 1);
    thr95_card = prctile(perm_z4, 95, 1);

    stat_axis = cluster_perm_fisherz(fz_axis_obs, perm_z8, thr95_axis);
    stat_card = cluster_perm_fisherz(fz_card_obs, perm_z4, thr95_card);
end


function z = fisherz_from_cov(Ssub)
    % 从协方差子矩阵 Ssub (k x k) 计算平均相关 r 的 Fisher z
    k = size(Ssub,1);
    if k < 2 || any(~isfinite(Ssub(:)))
        z = NaN; return;
    end
    d = sqrt(max(diag(Ssub), 0));
    denom = d*d.';
    R = Ssub ./ denom;
    R(denom==0) = NaN;
    mask = triu(true(k),1);
    rvals = R(mask);
    r_mean = mean(rvals, 'omitnan');
    if ~isfinite(r_mean)
        z = NaN; return;
    end
    r_mean = min(max(r_mean, -0.999999), 0.999999);
    z = atanh(r_mean);
end


function out = cluster_perm_fisherz(obs_z, perm_z, thr95)
    % 基于阈值向量 thr95（每时刻置换95%分位）的一侧簇检验，簇质量=簇内 Fisher z 之和
    [nPerm, ~] = size(perm_z);
    mu_perm = mean(perm_z, 1, 'omitnan');
    sd_perm = std(perm_z, 0, 1, 'omitnan');
    sd_perm(sd_perm<eps) = eps;
    % 单点p（经验分布，单侧：perm >= obs）
    p_time = (1 + sum(bsxfun(@ge, perm_z, obs_z), 1)) ./ (1 + nPerm);

    % 观测簇
    mask = obs_z >= thr95;
    [sIdx, eIdx] = find_runs(mask);
    nC = numel(sIdx);
    obs_mass = zeros(nC,1);
    for c = 1:nC
        rng = sIdx(c):eIdx(c);
        obs_mass(c) = nansum(obs_z(rng));
    end

    % 置换最大簇质量分布
    perm_max_mass = zeros(nPerm,1);
    for p = 1:nPerm
        m = perm_z(p,:) >= thr95;
        [ss, ee] = find_runs(m);
        mm = 0;
        for c = 1:numel(ss)
            rng = ss(c):ee(c);
            mm = max(mm, nansum(perm_z(p,rng)));
        end
        perm_max_mass(p) = mm;
    end

    % 簇p值（FWER）
    p_cluster = nan(nC,1);
    for c = 1:nC
        p_cluster(c) = (1 + sum(perm_max_mass >= obs_mass(c))) / (1 + nPerm);
    end

    out = struct();
    out.mu_perm   = mu_perm;
    out.sd_perm   = sd_perm;
    out.p_time    = p_time;
    out.clu_start = sIdx;
    out.clu_end   = eIdx;
    out.clu_mass  = obs_mass;
    out.clu_p     = p_cluster;
    out.perm_maxClu = perm_max_mass;
    out.thr95     = thr95;
end

function [sIdx, eIdx] = find_runs(mask)
    d = diff([false, mask, false]);
    sIdx = find(d==1);
    eIdx = find(d==-1)-1;
end

