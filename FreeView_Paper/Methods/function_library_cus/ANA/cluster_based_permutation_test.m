function [stat, details] = cluster_based_permutation_test(seriesA, seriesB, varargin)
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

    p = inputParser;
    p.addParameter('n_perm', 5000, @(x) isnumeric(x)&&isscalar(x)&&x>0);
    p.addParameter('alpha', 0.05, @(x) isnumeric(x)&&isscalar(x)&&x>0&&x<1);
    p.addParameter('tail', 'both', @(x) ischar(x) && any(strcmpi(x,{'both','left','right'})));
    p.parse(varargin{:});
    n_perm = round(p.Results.n_perm);
    alpha = p.Results.alpha;
    tail = lower(p.Results.tail);

    seriesA = double(seriesA);
    seriesB = double(seriesB);
    diffMat = seriesA - seriesB;
    [nSubj, nTime] = size(diffMat);

    stat = struct(); % initialize as non-empty struct
    stat.prob = ones(1, nTime);
    stat.t_obs = [];
    stat.cluster_inds = {};
    stat.cluster_p = [];
    details = struct();

    if nSubj < 2 || nTime < 1
        return;
    end

    meanDiff = mean(diffMat, 1, 'omitnan');
    sdDiff = std(diffMat, 0, 1, 'omitnan');
    sdDiff(sdDiff < eps) = Inf;
    t_obs = meanDiff ./ (sdDiff ./ sqrt(sum(~isnan(diffMat),1)));
    stat.t_obs = t_obs;

    df = max(sum(~isnan(diffMat),1) - 1, 1);

    % Set threshold and cluster mask according to tail
    switch tail
        case 'both'
            t_thr = tinv(1 - alpha/2, max(df));
            above = abs(t_obs) > t_thr;
        case 'right'
            t_thr = tinv(1 - alpha, max(df));
            above = t_obs > t_thr;
        case 'left'
            t_thr = tinv(alpha, max(df));
            above = t_obs < t_thr;
        otherwise
            error('tail must be "both", "left", or "right".');
    end

    [obsStarts, obsEnds] = find_runs(above);
    % Cluster mass: sum of t values (abs for both-tail, signed for one-tail)
    if strcmp(tail,'both')
        obsMass = arrayfun(@(s,e) sum(abs(t_obs(s:e)),'omitnan'), obsStarts, obsEnds);
    else
        obsMass = arrayfun(@(s,e) sum(t_obs(s:e),'omitnan'), obsStarts, obsEnds);
    end
    stat.cluster_inds = arrayfun(@(s,e) [s,e], obsStarts, obsEnds, 'uni', false);

    perm_max = zeros(n_perm,1);
    signTemplate = ones(nSubj,1);
    for pp = 1:n_perm
        flips = signTemplate;
        flips(rand(nSubj,1) > 0.5) = -1;
        permDiff = diffMat .* flips;
        permMean = mean(permDiff, 1, 'omitnan');
        permSD = std(permDiff, 0, 1, 'omitnan');
        permSD(permSD < eps) = Inf;
        permT = permMean ./ (permSD ./ sqrt(sum(~isnan(permDiff),1)));
        switch tail
            case 'both'
                permAbove = abs(permT) > t_thr;
            case 'right'
                permAbove = permT > t_thr;
            case 'left'
                permAbove = permT < t_thr;
        end
        [ps, pe] = find_runs(permAbove);
        if ~isempty(ps)
            if strcmp(tail,'both')
                permMass = arrayfun(@(s,e) sum(abs(permT(s:e)),'omitnan'), ps, pe);
            else
                permMass = arrayfun(@(s,e) sum(permT(s:e),'omitnan'), ps, pe);
            end
            perm_max(pp) = max(permMass);
        else
            perm_max(pp) = 0;
        end
    end

    cluster_p = ones(1, numel(obsMass));
    for cc = 1:numel(obsMass)
        cluster_p(cc) = (1 + sum(perm_max >= obsMass(cc))) / (1 + n_perm);
        rng = obsStarts(cc):obsEnds(cc);
        stat.prob(rng) = min(stat.prob(rng), cluster_p(cc));
    end
    stat.cluster_p = cluster_p;

    details.t_threshold = t_thr;
    details.perm_max = perm_max;
    details.obs_clusters = [obsStarts(:), obsEnds(:), obsMass(:)];
    details.tail = tail;
end
