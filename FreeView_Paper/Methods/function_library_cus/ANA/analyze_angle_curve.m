function [centers, groupMean, groupSE, subjCurves] = analyze_angle_curve(angles_FT, sub_FT, ses_FT, pairs_FT, nsbj_FT, binSize, foldPeriod, normMode,startAngle)
    % 对每个被试(会话)的所有 angles_FT 进行沿角度轴的滑动计数
    % foldPeriod: 360 或 45 或 90
    % binSize: bin 大小（deg）
    % normMode: 'zscore'|'minmax'|'sum1'|'demean'|'mean1'
    % startAngle: 起始角度偏移（deg），决定第一个 bin 的中心位置
    % 
    % 从 0-binSize/2 开始扫描，横轴范围为 [-binSize/2, foldPeriod-binSize/2]
    % 第一个数据点统计角度落在 (360-binSize/2, binSize/2] 范围内的点（跨越360°边界）
    
    
    % 确定扫描分辨率
    if foldPeriod == 360
        step = 0.01;
    else
        step = 0.1;
    end
    
    % 新的起始位置：0 - binSize/2
    % startAngle = -22.5/2;
    % 扫描范围：从 -binSize/2 到 foldPeriod - binSize/2
    centers = startAngle + (0:round(foldPeriod/step)-1)*step;
    nbins = numel(centers);
    subjCurves = zeros(nsbj_FT, nbins);

    for si = 1:nsbj_FT
        m = sub_FT==pairs_FT(si,1) & ses_FT==pairs_FT(si,2);
        angs = angles_FT(m);
        
        % 对每个 bin 中心，计算落在 [center-binSize/2, center+binSize/2) 内的点数
        for ii = 1:nbins
            center = centers(ii);
            leftEdge = center - binSize/2;
            rightEdge = center + binSize/2;
            
            % 处理环形周期性：将角度归一化到 [leftEdge, leftEdge+foldPeriod)
            % 这样可以正确处理跨越边界的情况
            angs_shifted = mod(angs - leftEdge, foldPeriod);
            
            % 统计落在 [0, binSize) 范围内的点（即原始 [leftEdge, rightEdge) 范围）
            subjCurves(si, ii) = sum(angs_shifted < binSize);
        end
    end

    % 被试内归一化
    if strcmpi(normMode, 'sum1')
        for si = 1:nsbj_FT
            m = sub_FT==pairs_FT(si,1) & ses_FT==pairs_FT(si,2);
            subjCurves(si,:) = subjCurves(si,:) ./ sum(m);
        end
    else
        subjCurves = normalize_by_dim(subjCurves, normMode);
    end

    groupMean = mean(subjCurves,1);
    groupSE   = std(subjCurves,0,1)/sqrt(nsbj_FT);
end
