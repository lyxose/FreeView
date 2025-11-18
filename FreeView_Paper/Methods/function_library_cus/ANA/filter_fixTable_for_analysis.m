function [fixTable, start_FT, dur_FT, angles_FT, xpos_FT, ypos_FT, sub_FT, ses_FT, tri_FT, dnfix_FT, subj_stats] = filter_fixTable_for_analysis(fixTable, R_max)
    % ---- 清洗筛选 ----
    % 1. 剔除起始fixation在中心1°以内的fixation
    fixTable.dropFix = false(height(fixTable),1);
    fixTable.dNfix = fixTable.Nfix;
    fixTable.dlNfix = fixTable.lNfix;
    while true
        dropMask = fixTable.r<1 & fixTable.dNfix==1;
        if sum(dropMask)==0
            break;
        end
        fixTable.dropFix(dropMask) = true;
        changedTrials = fixTable.TriID(dropMask);
        changedTrialsMask = ismember(fixTable.TriID, changedTrials);
        fixTable.dNfix = fixTable.dNfix - changedTrialsMask;
        fixTable.dlNfix(dropMask) = 0;
    end
    % 2. 剔除最末连续几个位于目标1°以内的fixation
    while true
        dropMask =  fixTable.degErr<1 & fixTable.dlNfix==1;
        if sum(dropMask)==0
            break;
        end
        fixTable.dropFix(dropMask) = true;
        changedTrials = fixTable.TriID(dropMask);
        changedTrialsMask = ismember(fixTable.TriID, changedTrials);
        fixTable.dlNfix = fixTable.dlNfix - changedTrialsMask;
        fixTable.dNfix(dropMask) = 0;
    end
    % 3. 剔除r>R_max（背景范围）的fixation
    bgSelMask = fixTable.r>R_max+1;
    fixTable.dropFix(bgSelMask) = true;
    fixTable.dNfix(bgSelMask) = 0;
    fixTable.dlNfix(bgSelMask) = 0;
    fixTable.saccLen(bgSelMask) = 0;
    fixTable.saccAngle(bgSelMask) = 0;
    % 4. 剔除最末几个指向目标方向（误差连续减小）的fixation
    uTri = unique(fixTable.TriID);
    for t = transpose(uTri)
        idxTrial = find(fixTable.TriID==t & ~fixTable.dropFix);
        if numel(idxTrial) < 2
            continue
        end
        err = fixTable.degErr(idxTrial);
        dErr = diff(err);
        c = 0;
        for j = numel(dErr):-1:1
            if dErr(j) < 0
                c = c + 1;
            else
                break
            end
        end
        if c > 0
            dropIdxLocal = idxTrial(end-c+1:end);
            fixTable.dropFix(dropIdxLocal) = true;
        end
    end
    % 依据最终dropFix重算剔除后的序号
    fixTable.dNfix(:)  = 0;
    fixTable.dlNfix(:) = 0;
    for t = transpose(uTri)
        keepIdx = find(fixTable.TriID==t & ~fixTable.dropFix);
        n = numel(keepIdx);
        if n>0
            fixTable.dNfix(keepIdx)  = (1:n).';
            fixTable.dlNfix(keepIdx) = (n:-1:1).';
        end
    end
    % 输出变量
    start_FT  = fixTable.startT(~fixTable.dropFix);
    dur_FT    = fixTable.dur(~fixTable.dropFix);
    angles_FT = fixTable.theta(~fixTable.dropFix);
    xpos_FT   = fixTable.xpos(~fixTable.dropFix);
    ypos_FT   = fixTable.ypos(~fixTable.dropFix);
    sub_FT    = fixTable.subID(~fixTable.dropFix);
    ses_FT    = fixTable.sessID(~fixTable.dropFix);
    tri_FT    = fixTable.TriID(~fixTable.dropFix);
    dnfix_FT  = fixTable.dNfix(~fixTable.dropFix);

    % 统计各个被试的有效注视点数、试次数、最大Ntrial、trial跨度（该被试的最大Ntrial-最小Ntrial）
    subjIDs = unique(fixTable.subID(~fixTable.dropFix));
    nSubj = numel(subjIDs);
    subj_stats = table('Size',[nSubj,7], ...
        'VariableTypes',{'double','double','double','double','double','double','double'}, ...
        'VariableNames',{'subID','nFix','nTrial','maxNtrial','trialSpan','maxTrialTime','meanTrialTime'});
    for i = 1:nSubj
        sid = subjIDs(i);
        mask = fixTable.subID == sid & ~fixTable.dropFix;
        subj_stats.subID(i) = sid;
        subj_stats.nFix(i) = sum(mask);
        trials = unique(fixTable.Ntrial(mask));
        subj_stats.nTrial(i) = numel(trials);
        subj_stats.maxNtrial(i) = max(trials);
        subj_stats.trialSpan(i) = max(trials) - min(trials) + 1;
        % 计算该被试每个trial的时长（最后一个fixation的startT+dur）
        trial_times = zeros(numel(trials),1);
        for j = 1:numel(trials)
            tmask = mask & fixTable.Ntrial == trials(j);
            if any(tmask)
                last_idx = find(tmask,1,'last');
                trial_times(j) = fixTable.startT(last_idx) + fixTable.dur(last_idx);
            end
        end
        subj_stats.maxTrialTime(i) = max(trial_times)/1000;
        subj_stats.meanTrialTime(i) = mean(trial_times)/1000;
    end

    % 计算均值和标准差
    stat_cols = subj_stats.Properties.VariableNames(2:end); % 排除subID
    col_mean = zeros(1,numel(stat_cols));
    col_sd   = zeros(1,numel(stat_cols));
    for k = 1:numel(stat_cols)
        vals = subj_stats.(stat_cols{k});
        col_mean(k) = mean(vals, 'omitnan');
        col_sd(k)   = std(vals, 0, 'omitnan');
    end

    % 标记超出均值±3sd的单元格
    outlier_mask = false(size(subj_stats));
    for k = 1:numel(stat_cols)
        vals = subj_stats.(stat_cols{k});
        outlier_mask(:,k+1) = vals < (col_mean(k)-3*col_sd(k)) | vals > (col_mean(k)+3*col_sd(k));
    end
    % 输出表头（列名），并保证列对齐（每列宽度16，左对齐，空位补齐）
    col_names = stat_cols;
    nCols = numel(col_names);
    header = sprintf('%8s', 'subID');
    for k = 1:nCols
        header = [header, sprintf(' %-16s', col_names{k})];
    end
    disp(header);

    % 格式定义（每列宽度16，左对齐）
    fmt_nFix      = '%-16d';
    fmt_nTrial    = '%-16d';
    fmt_maxTime   = '%-16.2f';

    % 在异常值前后加星号*进行标记，左侧为小于3sd，右侧为大于3sd，且保证列对齐
    for i = 1:nSubj
        outStr = sprintf('%8d', subj_stats.subID(i));
        for k = 1:nCols
            val = subj_stats.(col_names{k})(i);
            switch k
                case 1 % nFix
                    val_fmt = fmt_nFix;
                case {2,3,4} % nTrial, maxNtrial, trialSpan
                    val_fmt = fmt_nTrial;
                case {5,6} % maxTrialTime, meanTrialTime
                    val_fmt = fmt_maxTime;
            end
            val_str = sprintf(val_fmt, val);
            if outlier_mask(i,k+1)
                if val < (col_mean(k)-3*col_sd(k))
                    val_str = ['*', strtrim(val_str)];
                    val_str = sprintf('%-16s', val_str);
                else
                    val_str = [strtrim(val_str), '*'];
                    val_str = sprintf('%-16s', val_str);
                end
            end
            outStr = [outStr, ' ', val_str];
        end
        % 补齐到16列（subID+7列数据+9空列）
        nTotalCols = 16;
        nPad = nTotalCols - (1 + nCols);
        if nPad > 0
            outStr = [outStr, repmat(sprintf(' %-16s', ''), 1, nPad)];
        end
        disp(outStr);
    end

    % 绘制所有被试的trial时长分布（每个trial最后一个fixation的startT+dur即为该trial时长），密度曲线（填充曲线下面积）用不同颜色、半透明地叠加在一张图上，在每个曲线最大峰的位置用小字标注被试编号。
end

