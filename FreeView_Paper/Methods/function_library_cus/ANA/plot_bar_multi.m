function plot_bar_multi(data, colors, labels, varargin)
    % 多列输入的绘制（如16-bin），与原代码一致
    % data: nsbj x nBars
    % colors: nBars x 3
    % labels: cell array, nBars
    % varargin: 可选参数，支持 'ylabel', 'xlabel', 'title', 'showInd', 'showIndNum', 'baseline', 'baselineLabel'
    cfg = struct(varargin{:});
    [nSub, nB] = size(data);
    
    % 获取 baseline 参数，默认为 0
    if isfield(cfg, 'baseline')
        baseline = cfg.baseline;
    else
        baseline = 0;
    end
    
    % 获取 baselineLabel 参数，默认显示 baseline 的数值
    if isfield(cfg, 'baselineLabel')
        baselineLabel = cfg.baselineLabel;
    else
        baselineLabel = num2str(baseline);
    end
    
    if isempty(colors)
        colors = repmat([0.6 0.6 0.6], nB, 1);
    end
    m = mean(data,1,'omitnan');
    se = std(data,0,1,'omitnan') ./ sqrt(max(1, nSub));
    figure;
    hb = bar(m, 'FaceColor','flat', 'BaseValue', baseline);
    if size(colors,1) >= nB
        for ii = 1:nB
            hb.CData(ii,:) = colors(ii,:);
        end
    else
        hb.CData = repmat(colors(1,:), nB, 1);
    end
    hold on;
    xCenter = hb.XEndPoints;
    errorbar(xCenter, m, se, 'LineStyle','none','Color',[0.2 0.2 0.2],'LineWidth',1.2);
    
    if isfield(cfg,'xtickMode') && strcmp(cfg.xtickMode,'deg')
        maxTicks = 8;
        if nB > maxTicks
            step = ceil(nB / maxTicks);
            idx = 1:step:nB;
        else
            idx = 1:nB;
        end
        xticks(idx);
        xticklabels(labels(idx));
        xtickangle(0);
    elseif isfield(cfg,'xtickMode') && strcmp(cfg.xtickMode,'cat')
        xticks(1:nB);
        xticklabels(labels);
    else
        if nB <= 30
            xticks(1:nB);
            xticklabels(labels);
        else
            xticks(1:nB);
            step = ceil(nB / 30);
            idx = 1:step:nB;
            if idx(end) ~= nB, idx = [idx, nB]; end
            xticklabels(repmat({''}, 1, nB));
            xticklabels(idx, labels(idx));
        end
    end
    if isfield(cfg,'ylabel'), ylabel(cfg.ylabel); end
    if isfield(cfg,'xlabel'), xlabel(cfg.xlabel); end
    if isfield(cfg,'title'),  title(cfg.title);  end
    set(gca,'FontSize',12,'LineWidth',1.2); box off;
    bar_alpha_indv = 0.5;
    indYmax = nan(1,nB);
    
    if isfield(cfg,'showInd') && cfg.showInd
        hb.FaceAlpha = bar_alpha_indv;
        for ib = 1:nB
            x = xCenter(ib) * ones(nSub,1);
            scatter(x, data(:,ib), 20, repmat(colors(min(ib,size(colors,1)),:), nSub,1), ...
                'filled', 'MarkerFaceAlpha',0.6, 'MarkerEdgeAlpha',0.5);
            indYmax(ib) = max(data(:,ib));
            if isfield(cfg,'showIndNum') && cfg.showIndNum
                for si = 1:nSub
                    text(x(si), data(si,ib), num2str(si), 'FontSize',7, ...
                        'Color',[0.2 0.2 0.2], 'HorizontalAlign','center','VerticalAlign','middle');
                end
            end
        end
    end
    
    % 设置 y 轴刻度，确保 baseline 位置显示指定的标签
    yticks_current = get(gca, 'YTick');
    if ~ismember(baseline, yticks_current)
        yticks_new = sort([yticks_current, baseline]);
        set(gca, 'YTick', yticks_new);
    end
    yticklabels_current = get(gca, 'YTickLabel');
    yticks_vals = get(gca, 'YTick');
    baseline_idx = find(abs(yticks_vals - baseline) < eps*10);
    if ~isempty(baseline_idx)
        yticklabels_new = yticklabels_current;
        yticklabels_new{baseline_idx} = baselineLabel;
        set(gca, 'YTickLabel', yticklabels_new);
    end
    
    % 单样本 t 检验（对 baseline 做检验），并标注显著性符号
    pvals = nan(1, nB);
    for i = 1:nB
        try
            [~, pvals(i)] = ttest(data(:,i), baseline);
        catch
            pvals(i) = NaN;
        end
    end
    [~, ~, ~, adj_pvals] = fdr_bh(pvals, 0.05, 'pdep', 'yes');
    
    bump = diff(ylim) * 0.05; % 统一的bump，按y轴范围的5%
    for i = 1:nB
        s = sig_symbol(adj_pvals(i));
        if m(i) >= baseline
            y = m(i) + se(i) + bump;
            va = 'bottom';
        else
            y = m(i) - se(i) - bump;
            va = 'top';
        end
        text(xCenter(i), y, s, 'HorizontalAlign','center','VerticalAlign',va, ...
            'FontSize',12,'FontWeight','bold');
    end
    hold off;
end
