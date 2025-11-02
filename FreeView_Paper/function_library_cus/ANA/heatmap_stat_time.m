function plot_prop_violin(data, colors, labels, varargin)
    % 数据为比例，展示分布特征
    % data: nsbj x 1 或 nsbj x 2
    % colors: nBars x 3
    % labels: cell array, nBars
    % varargin: 可选参数，支持 'ylabel', 'xlabel', 'title', 'showInd', 'showIndNum', 'showIndLink', 'chanceLevel'
    cfg = struct(varargin{:});
    [nSub, nB] = size(data);
    if isempty(colors)
        colors = repmat([0.6 0.6 0.6], nB, 1);
    end
    if isfield(cfg,'chanceLevel')
        chanceLevel = cfg.chanceLevel;
    else
        chanceLevel = 0.5;
    end

    % --- 3sd异常值筛选 ---
    outlier_mask = false(size(data));
    for ib = 1:nB
        vals = data(:,ib);
        mu = mean(vals,'omitnan');
        sd = std(vals,'omitnan');
        outlier_mask(:,ib) = abs(vals-mu) > 3*sd;
    end
    data_valid = data;
    data_valid(outlier_mask) = NaN;

    m = mean(data_valid,1,'omitnan');
    se = std(data_valid,0,1,'omitnan') ./ sqrt(max(1, sum(~isnan(data_valid),1)));

    figure;
    hold on;
    x = 1:nB;

    % 绘制chance level基线
    yline(chanceLevel, '--', 'Color', [0.4 0.4 0.4], 'LineWidth', 1.2, 'Label', 'Gap', 'LabelHorizontalAlignment','left');

    % 绘制小提琴图
    violin_width = 0.3;
    indYmax = nan(1,nB);
    xi_all = [];
    for ib = 1:nB
        vals = data_valid(:,ib);
        vals = vals(~isnan(vals));
        if isempty(vals), continue; end

        % 核密度估计
        [f, xi] = ksdensity(vals, 'BoundaryCorrection','reflection', 'Support',[0 1]);
        f = f / max(f) * violin_width; % 归一化宽度
        xi_all = [xi_all, xi(:)];

        % 绘制小提琴形状（左右对称）
        patch([x(ib) - f, x(ib) + fliplr(f)], [xi, fliplr(xi)], colors(ib,:), ...
            'FaceAlpha', 0.3, 'EdgeColor', colors(ib,:), 'LineWidth', 1);

        % 绘制箱线图元素
        q = prctile(vals, [25 50 75]);
        iqr = q(3) - q(1);
        whisker_low = max(min(vals), q(1) - 1.5*iqr);
        whisker_high = min(max(vals), q(3) + 1.5*iqr);

        % 须
        plot([x(ib) x(ib)], [whisker_low q(1)], 'Color', colors(ib,:)*0.6, 'LineWidth', 1.2);
        plot([x(ib) x(ib)], [q(3) whisker_high], 'Color', colors(ib,:)*0.6, 'LineWidth', 1.2);

        % 箱体
        box_w = 0.08;
        rectangle('Position', [x(ib)-box_w/2, q(1), box_w, q(3)-q(1)], ...
            'FaceColor', [1 1 1 0.7], 'EdgeColor', colors(ib,:)*0.6, 'LineWidth', 1.5);

        % 中位数线
        plot([x(ib)-box_w/2 x(ib)+box_w/2], [q(2) q(2)], ...
            'Color', colors(ib,:)*0.4, 'LineWidth', 2);

        indYmax(ib) = max(vals);
    end

    % individual 数据点
    indYall = [];
    if isfield(cfg,'showInd') && cfg.showInd
        jitter_amp = 0.12;
        Xj = zeros(nSub,nB);
        for ib = 1:nB
            xj = x(ib) + (rand(nSub,1)-0.5)*jitter_amp;
            Xj(:,ib) = xj;
            % 正常点
            mask_in = ~outlier_mask(:,ib) & ~isnan(data(:,ib));
            scatter(xj(mask_in), data(mask_in,ib), 20, repmat(colors(ib,:), sum(mask_in),1), ...
                'filled', 'MarkerFaceAlpha',0.6, 'MarkerEdgeAlpha',0.4);
            % 异常点
            mask_out = outlier_mask(:,ib) & ~isnan(data(:,ib));
            scatter(xj(mask_out), data(mask_out,ib), 40, [0.5 0.5 0.5], 'x', ...
                'LineWidth',1.8,'MarkerEdgeAlpha',0.7);
            indYall = [indYall; data(mask_in,ib); data(mask_out,ib)];
        end

        if nB == 2 && isfield(cfg,'showIndLink') && cfg.showIndLink
            for si = 1:nSub
                if any(outlier_mask(si,:)), continue; end % 跳过有异常值的连线
                plot(Xj(si,:), data(si,:), '-', 'Color', [0.7 0.7 0.7], 'LineWidth', 0.6, 'HandleVisibility', 'off');
            end
        end

        if isfield(cfg,'showIndNum') && cfg.showIndNum
            for ib = 1:nB
                mask_in = ~outlier_mask(:,ib) & ~isnan(data(:,ib));
                for si = find(mask_in)'
                    text(Xj(si,ib), data(si,ib), num2str(si), 'FontSize',7, ...
                        'Color',[0.2 0.2 0.2], 'HorizontalAlign','center','VerticalAlign','middle');
                end
            end
        end
        y_span = ylim;
        bump = (y_span(2) - y_span(1)) * 0.05;
        % 显著性标注（组间配对 t 检验）
        if nB == 2
            mask_pair = ~any(outlier_mask,2) & all(~isnan(data),2);
            try
                [~, p_pair] = ttest(data(mask_pair,1), data(mask_pair,2));
            catch
                p_pair = NaN;
            end

            y1 =  max(xi_all(:,1)) + bump;
            y2 =  max(xi_all(:,2)) + bump;
            y_bar = max([y1, y2]) + bump*0.5;
            color_sig = ternary(p_pair < 0.05, [1 0 0], [0 0 0]);
            plot([x(1), x(1), x(2), x(2)], [y1, y_bar, y_bar, y2], '-', 'Color', color_sig, 'LineWidth', 1.8);
            text(mean(x), y_bar + bump*0.6, sprintf('p=%.3f', p_pair), ...
                'HorizontalAlignment','center', 'VerticalAlignment','bottom', ...
                'FontSize',12, 'FontWeight','bold', 'Color', color_sig);
        end
    end

    % 单样本 t 检验（对 chanceLevel 做检验），并标注显著性符号
    pvals = nan(1, nB);
    for i = 1:nB
        mask_in = ~outlier_mask(:,i) & ~isnan(data(:,i));
        try
            [~, pvals(i)] = ttest(data(mask_in,i), chanceLevel);
        catch
            pvals(i) = NaN;
        end
    end
    [~, ~, ~, adj_pvals] = fdr_bh(pvals, 0.05, 'pdep', 'yes');
    for i = 1:nB
        s = sig_symbol(adj_pvals(i));
        % y轴位置：小提琴顶部再往上偏移一点
        y = max(xi_all(:,i)) + bump;
        text(x(i), y, s, 'HorizontalAlign','center','VerticalAlign','bottom',...
            'FontSize',12,'FontWeight','bold');
    end

    % xticks/labels
    xticks(x);
    xticklabels(labels);
    if isfield(cfg,'ylabel'), ylabel(cfg.ylabel); end
    if isfield(cfg,'xlabel'), xlabel(cfg.xlabel); end
    if isfield(cfg,'title'),  title(cfg.title);  end
    set(gca,'FontSize',12,'LineWidth',1.2); box off;

    % 自动y轴范围，保证所有点（包括异常点）都在画面内
    y_min = min([chanceLevel, min(m-se), min(data_valid(:)), min(xi_all(:)), min(data(:))]);
    y_max = max([chanceLevel, max(m+se), max(data_valid(:)), max(xi_all(:)), max(data(:))]);
    if nB == 2 && exist('y_bar','var')
        y_max = max(y_max, y_bar + bump*1.6);
    end
    y_pad = 0.08 * (y_max - y_min + eps);
    ylim([max(0, y_min - y_pad), min(1, y_max + y_pad)]);
    hold off;
end