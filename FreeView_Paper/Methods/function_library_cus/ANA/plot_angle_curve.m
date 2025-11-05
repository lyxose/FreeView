function plot_angle_curve(data, centers, foldPeriod, binSize, CardColor, GapColor, ObliColor, normMode, varargin)
% ---- 扇区滑动绘图函数 ----
% data: Nsubj × nbins 原始数据
% centers: bin中心
% foldPeriod: 360/90/45
% binSize: bin宽度
% CardColor, GapColor, ObliColor: 颜色
% normMode: 标准化模式
% varargin: 可选参数，y轴范围 [ymin ymax]

    groupMean = mean(data, 1);
    groupSE = std(data, 0, 1) / sqrt(size(data, 1));
    nsbj_FT = size(data, 1);

    figure; hold on;
    fill([centers, fliplr(centers)], [groupMean+groupSE, fliplr(groupMean-groupSE)], ...
         [0.7 0.7 0.7], 'EdgeColor','none','FaceAlpha',0.35);
    plot(centers, groupMean, 'Color',[0.3 0.3 0.3], 'LineWidth',1.8);

    if ~isempty(varargin)
        yRange = varargin{1};
        if isnumeric(yRange) && numel(yRange)==2
            ylim(yRange);
        end
    end
    yl = ylim;

    xpad = max(1, 0.01 * (max(centers) - min(centers)));
    xlim([min(centers)-xpad, max(centers)+xpad]);

    if foldPeriod == 360
        axCols = [CardColor; ObliColor];
        xtickVals = 0:45:315;
        xticks(xtickVals);
        xticklabels(arrayfun(@(x) sprintf('%.1f°',mod(x,360)), xtickVals, 'uni', 0));
        for o = 0:7
            ang = o*45;
            if ang <= max(centers) && ang >= min(centers)
                line([ang,ang], yl+[0.01,-0.01], 'Color', axCols(mod(o,2)+1,:), 'LineStyle','--','LineWidth',0.8);
            end
        end
        for o = 0:7
            ang = o*45 + 22.5;
            if ang <= max(centers) && ang >= min(centers)
                line([ang,ang], yl+[0.01,-0.01], 'Color', GapColor, 'LineStyle','--','LineWidth',0.8);
            end
        end
        title(sprintf('Group Mean from fixTable (bin=%.1f°), n=%d', binSize, nsbj_FT));
    elseif foldPeriod == 45
        xtickVals = [0, 22.5];
        xticks(xtickVals);
        xticklabels(arrayfun(@(x) sprintf('%.1f°',x), xtickVals, 'uni', 0));
        line([0 0], yl+[0.01,-0.01], 'Color', CardColor, 'LineStyle','--','LineWidth',1.4);
        line([22.5 22.5], yl+[0.01,-0.01], 'Color', GapColor, 'LineStyle','--','LineWidth',1.4);
        title(sprintf('Group Mean from fixTable (folded 0-45°, bin=%.1f°), n=%d', binSize, nsbj_FT));
    elseif foldPeriod == 90
        xtickVals = [0, 22.5, 45, 67.5];
        xticks(xtickVals);
        xticklabels(arrayfun(@(x) sprintf('%.1f°', x), xtickVals, 'uni', 0));
        for a = [0, 22.5, 45, 67.5, 90]
            if a <= max(centers) && a >= min(centers)
                if a == 45
                    lc = ObliColor;
                elseif ismember(a, [22.5, 67.5])
                    lc = GapColor;
                else
                    lc = CardColor;
                end
                line([a a], yl + [0.01, -0.01], 'Color', lc, 'LineStyle', '--', 'LineWidth', 0.8);
            end
        end
        title(sprintf('Group Mean from fixTable (folded 0-90°, bin=%.1f°), n=%d', binSize, nsbj_FT));
        % 仅在90度条件下检验0度和45度的峰值差异显著性，并绘制显著性标记
        [~, idx0] = min(abs(centers - 0));
        [~, idx45] = min(abs(centers - 45));
        vals1 = data(:, idx0);
        vals2 = data(:, idx45);
        try
            [~, p_pair] = ttest(vals1, vals2);
        catch
            p_pair = NaN;
        end
        diffs = vals1 - vals2;
        d = mean(diffs) / std(diffs, 'omitnan');
        bump = 0.05 * range(ylim);
        y1 = groupMean(idx0) + groupSE(idx0) + bump;
        y2 = groupMean(idx45) + groupSE(idx45) + bump;
        y_bar = max([y1, y2, max(groupMean(idx0:idx45)+groupSE(idx0:idx45))]) + bump*0.5;
        color_sig = ternary(p_pair < 0.05, [1 0 0], [0 0 0]);
        plot([centers(idx0), centers(idx0), centers(idx45), centers(idx45)], ...
            [y1, y_bar, y_bar, y2], '-', 'Color', color_sig, 'LineWidth', 1.8);
        text(mean([centers(idx0), centers(idx45)]), y_bar + bump*0.6, ...
            sprintf('p=%.3f', p_pair), ...
            'HorizontalAlignment','center', 'VerticalAlignment','bottom', ...
            'FontSize',12, 'FontWeight','bold', 'Color', color_sig);
    end
    ylabel(['Count ', normMode]);
    xlabel('Angle (deg)');
    set(gca,'LineWidth',1.2); box off; hold off;
end
