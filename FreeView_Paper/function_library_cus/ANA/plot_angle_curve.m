function plot_angle_curve(centers, groupMean, groupSE, foldPeriod, binSize, nsbj_FT, AxisColor, GapColor, ObliColor, normMode)
% ---- 扇区滑动绘图函数 ----
    % 绘制按角度滑动计数曲线，带组均值和SE阴影

    figure; hold on;
    fill([centers, fliplr(centers)], [groupMean+groupSE, fliplr(groupMean-groupSE)], ...
         [0.7 0.7 0.7], 'EdgeColor','none','FaceAlpha',0.35);
    plot(centers, groupMean, 'Color',[0.3 0.3 0.3], 'LineWidth',1.8);
    yl = ylim;

    % 自动调整x范围为数据范围稍大一点
    xpad = max(1, 0.01 * (max(centers) - min(centers)));
    xlim([min(centers)-xpad, max(centers)+xpad]);

    if foldPeriod == 360
        axCols = [AxisColor; ObliColor];
        % 设置xticks与axis对齐（-22.5°, 22.5°, 67.5°, ..., 337.5°）
        xtickVals = 0:45:315;
        xticks(xtickVals);
        xticklabels(arrayfun(@(x) sprintf('%.1f°',mod(x,360)), xtickVals, 'uni', 0));
        % 绘制axis和oblique的辅助线（0°, 45°, 90°, ...）
        for o = 0:7
            ang = o*45;
            if ang <= max(centers) && ang >= min(centers)
                line([ang,ang], yl+[0.01,-0.01], 'Color', axCols(mod(o,2)+1,:), 'LineStyle','--','LineWidth',0.8);
            end
        end
        % 绘制gap的辅助线（22.5°, 67.5°, ...）
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
        line([0 0], yl+[0.01,-0.01], 'Color', AxisColor, 'LineStyle','--','LineWidth',1.4);
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
                    lc = AxisColor;
                end
                line([a a], yl + [0.01, -0.01], 'Color', lc, 'LineStyle', '--', 'LineWidth', 0.8);
            end
        end
        title(sprintf('Group Mean from fixTable (folded 0-90°, bin=%.1f°), n=%d', binSize, nsbj_FT));
    end
    ylabel(['Count ', normMode]);
    xlabel('Angle (deg)');
    set(gca,'LineWidth',1.2); box off; hold off;
end
