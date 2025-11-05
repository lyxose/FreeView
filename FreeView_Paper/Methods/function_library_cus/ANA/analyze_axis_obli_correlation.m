function analyze_axis_obli_correlation(Axis_Effect, Obli_Effect, varargin)
    % varargin: 可选控制变量结构体（如 Age, Gender, ...），每个结构体字段与Axis_Effect一致，值为nsubjx1
    % 支持参数 'doPerVersion'（逻辑），控制是否逐版本分析（默认false，合并分析）
    % 支持参数 'groupColors'（1x3 cell），指定每组的颜色（RGB向量），如 {[1 0 0],[0 1 0],[0 0 1]}
    % 支持参数 pairs_FTs（结构体，field name与版本一致，值为[nsubj x 2]的被试-会话对），用于标注被试编号
    % 合并分析时，跳过v1c，且v1_5数据取负号（翻转），保持配对关系
    % 可视化散点图、相关系数、回归线，输出统计结果

    % 1. 解析输入
    controlVars = {};
    controlNames = {};
    doPerVersion = false;
    groupColors = {[0.85 0.2 0.2],[0.2 0.5 0.85],[0.2 0.7 0.2]}; % 默认三组颜色
    pairs_FTs = [];
    for i = 1:numel(varargin)
        if isstruct(varargin{i}) && isfield(varargin{i}, 'v1') % 认为是pairs_FTs
            pairs_FTs = varargin{i};
        elseif isstruct(varargin{i})
            controlVars{end+1} = varargin{i};
            controlNames{end+1} = inputname(i+2);
        elseif ischar(varargin{i}) && strcmpi(varargin{i}, 'doPerVersion')
            if i+1 <= numel(varargin) && islogical(varargin{i+1})
                doPerVersion = varargin{i+1};
            else
                doPerVersion = true;
            end
        elseif islogical(varargin{i})
            doPerVersion = varargin{i};
        elseif ischar(varargin{i}) && strcmpi(varargin{i}, 'groupColors')
            if i+1 <= numel(varargin) && iscell(varargin{i+1})
                groupColors = varargin{i+1};
            end
        elseif iscell(varargin{i}) && numel(varargin{i})>=1 && size(varargin{i}{1},2)==3
            groupColors = varargin{i};
        end
    end

    % 2. 获取所有版本
    ver_names = intersect(fieldnames(Axis_Effect), fieldnames(Obli_Effect));
    ver_names = setdiff(ver_names, {'v1c'}); % 跳过v1c

    if doPerVersion
        % 3. 逐版本分析
        for vi = 1:numel(ver_names)
            ver = ver_names{vi};
            axis_vals = Axis_Effect.(ver)(:);
            obli_vals = Obli_Effect.(ver)(:,1) - Obli_Effect.(ver)(:,2);

            % v1_5翻转
            if strcmpi(ver, 'v1_5')
                axis_vals = -axis_vals;
                obli_vals = -obli_vals;
            end

            % 检查有效性
            valid_mask = isfinite(axis_vals) & isfinite(obli_vals);
            for ci = 1:numel(controlVars)
                ctrl = controlVars{ci}.(ver)(:);
                valid_mask = valid_mask & isfinite(ctrl);
            end

            axis_vals = axis_vals(valid_mask);
            obli_vals = obli_vals(valid_mask);
            nSubj = numel(axis_vals);

            % 4. 控制变量
            if numel(controlVars) > 0
                X_ctrl = zeros(nSubj, numel(controlVars));
                for ci = 1:numel(controlVars)
                    X_ctrl(:,ci) = controlVars{ci}.(ver)(valid_mask);
                end
                X = [axis_vals, X_ctrl];
                mdl = fitlm(X, obli_vals, 'Intercept', true);
                [r, p] = partialcorr(axis_vals, obli_vals, X_ctrl);
                fprintf('Version %s: Partial correlation (Axis vs Obli, controlling for %s): r = %.3f, p = %.4g\n', ...
                    map_labels(ver), strjoin(controlNames, ','), r, p);
            else
                [r, p] = corr(axis_vals, obli_vals, 'type', 'Pearson');
                mdl = fitlm(axis_vals, obli_vals, 'Intercept', true);
                fprintf('Version %s: Pearson correlation (Axis vs Obli): r = %.3f, p = %.4g\n', map_labels(ver), r, p);
            end

            % 5. 可视化
            figure;
            if vi <= numel(groupColors)
                thisColor = groupColors{vi};
            else
                thisColor = [0.5 0.5 0.5];
            end
            scatter(axis_vals, obli_vals, 40, thisColor, 'filled', 'MarkerFaceAlpha',0.7); hold on;

            % 被试编号标注（只取pairs_FTs第一列）
            if ~isempty(pairs_FTs) && isfield(pairs_FTs, ver)
                pairs_this = pairs_FTs.(ver);
                if size(pairs_this,1) == nSubj
                    for si = 1:nSubj
                        text(axis_vals(si), obli_vals(si), num2str(pairs_this(si,1)), ...
                            'Color', thisColor*0.7, 'FontSize', 9, 'FontWeight', 'bold', ...
                            'HorizontalAlignment','center', 'VerticalAlignment','middle');
                    end
                end
            end

            % 回归线与置信区间
            xfit = linspace(min(axis_vals), max(axis_vals), 100)';
            if numel(controlVars) > 0
                Xfit = [xfit, repmat(mean(X_ctrl,1), numel(xfit),1)];
                [yfit, yci] = predict(mdl, Xfit, 'Alpha', 0.05);
                yresid = obli_vals - predict(mdl, [axis_vals, X_ctrl]);
                yfit_std = std(yresid);
            else
                [yfit, yci] = predict(mdl, xfit, 'Alpha', 0.05);
                yresid = obli_vals - predict(mdl, axis_vals);
                yfit_std = std(yresid);
            end
            % 95%置信区间
            fill([xfit; flipud(xfit)], [yci(:,1); flipud(yci(:,2))], [0.7 0.7 0.7], ...
                'FaceAlpha',0.25, 'EdgeColor','none');
            % 回归线
            plot(xfit, yfit, 'k-', 'LineWidth',2);
            % 3倍标准残差范围
            plot(xfit, yfit+3*yfit_std, 'k--', 'LineWidth',1.2, 'Color',[0.5 0.5 0.5]);
            plot(xfit, yfit-3*yfit_std, 'k--', 'LineWidth',1.2, 'Color',[0.5 0.5 0.5]);

            % 显著性标注
            txt = sprintf('r = %.3f, p = %.4g', r, p);
            if numel(controlVars) > 0
                txt = [txt, sprintf('\nControls: %s', strjoin(controlNames, ','))];
            end
            text(mean(xfit), max(yfit), txt, 'FontSize',12, 'FontWeight','bold', ...
                'HorizontalAlignment','center', 'VerticalAlignment','bottom', 'Color','k');
            xlabel('Axis Effect'); ylabel('Obli Effect');
            title(sprintf('Axis vs Obli Effect — %s', ver));
            grid on; box off; hold off;
        end
    else
        % 合并分析
        axis_all = [];
        obli_all = [];
        group_all = [];
        ctrl_all = [];
        pairs_all = [];
        for vi = 1:numel(ver_names)
            ver = ver_names{vi};
            axis_vals = Axis_Effect.(ver)(:);
            obli_vals = Obli_Effect.(ver)(:,1) - Obli_Effect.(ver)(:,2);

            % v1_5翻转
            if strcmpi(ver, 'v1_5')
                axis_vals = -axis_vals;
                obli_vals = -obli_vals;
            end

            valid_mask = isfinite(axis_vals) & isfinite(obli_vals);
            for ci = 1:numel(controlVars)
                ctrl = controlVars{ci}.(ver)(:);
                valid_mask = valid_mask & isfinite(ctrl);
            end

            axis_vals = axis_vals(valid_mask);
            obli_vals = obli_vals(valid_mask);

            axis_all = [axis_all; axis_vals];
            obli_all = [obli_all; obli_vals];
            group_all = [group_all; repmat(vi, numel(axis_vals), 1)];
            % 控制变量合并
            for ci = 1:numel(controlVars)
                ctrl = controlVars{ci}.(ver)(valid_mask);
                if size(ctrl_all,2) < ci
                    ctrl_all(:,ci) = [];
                end
                ctrl_all = [ctrl_all; ctrl];
            end
            % 被试编号合并（只取pairs_FTs第一列）
            if ~isempty(pairs_FTs) && isfield(pairs_FTs, ver)
                pairs_this = pairs_FTs.(ver);
                pairs_all = [pairs_all; pairs_this(valid_mask,1)];
            end
        end

        nSubj = numel(axis_all);

        % 4. 控制变量
        if numel(controlVars) > 0
            X_ctrl = ctrl_all;
            X = [axis_all, X_ctrl];
            mdl = fitlm(X, obli_all, 'Intercept', true);
            [r, p] = partialcorr(axis_all, obli_all, X_ctrl);
            fprintf('Combined: Partial correlation (Axis vs Obli, controlling for %s): r = %.3f, p = %.4g\n', ...
                strjoin(controlNames, ','), r, p);
        else
            [r, p] = corr(axis_all, obli_all, 'type', 'Pearson');
            mdl = fitlm(axis_all, obli_all, 'Intercept', true);
            fprintf('Combined: Pearson correlation (Axis vs Obli): r = %.3f, p = %.4g\n', r, p);
        end

        % 5. 可视化
        figure;
        % 按组别着色
        for gi = 1:numel(ver_names)
            mask = group_all == gi;
            if gi <= numel(groupColors)
                thisColor = groupColors{gi};
            else
                thisColor = [0.5 0.5 0.5];
            end
            scatter(axis_all(mask), obli_all(mask), 40, thisColor, 'filled', 'MarkerFaceAlpha',0.7); hold on;
            % 被试编号标注（只取pairs_FTs第一列）
            if ~isempty(pairs_all)
                pairs_this = pairs_all(mask);
                xvals = axis_all(mask);
                yvals = obli_all(mask);
                for si = 1:numel(pairs_this)
                    text(xvals(si), yvals(si), num2str(pairs_this(si)), ...
                        'Color', thisColor*0.7, 'FontSize', 6, 'FontWeight', 'bold', ...
                        'HorizontalAlignment','center', 'VerticalAlignment','middle');
                end
            end
        end
        % 回归线与置信区间
        xfit = linspace(min(axis_all), max(axis_all), 100)';
        if numel(controlVars) > 0
            Xfit = [xfit, repmat(mean(X_ctrl,1), numel(xfit),1)];
            [yfit, yci] = predict(mdl, Xfit, 'Alpha', 0.05);
            yresid = obli_all - predict(mdl, [axis_all, X_ctrl]);
            yfit_std = std(yresid);
        else
            [yfit, yci] = predict(mdl, xfit, 'Alpha', 0.05);
            yresid = obli_all - predict(mdl, axis_all);
            yfit_std = std(yresid);
        end
        % 95%置信区间
        fill([xfit; flipud(xfit)], [yci(:,1); flipud(yci(:,2))], [0.7 0.7 0.7], ...
            'FaceAlpha',0.25, 'EdgeColor','none');
        % 回归线
        plot(xfit, yfit, 'k-', 'LineWidth',2);
        % 3倍标准残差范围
        plot(xfit, yfit+3*yfit_std, 'k--', 'LineWidth',1.2, 'Color',[0.5 0.5 0.5]);
        plot(xfit, yfit-3*yfit_std, 'k--', 'LineWidth',1.2, 'Color',[0.5 0.5 0.5]);

        txt = sprintf('r = %.3f, p = %.4g', r, p);
        if numel(controlVars) > 0
            txt = [txt, sprintf('\nControls: %s', strjoin(controlNames, ','))];
        end
        text(mean(xfit), max(yfit), txt, 'FontSize',12, 'FontWeight','bold', ...
            'HorizontalAlignment','center', 'VerticalAlignment','bottom', 'Color','k');
        xlabel('Axis Effect'); ylabel('Obli Effect');
        title('Axis vs Obli Effect — Combined Versions');
        grid on; box off; hold off;
        legend([map_labels(ver_names); {'95% CI'}], 'Location','best');
    end
end
