function plot_pie(prop_array, color_matrix, varargin)
% 绘制饼图，根据占比数组和颜色矩阵
% prop_array: n x 1 数组, 每个扇区的占比 (和为1)
% color_matrix: n x 3 数组, 每行是RGB颜色 (0~1)
% 可选参数:
%   'show_text' (默认true): 是否在扇区中央显示数值
%   'merge_groups' (默认{}): cell数组，每个元素为需要合并显示的扇区索引向量
%   'ax' (默认新建figure): 指定绘图的axes句柄
%   'legend' (默认{}): cell数组，指定每个扇区的legend文本

    p = inputParser;
    addParameter(p, 'show_text', true, @(x)islogical(x) || isnumeric(x));
    addParameter(p, 'merge_groups', {}, @(x)iscell(x));
    addParameter(p, 'ax', [], @(x) isempty(x) || isgraphics(x,'axes'));
    addParameter(p, 'legend', {}, @(x)iscell(x) && (isempty(x) || numel(x)==numel(prop_array)));
    parse(p, varargin{:});
    show_text = logical(p.Results.show_text);
    merge_groups = p.Results.merge_groups;
    ax = p.Results.ax;
    legend_text = p.Results.legend;

    n = numel(prop_array);
    if size(color_matrix,1) ~= n || size(color_matrix,2) ~= 3
        error('color_matrix尺寸应为 n x 3');
    end

    if isempty(ax)
        figure;
        ax = axes;
    end
    axes(ax); cla(ax);

    theta_list = [0; cumsum(prop_array(:))*2*pi];
    hold on;
    h_fill = gobjects(n,1);
    for i = 1:n
        theta1 = theta_list(i);
        theta2 = theta_list(i+1);
        t = linspace(theta1, theta2, 100);
        x = [0, cos(t)];
        y = [0, sin(t)];
        h_fill(i) = fill(x, y, color_matrix(i,:), 'EdgeColor', 'none', 'FaceAlpha', 0.4);
    end
    axis equal off;

    if show_text
        shown = false(1,n); % 标记哪些扇区已显示
        % 先处理merge_groups
        for g = 1:numel(merge_groups)
            idx = merge_groups{g};
            if any(idx<1) || any(idx>n) || any(diff(idx)~=1)
                error('merge_groups中的索引必须为相邻且合法的扇区编号');
            end
            theta1 = theta_list(idx(1));
            theta2 = theta_list(idx(end)+1);
            theta_mid = (theta1 + theta2)/2;
            prop_sum = sum(prop_array(idx));
            text(cos(theta_mid)*0.6, sin(theta_mid)*0.6, sprintf('%.3f', prop_sum), ...
                'HorizontalAlignment','center','FontSize',18,'FontWeight','bold');
            shown(idx) = true;
        end
        % 其它未合并的扇区单独显示
        for i = 1:n
            if ~shown(i)
                theta1 = theta_list(i);
                theta2 = theta_list(i+1);
                theta_mid = (theta1 + theta2)/2;
                text(cos(theta_mid)*0.6, sin(theta_mid)*0.6, sprintf('%.3f', prop_array(i)), ...
                    'HorizontalAlignment','center','FontSize',18,'FontWeight','bold');
            end
        end
    end

    % 增加legend
    if ~isempty(legend_text)
        legend(h_fill, legend_text, 'Location', 'bestoutside');
    end

    hold off;
end
