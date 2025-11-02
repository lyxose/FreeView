function [effect_all, group_all, ver_names] = CrossV_collect_data(EffectStruct, compare_groups, do_diff)
    % 从EffectStruct中提取指定组的效应数据，准备进行交互效应分析
    if nargin < 2 || isempty(compare_groups)
        ver_names = fieldnames(EffectStruct);
    else
        ver_names = compare_groups(:)'; % 保持compare_groups顺序
    end
    if nargin < 3, do_diff = false; end

    effect_all = [];
    group_all = [];
    for i = 1:numel(ver_names)
        eff = EffectStruct.(ver_names{i});
        if do_diff
            eff = eff(:,1) - eff(:,2);
        elseif size(eff,2)==2
            eff = eff(:,1);
        end
        effect_all = [effect_all; eff(:)];
        group_all = [group_all; repmat(i, numel(eff(:)), 1)];
    end
end