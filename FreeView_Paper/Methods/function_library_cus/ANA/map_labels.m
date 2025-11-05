function out_labels = map_labels(in_labels)
    label_map = containers.Map({'v1','v1_5','v2'}, {'v1.1','v1.2','v1.3'});
    label_map('v1.5') = 'v1.2';

    if ischar(in_labels)
        if isKey(label_map, in_labels)
            out_labels = label_map(in_labels);
        else
            out_labels = in_labels;
        end
    elseif iscell(in_labels) && numel(in_labels) == 1
        str = in_labels{1};
        if isKey(label_map, str)
            out_labels = label_map(str);
        else
            out_labels = str;
        end
    else
        out_labels = in_labels;
        for i = 1:numel(in_labels)
            if isKey(label_map, in_labels{i})
                out_labels{i} = label_map(in_labels{i});
            end
        end
    end
end