%% following the 
for p = select_sess
    % 加载数据
    sess = load(fullfile(dirs.mat, sprintf("Dat_Sub%.0f_Ses%.0f.mat", sub_ses_res(p,1), sub_ses_res(p,2))));
    
    resT = load(fullfile(dirs.fix, resfiles(p).name )).expT;
    data = resT.key2RT(52:92);
    data = data(~isnan(data));
    % 设置分组边界（强制最后一个bin宽度为1）
    binEdges = 0:1:16; % 关键修改：最后一个bin边界设为16
    adjusted_data = data;
    adjusted_data(data >=15) = 15; % 将所有≥15的值归入最后一个bin
    
    % 绘制直方图
    histogram(adjusted_data, 'BinEdges', binEdges, 'Normalization', 'probability');
    title(sprintf("Sub%.0f Ses%.0f %.3f min", sub_ses_res(p,1), sub_ses_res(p,2), double(sess.messages{end,1}-sess.messages{1,1})/1000000.0/60));
    
    % 计算并标记中位数（基于原始数据）
    median_val = median(data);
    hold on;
    xline(median_val, '--r', sprintf('Median=%.3f',median_val), 'LineWidth', 1.5, 'LabelVerticalAlignment', 'top');
    hold off;
    
    input('');
end