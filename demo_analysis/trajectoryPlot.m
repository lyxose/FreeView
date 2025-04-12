figure(2); hold on; axis equal;
colormap(parula);  % 选择颜色映射，例如 'parula', 'jet', 'hot'

maxTime = max(arrayfun(@(x) max(x.fix.startT), iresT.dat));  % 找到最大时间
caxis([0, maxTime]);  % 颜色范围基于时间

for i = 1:height(iresT)
    etime = iresT.dat(i).fix.endT;
    xpos = iresT.dat(i).fix.xpos(etime>0);
    ypos = iresT.dat(i).fix.ypos(etime>0);
    time = iresT.dat(i).fix.startT(etime>0);

    % 归一化时间值到 [0,1]，用于颜色映射
    normTime = time / maxTime;

    % 绘制渐变线段
    for j = 1:length(xpos)-1
        c = normTime(j);  % 当前时间的归一化值
        color = parula(256);  % 获取 colormap 颜色
        idx = max(1, round(c * 255));  % 找到颜色索引
        plot(xpos(j:j+1), ypos(j:j+1), '-', 'Color', color(idx, :), 'LineWidth', 1.5);
    end
end

colorbar;  % 显示颜色条
title('Fixation Path with Time Gradient');
xlabel('X Position');
ylabel('Y Position');
