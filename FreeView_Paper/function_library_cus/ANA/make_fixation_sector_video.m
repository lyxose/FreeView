function make_fixation_sector_video(start_FT, dur_FT, xpos_FT, ypos_FT, win_left, win_right, frame_rate, video_width, video_height, R_max, cmap16_FT, sector_edges, ut)
% 生成注视点视频，背景为16扇区，每帧显示当前注视点及扇区计数
% 输入参数：
%   start_FT, dur_FT, xpos_FT, ypos_FT : 注视点的开始时间、持续时间、X、Y坐标（向量，单位ms/pixel）
%   win_left, win_right                : 视频时间窗（ms）
%   frame_rate                         : 帧率（fps）
%   video_width, video_height          : 视频尺寸（pixel）
%   R_max                              : 扇区半径（deg）
%   cmap16_FT                          : 16×3颜色映射（可选，默认hsv）
%   sector_edges                       : 17元素扇区边界（可选）

if nargin < 11 || isempty(cmap16_FT) || size(cmap16_FT,1)~=16
    cmap16_FT = hsv(16);
end
if nargin < 12 || isempty(sector_edges) || numel(sector_edges)~=17
    sector_edges = linspace(-11.25,360-11.25,17);
end

frame_time_res = 1000 / frame_rate;
tVec = win_left:frame_time_res:win_right;
n_frames = numel(tVec);
center_pix = [video_width, video_height] / 2;
light_cmap = min(1, 0.15 + 0.85*cmap16_FT);
RpatchPix = max(1, round(ut.deg2pix(R_max)));
theta_patch_res = 80;
sectorPolys = cell(16,1);
labelPos = zeros(16,2);
for si = 1:16
    a1 = sector_edges(si);
    a2 = sector_edges(si+1);
    if a2 < a1, a2 = a2 + 360; end
    aa = linspace(a1, a2, theta_patch_res);
    xp = center_pix(1) + RpatchPix * cosd(aa);
    yp = center_pix(2) + RpatchPix * sind(aa);
    sectorPolys{si} = [center_pix; [xp(:), yp(:)]];
    midA = mod((a1 + (a2 - a1)/2), 360);
    labR = RpatchPix * 1.08;
    labelPos(si,:) = center_pix + labR * [cosd(midA), sind(midA)];
end

vidName = sprintf('FixPoints_%d_%dms_%.0ffps_Rmax%.1fdeg.avi', win_left, win_right, frame_rate, R_max);
vw_fp = VideoWriter(vidName);
vw_fp.FrameRate = frame_rate;
open(vw_fp);

fig = figure('Visible','off','Position',[100 100 video_width video_height]);
ax = axes(fig,'Position',[0 0 1 1]);
axis(ax,[1 video_width 1 video_height]);
set(ax,'YDir','reverse','XTick',[],'YTick',[]); hold(ax,'on'); axis(ax,'ij');

for si = 1:16
    P = sectorPolys{si};
    patch('XData',P(:,1),'YData',P(:,2), ...
          'FaceColor',light_cmap(si,:), ...
          'EdgeColor','none','FaceAlpha',0.28,'Parent',ax);
end
for si = 1:16
    lbl = sprintf('S%02d', si);
    text(labelPos(si,1), labelPos(si,2), lbl, ...
         'Parent', ax, 'HorizontalAlignment','center', 'VerticalAlignment','middle', ...
         'FontWeight','bold', 'FontSize', 10, 'Color', [0.06 0.06 0.06]);
end

hScat = scatter(ax, nan, nan, 30, 'k','filled','MarkerFaceAlpha',0.85,'MarkerEdgeAlpha',0.85);
hTxt  = text(ax, 15, 30, '', 'Parent',ax, 'Color','k','FontSize',12,'FontWeight','bold', ...
             'VerticalAlignment','top','BackgroundColor',[1 1 1 0.65], 'Margin',4);
for fi = 1:n_frames
    tNow = tVec(fi);
    active = (start_FT <= tNow) & (start_FT + dur_FT > tNow);
    if any(active)
        xA = xpos_FT(active);
        yA = ypos_FT(active);
        angA = mod(atan2d(yA - center_pix(2), xA - center_pix(1)), 360);
        sectorIdx = zeros(numel(angA),1,'uint8');
        for si = 1:16
            l = mod(sector_edges(si),360);
            r = mod(sector_edges(si+1),360);
            if l < r
                mask = angA >= l & angA < r;
            else
                mask = angA >= l | angA < r;
            end
            sectorIdx(mask) = si;
        end
        counts = accumarray(double(sectorIdx),1,[16 1]);
        ptColors = cmap16_FT(double(sectorIdx),:);
        set(hScat,'XData',xA,'YData',yA,'CData',ptColors);
        totalN = numel(xA);
    else
        set(hScat,'XData',nan,'YData',nan);
        counts = zeros(16,1);
        totalN = 0;
    end
    % 生成每组的文本和颜色
    axisIdx = [1,3,5,7,9,11,13,15];
    gapIdx  = [2,4,6,8,10,12,14,16];
    oridIdx = [3,7,11,15];

    axisSum = sum(counts(axisIdx));
    gapSum  = sum(counts(gapIdx));
    oridSum = sum(counts(oridIdx));

    % 文本内容
    infoStr1 = sprintf('Time: %d ms   Total Fix: %d', round(tNow), totalN);

    % 每扇区文本（分两行）
    line1 = '';
    for i = 1:8
        line1 = [line1, sprintf('S%02d:%d  ', i, counts(i))];
    end
    line2 = '';
    for i = 9:16
        line2 = [line2, sprintf('S%02d:%d  ', i, counts(i))];
    end

    axisText = sprintf('Axis:%d', axisSum);
    gapText  = sprintf('Gap:%d', gapSum);
    oridText = sprintf('Orid:%d', oridSum);

    % 清除旧文本
    delete(findall(ax, 'Tag', 'fixText'));

    % 绘制文本（分行分色）
    y0 = 40; dy = 22;
    x0 = 15;

    % 第一行：时间和总fix
    text(x0, y0, infoStr1, 'Parent', ax, 'Color', [0 0 0], ...
        'FontSize', 13, 'FontWeight', 'bold', 'Tag', 'fixText', ...
        'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');

    % 第二行：S01-S08
    x1 = x0; y1 = y0 + dy;
    xstep = 90;
    for i = 1:8
        str = sprintf('S%02d:%d', i, counts(i));
        text(x1 + (i-1)*xstep, y1, str, 'Parent', ax, ...
            'Color', cmap16_FT(i,:), 'FontSize', 12, 'FontWeight', 'bold', ...
            'Tag', 'fixText', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');
    end

    % 第三行：S09-S16
    y2 = y1 + dy;
    for i = 9:16
        str = sprintf('S%02d:%d', i, counts(i));
        text(x1 + (i-9)*xstep, y2, str, 'Parent', ax, ...
            'Color', cmap16_FT(i,:), 'FontSize', 12, 'FontWeight', 'bold', ...
            'Tag', 'fixText', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');
    end

    % 第四行：Axis/Gap/Orid
    y3 = y2 + dy;
    text(x0, y3, axisText, 'Parent', ax, 'Color', cmap16_FT(1,:), ...
        'FontSize', 13, 'FontWeight', 'bold', 'Tag', 'fixText', ...
        'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');
    text(x0+120, y3, gapText, 'Parent', ax, 'Color', cmap16_FT(2,:), ...
        'FontSize', 13, 'FontWeight', 'bold', 'Tag', 'fixText', ...
        'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');
    text(x0+240, y3, oridText, 'Parent', ax, 'Color', cmap16_FT(3,:), ...
        'FontSize', 13, 'FontWeight', 'bold', 'Tag', 'fixText', ...
        'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');

    drawnow;
    writeVideo(vw_fp, getframe(fig));
    end

close(vw_fp);
close(fig);
fprintf('视频已保存: %s\n', vidName);
end

