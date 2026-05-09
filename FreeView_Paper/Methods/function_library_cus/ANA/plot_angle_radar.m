function plot_angle_radar(data, centers, foldPeriod, binSize, CardColor, GapColor, ObliColor, normMode, varargin)
% PLOT_ANGLE_RADAR Draw scan curve as a radar/polar chart
% data: Nsubj x nbins or 1 x nbins
% centers: bin centers (degrees)
% foldPeriod: 360/90/45
% binSize: bin width (degrees)
% Colors: as in plot_angle_curve
% normMode: label for radius

if isvector(data)
    data = reshape(data, 1, []);
end
nsbj = size(data,1);
groupMean = mean(data,1);
if nsbj > 1
    groupSE = std(data,0,1)/sqrt(nsbj);
else
    groupSE = zeros(size(groupMean));
end

% close the loop for plotting
theta = deg2rad(centers);
theta_closed = [theta, theta(1)];
mean_closed = [groupMean, groupMean(1)];
se_hi = groupMean + groupSE;
se_lo = groupMean - groupSE;
se_hi_closed = [se_hi, se_hi(1)];
se_lo_closed = [se_lo, se_lo(1)];

figure('Color','w'); hold on;
% draw CI band when multiple subjects
if nsbj > 1
    % build patch polygon in Cartesian coords
    [x_outer, y_outer] = pol2cart(theta_closed, se_hi_closed);
    [x_inner, y_inner] = pol2cart(fliplr(theta_closed), fliplr(se_lo_closed));
    x_poly = [x_outer, x_inner];
    y_poly = [y_outer, y_inner];
    patch(x_poly, y_poly, [0.8 0.8 0.8], 'EdgeColor','none', 'FaceAlpha', 0.3);
end
% plot mean
[xm, ym] = pol2cart(theta_closed, mean_closed);
plot(xm, ym, '-k', 'LineWidth', 1.6);
% markers
[xm_m, ym_m] = pol2cart(theta, groupMean);
plot(xm_m, ym_m, 'o', 'MarkerFaceColor',[0.2 0.2 0.2], 'MarkerEdgeColor','none');

% formatting: draw radial grid lines and angle ticks
maxr = max(mean_closed) * 1.15;
if maxr == 0, maxr = 1; end
% draw concentric circles
nrings = 4;
rs = linspace(0, maxr, nrings+1);
for r = rs
    th = linspace(0, 2*pi, 200);
    [xc, yc] = pol2cart(th, r);
    plot(xc, yc, ':', 'Color', [0.6 0.6 0.6]);
end
% angle spokes and labels every 45 deg when 360
if foldPeriod == 360
    angs = 0:45:315;
else
    angs = centers;
end
for a = angs
    [xs, ys] = pol2cart(deg2rad(a), [0, maxr]);
    plot(xs, ys, '-', 'Color', [0.85 0.85 0.85]);
    [xl, yl] = pol2cart(deg2rad(a), maxr*1.05);
    txt = sprintf('%.0f°', mod(a,360));
    text(xl, yl, txt, 'HorizontalAlignment','center', 'VerticalAlignment','middle');
end
axis equal off;
% legend/title
title(sprintf('Angle Radar (n=%d)', nsbj));
if nargin>=8 && ischar(normMode)
    ylabel_str = normMode;
else
    ylabel_str = '';
end
% draw a small axis indicator
text(-maxr*0.02, -maxr*0.02, ylabel_str, 'HorizontalAlignment','left');
hold off;
end
