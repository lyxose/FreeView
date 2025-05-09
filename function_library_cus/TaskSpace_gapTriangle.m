function [Eccent, Orient, x, y] = TaskSpace_gapTriangle(Side_region, Width_region, Gap_region, trialNum, noiseP, scWidth, scHeight, bgWidth, rot_ang, tgSeed)
% right down as positive! 
% Side_region:  horizontal range of the equilateral triangle base side
% Width_region: vertical range of base side line width
if nargin < 7
    tgSeed = randi(10000);
end
ut = UT(54, scWidth, 60, false); % assumption just for the target location generation
canvas = zeros(scHeight,scWidth);
[X, Y] = meshgrid(1:scWidth, 1:scHeight);

bg_circle = (X - scWidth/2).^2 + (Y - scHeight/2).^2 <= (ut.deg2pix(bgWidth)/2)^2;

xc = ut.deg2pix(Side_region);%[-7 9]);
yc = ut.deg2pix(Width_region);%[1 1.1]);
xcg = ut.deg2pix(Gap_region);%[-1 3]);

xshift = sum(xc)/2;
yshift = -diff(xc)/2/sqrt(3) + sum(yc)/2;

xc = xc - xshift;
yc = yc - yshift;
xcg = xcg - xshift;


vertices_x = [xc(1) xc(1) xc(2) xc(2)]+scWidth/2;
vertices_y = [yc(2) yc(1) yc(1) yc(2)]+scHeight/2;

ROI_matrix = poly2mask(vertices_x, vertices_y, scHeight, scWidth);

% the gap region
% yc = yc;
vertices_x = [xcg(1) xcg(1) xcg(2) xcg(2)]+scWidth/2;
vertices_y = [yc(2) yc(1) yc(1) yc(2)]+scHeight/2;

gap_matrix = poly2mask(vertices_x, vertices_y, scHeight, scWidth);
ROI_matrix = ROI_matrix & ~gap_matrix;

canvas(ROI_matrix) = 1;

rotated_canvas1 = imrotate(canvas, 240, 'bilinear', 'crop');
rotated_canvas2 = imrotate(canvas, 120, 'bilinear', 'crop');
canvas(bg_circle) = canvas(bg_circle) + rotated_canvas1(bg_circle) + rotated_canvas2(bg_circle);
canvas(canvas>1) = 1;

% round corner
r = diff(yc)/2;
loc1 = [scWidth/2, scHeight/2] + [0,-diff(xc)/sqrt(3)];
loc2 = [scWidth/2, scHeight/2] + [-diff(xc)/2,diff(xc)/2/sqrt(3)];
loc3 = [scWidth/2, scHeight/2] + [diff(xc)/2, diff(xc)/2/sqrt(3)];
circle1 = (X - loc1(1)).^2 + (Y - loc1(2)).^2 <= r^2;
circle2 = (X - loc2(1)).^2 + (Y - loc2(2)).^2 <= r^2;
circle3 = (X - loc3(1)).^2 + (Y - loc3(2)).^2 <= r^2;
canvas(circle1) = 1;
canvas(circle2) = 1;
canvas(circle3) = 1;

canvas = imtranslate(canvas, [xshift, yshift], 'FillValues', 0);

if any(canvas(~bg_circle))
    warning('ROI exceed the background range, please check carefully!!')
end
% probability
ROI_matrix = canvas>0;
canvas(ROI_matrix) = (1-noiseP) / sum(ROI_matrix(:));

% background noise 
bgIdx = bg_circle & ~ROI_matrix;
canvas(bgIdx) = noiseP / sum(bgIdx(:));

% rotate
rotated_canvas = imrotate(canvas, rot_ang, 'bilinear', 'crop');
canvas(bg_circle) = rotated_canvas(bg_circle);


full_pdf = canvas ./ sum(canvas(:));            % 
% figure()
% imagesc(full_pdf)
% disp(max(canvas(:)))
% sampling
prob_vector = full_pdf(:);
original_rng = rng;    % original seed
rng(tgSeed);           % use customized seed for each subject
indices = randsample(1:numel(prob_vector), trialNum, true, prob_vector);
rng(original_rng);     % 

[y, x] = ind2sub(size(full_pdf), indices);
PolarCoor=ut.Rect2Pol(transpose([x-scWidth/2; -y+scHeight/2]));
Orient = PolarCoor(:,2);
Eccent = PolarCoor(:,1);

