function [Eccent, Orient, x, y] = TaskSpace_gapRect(Hori_region, Verti_region, Gap_region, trialNum, noiseP, scWidth, scHeight, bgWidth, rot_ang, tgSeed)
if nargin < 7
    tgSeed = randi(10000);
end
ut = UT(54, scWidth, 60, true); % assumption just for the target location generation
canvas = zeros(scHeight,scWidth);
[X, Y] = meshgrid(1:scWidth, 1:scHeight);

bg_circle = (X - scWidth/2).^2 + (Y - scHeight/2).^2 <= (ut.deg2pix(bgWidth)/2)^2;

xc = ut.deg2pix(Hori_region);%[-7 9]);
yc = ut.deg2pix(Verti_region);%[1 1.1]);
vertices_x = [xc(1) xc(1) xc(2) xc(2)]+scWidth/2;
vertices_y = [yc(2) yc(1) yc(1) yc(2)]+scHeight/2;

ROI_matrix = poly2mask(vertices_x, vertices_y, scHeight, scWidth);

% the gap region
xc = ut.deg2pix(Gap_region);%[-1 3]);
% yc = yc;
vertices_x = [xc(1) xc(1) xc(2) xc(2)]+scWidth/2;
vertices_y = [yc(2) yc(1) yc(1) yc(2)]+scHeight/2;

gap_matrix = poly2mask(vertices_x, vertices_y, scHeight, scWidth);
ROI_matrix = ROI_matrix & ~gap_matrix;
if any(ROI_matrix(~bg_circle))
    warning('ROI exceed the background range, please check carefully!!')
end
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

