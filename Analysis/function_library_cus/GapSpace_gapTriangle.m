function sumCounts = GapSpace_gapTriangle(fixPos, Side_region, Width_region, Gap_region, scWidth, scHeight, bgWidth, rot_ang)
% right down as positive! 
% Side_region:  horizontal range of the equilateral triangle base side
% Width_region: vertical range of base side line width
% if nargin < 7
%     tgSeed = randi(10000);
% end

ut = UT(54, scWidth, 60, false); % assumption just for the target location generation
canvas = zeros(scHeight,scWidth);
[X, Y] = meshgrid(1:scWidth, 1:scHeight);

bg_circle = (X - scWidth/2).^2 + (Y - scHeight/2).^2 <= (ut.deg2pix(bgWidth)/2)^2;

xc = ut.deg2pix(Side_region);%[-7 9]);
yc = ut.deg2pix(Width_region);%[1 1.1]);
xcg = ut.deg2pix(Gap_region);%[-1 3]);

xshift = sum(xc)/2;
yshift = -diff(xc)/2/sqrt(3) + sum(yc)/2;

xcg = xcg - xshift;


% vertices_x = [xc(1) xc(1) xc(2) xc(2)]+scWidth/2;
% vertices_y = [yc(2) yc(1) yc(1) yc(2)]+scHeight/2;
% 
% ROI_matrix = poly2mask(vertices_x, vertices_y, scHeight, scWidth);

% the gap region
% yc = yc;
vertices_x = [xcg(1) xcg(1) xcg(2) xcg(2)]+scWidth/2;
vertices_y = [yc(2) yc(1) yc(1) yc(2)]+scHeight/2;

gap_matrix = poly2mask(vertices_x, vertices_y, scHeight, scWidth);
% ROI_matrix = ROI_matrix & ~gap_matrix;

canvas(gap_matrix) = 1;

rotated_canvas1 = imrotate(canvas, 240, 'bilinear', 'crop');
rotated_canvas2 = imrotate(canvas, 120, 'bilinear', 'crop');
canvas(bg_circle) = canvas(bg_circle) + rotated_canvas1(bg_circle) + rotated_canvas2(bg_circle);
canvas(canvas>1) = 1;

% round corner
% r = diff(yc)/2;
% loc1 = [scWidth/2, scHeight/2] + [0,-diff(xc)/sqrt(3)];
% loc2 = [scWidth/2, scHeight/2] + [-diff(xc)/2,diff(xc)/2/sqrt(3)];
% loc3 = [scWidth/2, scHeight/2] + [diff(xc)/2, diff(xc)/2/sqrt(3)];
% circle1 = (X - loc1(1)).^2 + (Y - loc1(2)).^2 <= r^2;
% circle2 = (X - loc2(1)).^2 + (Y - loc2(2)).^2 <= r^2;
% circle3 = (X - loc3(1)).^2 + (Y - loc3(2)).^2 <= r^2;
% canvas(circle1) = 1;
% canvas(circle2) = 1;
% canvas(circle3) = 1;

% region of each angle elements 
% locs = {loc1, loc2, loc3};
% if 0<min(xcg) || 0>max(xcg)
%     warning('Current method of deviding angle elements is not suitable for the gap which is not crossing the middle of each side，please check!');
% end
% r_angle = diff(xc)/2;
% elementMasks={[],[],[]};
% for i =1:3
%     loci = locs{i};
%     elementMasks{i} = canvas & (X - loci(1)).^2 + (Y - loci(2)).^2 <= r_angle^2;
%     elementMasks{i} = imtranslate(elementMasks{i}, [xshift, yshift], 'FillValues', 0);
%     elementMasks{i} = imrotate(elementMasks{i}, rot_ang, 'bilinear', 'crop');
% end
canvas = imtranslate(canvas, [xshift, yshift], 'FillValues', 0);

if any(canvas(~bg_circle))
    warning('ROI exceed the background range, please check carefully!!')
end

% rotate
rotated_canvas = imrotate(canvas, rot_ang, 'bilinear', 'crop');
canvas(bg_circle) = rotated_canvas(bg_circle);

edges_x = 0:scWidth;% linspace(1, img_width,  round(img_width/binSize));  % 网格的X边界
edges_y = 0:scHeight;% linspace(1, img_height, round(img_height/binSize)); % 网格的Y边界
counts = hist3(fixPos, 'Edges', {edges_x(1:end-1), edges_y(1:end-1)});
sumCounts = sum(canvas.*transpose(counts),"all");
% full_pdf = canvas ./ sum(canvas(:));            % 
% figure()
% imagesc(full_pdf)
% disp(max(canvas(:)))
% sampling
% prob_vector = full_pdf(:);
% original_rng = rng;    % original seed
% rng(tgSeed);           % use customized seed for each subject
% indices = randsample(1:numel(prob_vector), trialNum, true, prob_vector);
% rng(original_rng);     % 
% 
% [y, x] = ind2sub(size(full_pdf), indices);
% % tag all samples
% mask1 = elementMasks{1}(indices);
% mask2 = elementMasks{2}(indices);
% mask3 = elementMasks{3}(indices);
% elementClusterTags = mask1 * 1 + mask2.*(~mask1)*2 + mask3.*(~mask1 & ~mask2)*3;
% 
% PolarCoor=ut.Rect2Pol(transpose([x-scWidth/2; -y+scHeight/2]));
% Orient = PolarCoor(:,2);
% Eccent = PolarCoor(:,1);

