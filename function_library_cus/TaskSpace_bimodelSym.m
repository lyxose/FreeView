function [Eccent, Orient, x, y] = TaskSpace_bimodelSym(trialNum, scWidth, scHeight, bgWidth, tgSeed, visualiztion)
%TASKSPACE_BIMODELSYM Generate bimodal probability distribution and sample points
%   Outputs: Eccent - Eccentricity, Orient - Orientation angle
%            x,y - Sampled coordinates in pixel space

%% Parameter initialization
if nargin < 6
    visualiztion=false;
end
if nargin < 5
    tgSeed = randi(10000); % Generate random seed if not provided
end

ut = UT(54, scWidth, 60, true); % Coordinate system converter
canvas = zeros(scHeight, scWidth);
[X, Y] = meshgrid(1:scWidth, 1:scHeight);

%% Create circular background region
bg_circle = (X - scWidth/2).^2 + (Y - scHeight/2).^2 <= (ut.deg2pix(bgWidth)/2)^2;

%% Calculate 45° diagonal centers
bgCenterX = scWidth/2;
bgCenterY = scHeight/2;
radius = ut.deg2pix(bgWidth)/2 * 0.6   ; % 60% of max radius

% Convert polar to Cartesian coordinates (135° and 225°)
[dx1, dy1] = pol2cart(deg2rad(135), radius);
[dx2, dy2] = pol2cart(deg2rad(225), radius);

% Adjust for screen coordinates (Y-axis inverted)
mu_x1 = bgCenterX + dx1;
mu_y1 = bgCenterY - dy1; 
mu_x2 = bgCenterX + dx2;
mu_y2 = bgCenterY - dy2;

%% Generate bimodal Gaussian distribution
sigma = radius * 0.15; % 15% of radius as sigma
gauss2D = @(x,y,mu_x,mu_y) exp(-((x-mu_x).^2 + (y-mu_y).^2)/(2*sigma^2));

gauss1 = gauss2D(X, Y, mu_x1, mu_y1);
gauss2 = gauss2D(X, Y, mu_x2, mu_y2);
combined_gauss = gauss1 + gauss2;
canvas(bg_circle) = combined_gauss(bg_circle);

%% Create gray semicircle in right half
% theta = linspace(pi/2, 3*pi/2, 100);
% r = linspace(0, ut.deg2pix(bgWidth)/2, 50);
% [R, Theta] = meshgrid(r, theta);
% 
% X_semi = R .* cos(Theta) + bgCenterX;
% Y_semi = R .* sin(Theta) + bgCenterY;
% semi_mask = poly2mask(X_semi(:), Y_semi(:), scHeight, scWidth);
leftSpace = canvas(X<=scWidth/2 & bg_circle);
leftWeight = sum(leftSpace,"all")/length(leftSpace);
canvas(X>scWidth/2 & bg_circle) = leftWeight;


%% Normalize probability distribution and sample points
prob_map = canvas ./ sum(canvas(:));
rng(tgSeed);
indices = randsample(numel(prob_map), trialNum, true, prob_map(:));
[y, x] = ind2sub(size(prob_map), indices);

%% visualiztion
if visualiztion
% Visualize probability density distribution
figure('Name','Probability Density Map');
imagesc(canvas);
colormap(hot);
colorbar;
title('Normalized Probability Density Distribution');
axis equal tight;

% Visualize sampling distribution with 1:1 aspect
figure('Name','Sampling Distribution');
scatter(x, y, 20, 'filled', 'MarkerFaceAlpha',0.4);
hold on;

% Plot 45° reference line
plot([bgCenterX-radius bgCenterX+radius],...
     [bgCenterY+radius bgCenterY-radius],...
     '--', 'Color',[0.7 0.7 0.7], 'LineWidth',1.5);

% Mark Gaussian centers
plot(mu_x1, mu_y1, 'pentagram',...
    'MarkerSize',14,...
    'MarkerFaceColor','r',...
    'MarkerEdgeColor','k');
plot(mu_x2, mu_y2, 'pentagram',...
    'MarkerSize',14,...
    'MarkerFaceColor','r',...
    'MarkerEdgeColor','k');

% Axis formatting
daspect([1 1 1]); % Enforce 1:1 aspect ratio
xlim([1 scWidth]);
ylim([1 scHeight]);
set(gca, 'YDir','reverse',...
    'XAxisLocation','top',...
    'FontSize',11,...
    'GridColor',[0.9 0.9 0.9]);
title(sprintf('Sampled Points (n=%d)', trialNum));
grid on;
end
%% Coordinate conversion
PolarCoor = ut.Rect2Pol([x - scWidth/2, -(y - scHeight/2)]);
Eccent = PolarCoor(:,1);
Orient = mod(PolarCoor(:,2), 360);
end