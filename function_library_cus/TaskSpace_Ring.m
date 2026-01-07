function [Eccent, Orient, x, y, elementClusterTags] = TaskSpace_Ring(trialNum, scWidth, scHeight, R_max, R_min, tgSeed)
%TASKSPACE_RING Sample targets uniformly in angle and within an annulus.
%   trialNum: number of samples
%   scWidth, scHeight: screen resolution in pixels
%   R_max, R_min: outer/inner radius in deg (defaults 7 / 2)
%   tgSeed: optional RNG seed for reproducibility

if nargin < 6 || isempty(tgSeed)
    tgSeed = randi(10000);
end
if nargin < 5 || isempty(R_min)
    R_min = 2;
end
if nargin < 4 || isempty(R_max)
    R_max = 7;
end

ut = UT(54, scWidth, 60, true); % deg/pix converter
original_rng = rng;
rng(tgSeed);
Orient = rand(trialNum, 1) * 360;
Eccent = sqrt(R_min.^2 + (R_max.^2 - R_min.^2) .* rand(trialNum, 1));
rng(original_rng);

XY = ut.Pol2Rect([Eccent, Orient]);
x = XY(:, 1) + scWidth/2;
y = -XY(:, 2) + scHeight/2;

if nargout > 4
    elementClusterTags = ones(1,trialNum);
end
end
