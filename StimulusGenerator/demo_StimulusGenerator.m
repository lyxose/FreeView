%% 刺激生成器Demo - Stimulus Generator Demo
% 本脚本演示如何使用StimulusGenerator工具包生成自定义视觉刺激
% This script demonstrates how to use the StimulusGenerator toolkit to generate custom visual stimuli
%
% 作者 / Author: Yuxin Lu, IPCAS
% 日期 / Date: 2025.1.7

clear; close all; clc;

%% ========== 步骤1: 配置屏幕参数 / Step 1: Configure Screen Parameters ==========

% 屏幕物理参数 / Screen physical parameters
screenWidth_cm = 53.5;      % 屏幕宽度 (厘米) / Screen width (cm)
screenWidth_px = 1920;      % 屏幕宽度 (像素) / Screen width (pixels)
screenHeight_px = 1080;     % 屏幕高度 (像素) / Screen height (pixels)
viewingDistance_cm = 68;    % 观察距离 (厘米) / Viewing distance (cm)

% 创建单位转换器 / Create unit transformer
ut = UT(screenWidth_cm, screenWidth_px, viewingDistance_cm);

fprintf('===== 屏幕配置 / Screen Configuration =====\n');
fprintf('屏幕尺寸 / Screen size: %d x %d 像素 / pixels\n', screenWidth_px, screenHeight_px);
fprintf('屏幕宽度 / Screen width: %.1f cm\n', screenWidth_cm);
fprintf('观察距离 / Viewing distance: %.1f cm\n', viewingDistance_cm);
fprintf('每厘米像素数 / Pixels per cm: %.2f\n', ut.ppcm);
fprintf('1度视角 = %.2f 像素 / 1 degree = %.2f pixels\n\n', ut.deg2pix(1));

%% ========== 步骤2: 设置刺激参数 / Step 2: Set Stimulus Parameters ==========

% 背景参数 / Background parameters
bgWidth_deg = 15;           % 背景圆形区域直径 (度) / Background circular area diameter (degrees)
bgContrast = 0.2;           % 背景对比度 (0-1) / Background contrast (0-1)

% Gabor目标参数 / Gabor target parameters
tgContrast = 0.3;           % 目标对比度 (0-1) / Target contrast (0-1)
GaborSF = 6;                % 空间频率 (cycles/degree) / Spatial frequency (cycles/degree)
GaborCyc = 2;               % Gabor周期数 / Number of Gabor cycles
GaborWidth_deg = GaborCyc / GaborSF;  % Gabor宽度 (度) / Gabor width (degrees)
GaborOrient = -45;          % Gabor方向 (度) / Gabor orientation (degrees)

% 目标位置 (视角度数) / Target position (degrees)
tgEccent_deg = 4;           % 偏心度 (度) / Eccentricity (degrees)
tgAngle_deg = 45;           % 角度 (度) / Angle (degrees)

% 随机种子 / Random seed
seed = 12345;

fprintf('===== 刺激参数 / Stimulus Parameters =====\n');
fprintf('背景直径 / Background diameter: %.1f 度 / degrees\n', bgWidth_deg);
fprintf('背景对比度 / Background contrast: %.2f\n', bgContrast);
fprintf('目标对比度 / Target contrast: %.2f\n', tgContrast);
fprintf('Gabor空间频率 / Gabor SF: %.1f cycles/degree\n', GaborSF);
fprintf('Gabor宽度 / Gabor width: %.2f 度 / degrees\n', GaborWidth_deg);
fprintf('Gabor方向 / Gabor orientation: %.0f 度 / degrees\n', GaborOrient);
fprintf('目标位置 / Target position: 偏心度 / eccentricity %.1f°, 角度 / angle %.0f°\n\n', tgEccent_deg, tgAngle_deg);

%% ========== 步骤3: 计算目标位置 / Step 3: Calculate Target Position ==========

% 将极坐标转换为直角坐标 (度数) / Convert polar to rectangular coordinates (degrees)
tgX_deg = tgEccent_deg * cosd(tgAngle_deg);
tgY_deg = tgEccent_deg * sind(tgAngle_deg);

% 转换为像素坐标 (图像坐标系: 左上角为原点，右下为正方向)
% Convert to pixel coordinates (image coordinate system: upper-left origin, right-down positive)
screenCenter_px = [screenWidth_px/2, screenHeight_px/2];
tgCenter_px = ut.deg2pix([tgX_deg, tgY_deg]);
tgCenter_px = tgCenter_px .* [1, -1] + screenCenter_px;  % 转换坐标系 / Convert coordinate system

fprintf('===== 坐标计算 / Coordinate Calculation =====\n');
fprintf('目标位置 (视角) / Target position (degrees): (%.2f, %.2f)\n', tgX_deg, tgY_deg);
fprintf('目标位置 (像素) / Target position (pixels): (%.0f, %.0f)\n', tgCenter_px(1), tgCenter_px(2));
fprintf('屏幕中心 (像素) / Screen center (pixels): (%.0f, %.0f)\n\n', screenCenter_px(1), screenCenter_px(2));

%% ========== 步骤4: 生成刺激 / Step 4: Generate Stimulus ==========

% 生成实验用刺激
stimulus = genStim([screenWidth_px, screenHeight_px], ut, bgContrast, tgContrast, ...
                   tgCenter_px, GaborSF, GaborWidth_deg, GaborOrient, bgWidth_deg, seed);

figure();
imshow(stimulus);


% 单独生成背景图像
background = tPinkNoise(screenWidth_px, seed, bgContrast);
background = background(1:screenHeight_px, 1:screenWidth_px);
background_circular = winOverlap(zeros([screenHeight_px, screenWidth_px]) + 0.5, ...
                                  background, ut.deg2pix(bgWidth_deg), ...
                                  screenCenter_px, 'hard');
figure();
imshow(background_circular);

% 单独生成圆形小gabor
lambda = ut.deg2pix(1/GaborSF);
Texture = grating([screenHeight_px, screenWidth_px], tgCenter_px, ...
                    1/lambda, GaborOrient, tgContrast);
gabor_circular = winOverlap(zeros([screenHeight_px, screenWidth_px]) + 0.5, ...
                              Texture, ut.deg2pix(GaborWidth_deg), ...
                              tgCenter_px, 'cos');
figure();
imshow(gabor_circular);