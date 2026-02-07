# 刺激生成器工具包 / Stimulus Generator Toolkit

## 📋 简介 / Introduction

这是一个用于生成视觉心理物理学实验刺激的MATLAB工具包。可以在1/f粉噪声背景上生成Gabor目标刺激。

This is a MATLAB toolkit for generating visual psychophysics experimental stimuli. It generates Gabor target stimuli on 1/f pink noise backgrounds.

**作者 / Author:** Yuxin Lu, IPCAS  
**邮箱 / Email:** luyx@psych.ac.cn  
**日期 / Date:** 2025.1.7

---

## 📦 工具包内容 / Package Contents

### 核心函数 / Core Functions

1. **`genStim.m`** - 主刺激生成函数 / Main stimulus generation function
2. **`UT.m`** - 单位转换类 / Unit transformation class
3. **`grating.m`** - 正弦光栅生成 / Sinusoidal grating generation
4. **`tPinkNoise.m`** - 1/f粉噪声生成 / 1/f pink noise generation
5. **`winOverlap.m`** - 窗口混合函数 / Window blending function

### 示例文件 / Example Files

- **`demo_StimulusGenerator.m`** - 完整使用示例 / Complete usage example
- **`README.md`** - 本文档 / This document

---

## 🚀 快速开始 / Quick Start

### 1. 安装 / Installation

将 `StimulusGenerator` 文件夹复制到你的MATLAB路径中，或在脚本中添加：

Copy the `StimulusGenerator` folder to your MATLAB path, or add in your script:

```matlab
addpath('path/to/StimulusGenerator');
```

### 2. 基本使用 / Basic Usage

```matlab
% 步骤1: 配置屏幕参数
screenWidth_cm = 53.5;
screenWidth_px = 1920;
screenHeight_px = 1080;
viewingDistance_cm = 68;
ut = UT(screenWidth_cm, screenWidth_px, viewingDistance_cm);

% 步骤2: 设置刺激参数
bgWidth_deg = 15;
bgContrast = 0.2;
tgContrast = 0.3;
GaborSF = 6;
GaborCyc = 2;
GaborWidth_deg = GaborCyc / GaborSF;
GaborOrient = -45;
tgEccent_deg = 4;   % 偏心度
tgAngle_deg = 45;   % 角度
seed = 12345;

% 步骤3: 计算目标位置
tgX_deg = tgEccent_deg * cosd(tgAngle_deg);
tgY_deg = tgEccent_deg * sind(tgAngle_deg);
screenCenter_px = [screenWidth_px/2, screenHeight_px/2];
tgCenter_px = ut.deg2pix([tgX_deg, tgY_deg]);
tgCenter_px = tgCenter_px .* [1, -1] + screenCenter_px;

% 步骤4: 生成刺激
stimulus = genStim([screenWidth_px, screenHeight_px], ut, bgContrast, tgContrast, ...
                   tgCenter_px, GaborSF, GaborWidth_deg, GaborOrient, bgWidth_deg, seed);

imshow(stimulus, []); colorbar;
```

### 3. 运行Demo / Run Demo

```matlab
cd StimulusGenerator
demo_StimulusGenerator
```

---

## 📚 详细文档 / Detailed Documentation

### UT 类 / UT Class

单位转换工具类，用于视角度数、像素、厘米之间的转换。

Unit conversion utility class for converting between visual degrees, pixels, and centimeters.

#### 构造函数 / Constructor

```matlab
ut = UT(screenWidth_cm, screenWidth_px, viewingDistance_cm, rndPix)
```

**参数 / Parameters:**
- `screenWidth_cm`: 屏幕宽度(厘米) / Screen width (cm)
- `screenWidth_px`: 水平像素数 / Horizontal pixels
- `viewingDistance_cm`: 观察距离(厘米) / Viewing distance (cm)
- `rndPix`: 是否四舍五入像素值(默认true) / Whether to round pixels (default true)

#### 主要方法 / Main Methods

| 方法 / Method | 功能 / Function | 示例 / Example |
|--------------|----------------|---------------|
| `deg2pix(deg)` | 度数→像素 / Degrees to pixels | `ut.deg2pix(5)` |
| `pix2deg(pix)` | 像素→度数 / Pixels to degrees | `ut.pix2deg(100)` |
| `cm2pix(cm)` | 厘米→像素 / cm to pixels | `ut.cm2pix(10)` |
| `pix2cm(pix)` | 像素→厘米 / Pixels to cm | `ut.pix2cm(100)` |
| `Pol2Rect(polar)` | 极坐标→直角坐标 / Polar to Rect | `ut.Pol2Rect([5, 45])` |
| `Rect2Pol(rect)` | 直角坐标→极坐标 / Rect to Polar | `ut.Rect2Pol([10, 10])` |

---

### genStim 函数 / genStim Function

生成完整的视觉刺激（Gabor目标 + 1/f噪声背景）。

Generate complete visual stimulus (Gabor target + 1/f noise background).

#### 语法 / Syntax

```matlab
stimulus = genStim(winRect, ut, bgContrast, tgContrast, tgCenter, ...
                   GaborSF, GaborWidth, GaborOrient, bgWidth, seed)
```

#### 参数 / Parameters

| 参数 / Parameter | 类型 / Type | 说明 / Description |
|-----------------|------------|-------------------|
| `winRect` | 数组 / Array | 窗口尺寸 `[width, height]` 或 `[x0, y0, width, height]` |
| `ut` | UT对象 / UT object | 单位转换器实例 / Unit transformer instance |
| `bgContrast` | 标量 / Scalar | 背景对比度 (0-1) / Background contrast |
| `tgContrast` | 标量 / Scalar | 目标对比度 (0-1) / Target contrast |
| `tgCenter` | 数组 / Array | 目标中心像素坐标 `[x, y]` / Target center in pixels |
| `GaborSF` | 标量 / Scalar | 空间频率 (cpd) / Spatial frequency (cycles/degree) |
| `GaborWidth` | 标量 / Scalar | Gabor宽度(度, FWHM) / Gabor width (degrees, FWHM) |
| `GaborOrient` | 标量 / Scalar | Gabor方向(度) / Gabor orientation (degrees) |
| `bgWidth` | 标量 / Scalar | 背景圆形直径(度) / Background diameter (degrees) |
| `seed` | 整数 / Integer | 随机种子 / Random seed |

#### 返回值 / Returns

- `stimulus`: 刺激图像矩阵，范围 `[0, 1]` / Stimulus matrix, range `[0, 1]`

---

### grating 函数 / grating Function

生成正弦光栅图像。

Generate sinusoidal grating image.

#### 语法 / Syntax

```matlab
grating = grating(sizePix, centerLoc, freqPix, orientDeg, contrast, phaseRad, show)
```

#### 参数 / Parameters

- `sizePix`: 图像尺寸 `[height, width]` / Image size
- `centerLoc`: 中心位置 `[x, y]`像素 / Center location in pixels
- `freqPix`: 空间频率 (cycles/pixel) / Spatial frequency
- `orientDeg`: 方向(度) / Orientation (degrees)
- `contrast`: 对比度 (0-1) / Contrast
- `phaseRad`: 相位(弧度, 默认0) / Phase (radians, default 0)
- `show`: 是否显示(默认false) / Whether to display (default false)

---

### tPinkNoise 函数 / tPinkNoise Function

生成2D 1/f粉噪声纹理。

Generate 2D 1/f pink noise texture.

#### 语法 / Syntax

```matlab
noise = tPinkNoise(N, seed, contrast, show)
```

#### 参数 / Parameters

- `N`: 图像尺寸 (N×N) / Image size (N×N)
- `seed`: 随机种子 / Random seed
- `contrast`: 对比度(0-1, 默认1) / Contrast (0-1, default 1)
- `show`: 是否显示(默认false) / Whether to display (default false)

---

### winOverlap 函数 / winOverlap Function

使用窗口混合两个图像。

Blend two images using a window.

#### 语法 / Syntax

```matlab
img = winOverlap(background, source, widthPix, centerLoc, windowType, show)
```

#### 参数 / Parameters

- `background`: 背景图像 / Background image
- `source`: 源图像 / Source image
- `widthPix`: 窗口宽度(像素, FWHM) / Window width (pixels, FWHM)
- `centerLoc`: 窗口中心 `[x, y]` / Window center
- `windowType`: 窗口类型: `'cos'`, `'Gaussian'`, `'linear'`, `'hard'` (默认`'cos'`)
- `show`: 是否显示(默认false) / Whether to display (default false)

---

## 💡 使用技巧 / Usage Tips

### 1. 坐标系统说明 / Coordinate System

本工具包使用**图像坐标系**：
- 原点 (0, 0) 在**左上角** / Origin at **upper-left**
- X轴向**右**为正 / X-axis: **right** is positive
- Y轴向**下**为正 / Y-axis: **down** is positive

视觉坐标系转换：
```matlab
% 视觉坐标 (右上为正) → 图像坐标 (右下为正)
% Visual coords (right-up positive) → Image coords (right-down positive)
screenCenter = [width/2, height/2];
tgCenter_image = [tgX, -tgY] + screenCenter;
```

### 2. 参数选择建议 / Parameter Selection Recommendations

| 参数 / Parameter | 推荐范围 / Recommended Range | 说明 / Notes |
|-----------------|------------------------------|-------------|
| 背景对比度 / Background contrast | 0.1 - 0.3 | 过高会掩蔽目标 / Too high masks target |
| 目标对比度 / Target contrast | 0.2 - 0.5 | 根据实验需求调整 / Adjust per experiment |
| Gabor空间频率 / SF | 4 - 8 cpd | 标准foveal范围 / Standard foveal range |
| Gabor宽度 / Width | 1 - 3度 / degrees | 2周期为常用值 / 2 cycles is common |
| 背景直径 / Background diameter | 10 - 20度 / degrees | 覆盖中心视野 / Cover central vision |

### 3. 性能优化 / Performance Optimization

```matlab
% 预先计算多个刺激时，重复使用UT对象
% Reuse UT object when pre-computing multiple stimuli
ut = UT(screenWidth_cm, screenWidth_px, viewingDistance_cm);

% 批量生成
% Batch generation
stimuli = cell(1, nTrials);
for i = 1:nTrials
    stimuli{i} = genStim(...);  % 复用ut / Reuse ut
end
```

### 4. 常见问题排查 / Troubleshooting

**问题: 生成的刺激全黑或全白**  
**Problem: Generated stimulus is all black or white**

→ 检查对比度参数范围 (0-1)  
→ Check contrast parameter range (0-1)

**问题: Gabor看不见**  
**Problem: Gabor not visible**

→ 增加目标对比度 `tgContrast`  
→ 减小背景对比度 `bgContrast`  
→ Increase target contrast `tgContrast`  
→ Decrease background contrast `bgContrast`

**问题: 空间频率警告**  
**Problem: Spatial frequency warning**

→ 降低 `GaborSF` 值  
→ 增加屏幕分辨率  
→ Decrease `GaborSF` value  
→ Increase screen resolution

---

## 🔬 实验整合 / Experimental Integration

### 与Psychtoolbox整合 / Integration with Psychtoolbox

```matlab
% 1. 预先生成所有刺激
% Pre-generate all stimuli
stimuli = cell(1, nTrials);
for i = 1:nTrials
    stimuli{i} = genStim(...);
end

% 2. 在实验中使用
% Use in experiment
Screen('Preference', 'SkipSyncTests', 1);
[window, rect] = PsychImaging('OpenWindow', 0, 128);

for trial = 1:nTrials
    % 创建纹理
    % Create texture
    tex = Screen('MakeTexture', window, stimuli{trial} * 255);
    
    % 显示
    % Display
    Screen('DrawTexture', window, tex);
    Screen('Flip', window);
    
    % ... 实验逻辑 / Experiment logic
    
    Screen('Close', tex);
end

sca;
```

---

## 📊 示例输出 / Example Output

运行 `demo_StimulusGenerator.m` 将生成：

Running `demo_StimulusGenerator.m` will generate:


1. **单刺激展示**：
    - `stimulus`：完整刺激（Gabor+背景）
    - `background_circular`：圆形1/f噪声背景
    - `gabor_circular`：圆形Gabor目标

2. （如需批量或参数展示，可参考demo脚本扩展）

---

## 📝 引用 / Citation

如果在研究中使用本工具包，请引用：

If you use this toolkit in your research, please cite:

```
Lu, Y. (2025). Stimulus Generator Toolkit for Visual Psychophysics. 
Institute of Psychology, Chinese Academy of Sciences.
```

---

## 📧 联系与支持 / Contact & Support

- **作者 / Author:** Yuxin Lu  
- **邮箱 / Email:** luyx@psych.ac.cn  
- **机构 / Institution:** IPCAS (Institute of Psychology, Chinese Academy of Sciences)

有问题或建议？请发送邮件！  
Questions or suggestions? Please send an email!

---

## 📜 许可证 / License

本工具包供学术研究使用。  
This toolkit is for academic research use.

---

## 🔄 更新日志 / Changelog

### v1.0 (2025.1.7)
- ✅ 初始版本发布 / Initial release
- ✅ 核心刺激生成功能 / Core stimulus generation
- ✅ 完整中英文文档 / Complete bilingual documentation
- ✅ Demo示例脚本 / Demo example script

---

**祝实验顺利！/ Good luck with your experiments! 🎯**
