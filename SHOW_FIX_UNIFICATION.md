# show_fix 函数统一更新说明

## 总体修改

将 `show_fix.m` 函数升级以支持 **Tobii（TITTA）** 和 **EyeLink** 两种眼动仪设备，并统一在两个实验脚本中使用该函数来呈现注视点，确保实验逻辑一致。

## 修改内容

### 1. show_fix.m 函数改进

**新增功能：**
- 支持两种眼动仪设备（Tobii 和 EyeLink）
- 自动检测眼动仪类型
- 向后兼容旧的函数签名

**新的函数签名：**
```matlab
[startT, headDist] = show_fix(wpnt, x, y, fixTime, fixClrs, winRect, varargin)
```

**参数（使用 name-value 对）：**
- `'eyeTrackerType'`: 'Tobii' | 'EyeLink' | 'None'（自动检测）
- `'EThndl'`: Titta 对象（用于 Tobii）
- `'el'`: EyeLink 配置结构体（用于 EyeLink）
- `'tobiiFreq'`: Tobii 采样频率（默认 250 Hz）
- `'monWidth'`, `'monHeight'`: 显示器物理尺寸（默认 51.1 × 28.7 cm）

**设备特定行为：**

| 设备 | 行为 |
|------|------|
| Tobii | 显示注视点，验证眼睛是否正确固定在目标位置，如未固定则重试 |
| EyeLink | 显示注视点，等待指定时间（可选：添加注视验证） |
| 无设备 | 显示注视点，等待指定时间 |

### 2. FreeViewExp_PTB_Eyelink.m 更新

**修改位置 1：练习试次注视点呈现**
```matlab
% 旧方法：
Screen('DrawDots', wpnt, fixCenter, 14, fixClrs(1), [], 3);
Screen('DrawDots', wpnt, fixCenter, 7, fixClrs(2), [], 3);
startT = Screen('Flip', wpnt);
WaitSecs(fixTime);

% 新方法：
[startT, headDist] = show_fix(wpnt, fixCenter(1), fixCenter(2), fixTime, fixClrs, winRect, ...
                               'eyeTrackerType', 'EyeLink', 'el', el, ...
                               'monWidth', monWidth, 'monHeight', monHeight);
```

**修改位置 2：正式试次注视点呈现**
```matlab
% 统一使用 show_fix 函数替换直接绘制和等待代码
```

### 3. FreeViewExp_PTB_TITTA.m 更新

**从旧签名更新为新签名：**
```matlab
% 旧方法：
[startT, headDist] = show_fix(wpnt, fixCenter(1), fixCenter(2), fixTime, fixClrs, winRect, ...
                               true, tobiiFreq, EThndl, monWidth, monHeight);

% 新方法：
[startT, headDist] = show_fix(wpnt, fixCenter(1), fixCenter(2), fixTime, fixClrs, winRect, ...
                               'eyeTrackerType', 'Tobii', 'EThndl', EThndl, ...
                               'tobiiFreq', tobiiFreq, 'monWidth', monWidth, 'monHeight', monHeight);
```

修改了两处调用位置：
1. 练习试次（line ~189）
2. 正式实验循环（line ~312）

## 优势

### 代码一致性
- 两个设备现在使用相同的注视点呈现函数
- 统一的参数接口和返回值
- 易于后续添加新设备支持

### 可维护性
- 设备特定逻辑集中在 `show_fix.m`
- 实验脚本中无需重复设备检查代码
- 更新注视点呈现逻辑只需修改一个函数

### 可扩展性
- 易于为 EyeLink 添加更复杂的注视点验证
- 支持混合设备配置
- 名字-值对参数允许无痛扩展

## 使用示例

### 对于 Tobii（TITTA）：
```matlab
[startT, headDist] = show_fix(wpnt, centerX, centerY, fixTime, fixClrs, winRect, ...
                               'eyeTrackerType', 'Tobii', ...
                               'EThndl', EThndl, ...
                               'tobiiFreq', 250, ...
                               'monWidth', 51.1, 'monHeight', 28.7);
```

### 对于 EyeLink：
```matlab
[startT, headDist] = show_fix(wpnt, centerX, centerY, fixTime, fixClrs, winRect, ...
                               'eyeTrackerType', 'EyeLink', ...
                               'el', el, ...
                               'monWidth', 51.1, 'monHeight', 28.7);
```

### 自动检测（推荐）：
```matlab
% 如果 EThndl 或 el 有效，会自动选择正确的设备
[startT, headDist] = show_fix(wpnt, centerX, centerY, fixTime, fixClrs, winRect, ...
                               'EThndl', EThndl, 'el', el, ...
                               'tobiiFreq', 250, ...
                               'monWidth', 51.1, 'monHeight', 28.7);
```

## 向后兼容性

函数仍然支持旧的位置参数调用方式（如 `show_fix(wpnt, x, y, fixTime, fixClrs, winRect, true, tobiiFreq, EThndl, ...)`），但 **强烈建议使用新的 name-value 格式** 以获得更好的可读性和灵活性。

## 未来改进

可以考虑为 EyeLink 添加：
- 注视点有效性检查（类似 Tobii 的中心固定验证）
- 重新校准触发器（如果检测到漂移）
- 头部位置反馈
- 实时眼动数据日志记录

## 文件列表

修改的文件：
- `function_library_cus/show_fix.m` - 主函数更新
- `FreeViewExp_PTB_Eyelink.m` - 两处调用位置更新
- `FreeViewExp_PTB_TITTA.m` - 两处调用位置更新
