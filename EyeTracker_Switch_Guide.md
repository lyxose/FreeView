# 眼动仪切换使用指南

## 概述
本项目已重构以支持两种眼动仪设备：
- **Tobii Pro Fusion** (通过 TITTA 工具箱)
- **EyeLink** (通过 EyeLink Toolbox)

## 文件结构

### 主程序
- `FreeView_main.m` - 主实验控制脚本，包含公共参数和逻辑

### 设备特定脚本
- `FreeViewExp_PTB_TITTA.m` - Tobii 眼动仪实验脚本
- `FreeViewExp_PTB_Eyelink.m` - EyeLink 眼动仪实验脚本

## 切换眼动仪设备

### 方法：修改主程序参数

打开 `FreeView_main.m`，找到第 54 行左右的眼动仪选择部分：

```matlab
%% Eye tracker selection
% Set eyeTrackerType to 'Tobii' or 'EyeLink'
eyeTrackerType = 'Tobii';  % Change to 'EyeLink' to use EyeLink eye tracker
```

**使用 Tobii:**
```matlab
eyeTrackerType = 'Tobii';
```

**使用 EyeLink:**
```matlab
eyeTrackerType = 'EyeLink';
```

## 数据保存说明

### 共同数据（两种设备都会保存）
- `./Data/Formal/Result_Sub[ID]_Ses[Session]_[Location]_[Name]_[DateTime].mat` - 实验结果数据
- `./Data/Formal/EXP_Sub[ID]_Ses[Session]_[Location]_[Name]_[DateTime].mat` - 完整实验数据 (如果 saveRaw = true)

### Tobii 特定数据
- `./Data/Formal/Dat_Sub[ID]_Ses[Session].mat` - Tobii 眼动数据（通过 TITTA 保存）

### EyeLink 特定数据
- `./Data/Formal/S[ID]S[Session]_Sub[ID]_Ses[Session]_[DateTime].edf` - EyeLink 原始数据文件

## 设备特定参数

### Tobii (在 FreeViewExp_PTB_TITTA.m 中)
```matlab
tobiiFreq = 250;              % 采样率 250 Hz
tobiiMod = 'Tobii Pro Fusion'; % 设备型号
```

### EyeLink (在 FreeViewExp_PTB_Eyelink.m 中)
```matlab
eyelinkFreq = 1000;  % 采样率，可选 250, 500, 1000, 2000 Hz
dummymode = 0;       % 0 = 真实连接，1 = 测试模式
```

## 注意事项

### 显示器尺寸配置
如果使用 EyeLink，请在 `FreeViewExp_PTB_Eyelink.m` 中调整显示器物理尺寸：

```matlab
monWidth = 52.7;   % 显示器宽度（厘米）
monHeight = 29.7;  % 显示器高度（厘米）
```

### EyeLink 文件名限制
- EDF 文件名最多 8 个字符
- 只能包含字母、数字和下划线
- 当前格式：`S[SubjID]S[Session]` (例如：S1S1.edf)

### 校准设置
- Tobii: 支持动画校准点和双眼分别校准
- EyeLink: 使用 HV9 校准模式（9点水平-垂直）

## 实验流程对比

| 步骤 | Tobii | EyeLink |
|------|-------|---------|
| 初始化 | Titta.getDefaults() | EyelinkInit() |
| 校准 | EThndl.calibrate() | EyelinkDoTrackerSetup() |
| 开始记录 | EThndl.buffer.start() | Eyelink('StartRecording') |
| 获取数据 | EThndl.buffer.peekTimeRange() | Eyelink('NewestFloatSample') |
| 发送消息 | EThndl.sendMessage() | Eyelink('Message') |
| 停止记录 | 自动 | Eyelink('StopRecording') |
| 保存数据 | EThndl.saveData() | Eyelink('ReceiveFile') |
| 关闭 | EThndl.deInit() | Eyelink('Shutdown') |

## 故障排除

### Tobii 连接问题
1. 确认 TITTA 工具箱已安装
2. 检查 Tobii Pro Eye Tracker Manager 是否运行
3. 验证设备型号名称是否正确

### EyeLink 连接问题
1. 确认 EyeLink Developers Kit 已安装
2. 检查以太网连接（默认 IP: 100.1.1.1）
3. 如需测试，设置 `dummymode = 1`

### 校准失败
- Tobii: 检查光照条件，确保被试正确定位
- EyeLink: 调整相机图像，确保瞳孔和角膜反射清晰

## 未来扩展

如需添加其他眼动仪设备：
1. 创建新的设备脚本（如 `FreeViewExp_PTB_[DeviceName].m`）
2. 在 `FreeView_main.m` 的 eyeTrackerType 选择部分添加新分支
3. 实现相同的实验流程和数据接口

## 联系与支持

如有问题或建议，请查阅：
- TITTA 文档: https://github.com/dcnieho/Titta
- EyeLink 文档: SR Research Support Forum
