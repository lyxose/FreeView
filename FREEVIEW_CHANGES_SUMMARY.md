# FreeView EyeLink - 被试信息管理和输出位置调整总结

## 做出的主要更改

### 1. 创建新函数：`InformationBox_FV.m`
- **位置**：`function_library_cus/InformationBox_FV.m`
- **功能**：专门用于 FreeView 实验的被试信息管理
- **数据存储**：`./Data/SubjInfo_FV.csv`

**主要特性**：
- 支持两种输入模式：
  - 若提供 `subjID` 和 `session`，则只弹窗请求其他信息（location, name, gender, age, date）
  - 若未提供，则从 CSV 读取并自动建议下一个被试号，然后弹窗请求全部信息
- 自动创建/更新 `SubjInfo_FV.csv`，字段包括：
  ```
  subjID, session, location, subjName, subjGender, subjAge, date, Note
  ```
- 若被试已存在，自动建议下一 session 号

### 2. 修改 `run_freeview_eyelink.m`

#### **Standalone 模式调整**：
```matlab
% 旧逻辑：直接 inputdlg 获取 subjID, session, suffix
% 新逻辑：
1. input() 获取 subjID（命令行）
2. ListenChar(0) - 释放键盘监听
3. 调用 InformationBox_FV(subjID, []) 获取其他信息
4. ListenChar(2) - 恢复键盘监听
5. inputdlg 获取 FV EDF suffix
```

#### **非 Standalone 模式调整**：
- 新增被试信息检查逻辑（第 100-110 行）
- 若 `subjID` 或 `session` 为空，自动调用 `InformationBox_FV()` 弹窗
- 使用 `ListenChar()` 控制键盘监听，确保 dialog 能正常输入

#### **输出位置调整**：
- **旧位置**：`./Data/Formal/`
- **新位置**：`./Data/FV/`
- 自动创建 `./Data/FV/` 目录

### 3. 键盘监听管理
在弹窗输入前后添加了 `ListenChar()` 调用：
- `ListenChar(0)`：释放键盘监听，允许 inputdlg 和 input() 接收用户输入
- `ListenChar(2)`：恢复键盘监听（PsychToolbox 实验模式），用于刺激呈现

这确保了在 dialog 弹窗期间用户能正常输入，而进入实验刺激程序后键盘被主程序控制。

## 使用方式

### Standalone 调用
```matlab
run_freeview_eyelink()
% 依次输入：
% 1. 被试 ID（命令行）
% 2. Session, location, name, gender, age, date（dialog）
% 3. EDF suffix（dialog）
```

### 集成调用（从其他实验中调用）
```matlab
% 被试信息完整的调用
run_freeview_eyelink(wpnt, winRect, el, edfFile, subjID, session, ...)

% 被试信息不完整的调用（会自动弹窗）
run_freeview_eyelink(wpnt, winRect, el, edfFile, [], [], ...)
```

## 数据文件格式

### `./Data/SubjInfo_FV.csv`
CSV 格式，包含所有 FreeView 被试的汇总信息：
```
subjID,session,location,subjName,subjGender,subjAge,date,Note
1,1,IP,MingZipinying,M,22,20260120,
2,1,IP,ZhangMingming,F,24,20260120,
```

## 向后兼容性
- 若未在 standalone 模式下使用，现有代码继续支持直接传入 `subjID` 和 `session` 参数
- EDF 文件输出位置从 `./Data/Formal` 改为 `./Data/FV`，需要相应调整数据分析脚本
