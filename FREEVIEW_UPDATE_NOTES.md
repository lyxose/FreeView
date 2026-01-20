# FreeView 被试信息管理更新说明

## 更新日期
2026-01-20

## 主要改进

### 1. 解决按键堆积问题
**问题**：在 PTB 界面后使用 `input()` 会导致之前的按键堆积到输入框中，无法正确输入被试编号。

**解决方案**：
- 移除所有 `input()` 命令行输入
- 统一使用 `inputdlg()` 对话框输入
- 在弹出对话框前调用 `ListenChar(0)` 释放键盘监听
- 对话框关闭后调用 `ListenChar(2)` 恢复键盘监听

### 2. 改进传参模式
**旧模式**：
```matlab
run_freeview_eyelink(wpnt, winRect, el, edfFile, subjID, session, ...)
```

**新模式**：
```matlab
% 创建被试信息结构体
subjInfo = struct('subjID', subjID, 'session', session, ...
                  'location', location, 'subjName', subjName, ...
                  'subjGender', subjGender, 'subjAge', subjAge);
                  
% 调用时传入结构体
run_freeview_eyelink(wpnt, winRect, el, edfFile, subjInfo, ...)
```

**优点**：
- 参数更清晰，易于扩展
- 支持部分信息传入（字段可选）
- 避免参数顺序错误

### 3. 强制确认机制
**行为**：无论是否传入被试信息，始终弹出对话框让主试确认或修改。

**特性**：
- 传入完整信息 → 对话框预填充，主试可直接确认或修改
- 传入部分信息 → 预填充已有信息，其余字段使用默认值或历史数据
- 不传入信息 → 根据 CSV 历史智能建议（新被试号、递增 session 等）

## 修改的文件

### 1. `run_freeview_eyelink.m`
**主要改动**：
- 函数签名：`subjID, session` → `subjInfo` (struct)
- Standalone 模式：移除 `input()`，改用 `InformationBox_FV(struct())`
- 非 standalone 模式：始终调用 `InformationBox_FV(subjInfo)` 进行确认

### 2. `InformationBox_FV.m`
**主要改动**：
- 接受 `subjInfo` 结构体参数（支持空结构体或部分字段）
- 始终显示对话框（包含所有字段：subjID, session, location, name, gender, age, date）
- 智能预填充逻辑：
  - 若 subjID 存在于 CSV → 加载历史信息，session 自动 +1（除非传入了 session）
  - 若 subjID 不存在 → 使用传入值或默认值
  - 若未传入 subjID → 建议下一个编号

### 3. `FreeViewExp_PTB_Eyelink.m`
**主要改动**：
```matlab
% 旧代码
run_freeview_eyelink(wpnt, winRect, el, edfFile, [], [], ...)

% 新代码
subjInfo = struct('subjID', subjID, 'session', session, ...
                  'location', location, 'subjName', subjName, ...
                  'subjGender', subjGender, 'subjAge', subjAge);
run_freeview_eyelink(wpnt, winRect, el, edfFile, subjInfo, ...)
```

## 使用示例

### 示例 1: Standalone 模式
```matlab
run_freeview_eyelink()
% 自动弹出对话框，主试输入所有信息
```

### 示例 2: 从主实验调用（完整信息）
```matlab
subjInfo = struct('subjID', 1, 'session', 2, ...
                  'location', 'IP', 'subjName', 'ZhangSan', ...
                  'subjGender', 'M', 'subjAge', 25);
run_freeview_eyelink(wpnt, winRect, el, edfFile, subjInfo, ...);
% 弹出对话框，预填充所有信息，主试确认即可
```

### 示例 3: 从主实验调用（部分信息）
```matlab
subjInfo = struct('subjID', 1);
run_freeview_eyelink(wpnt, winRect, el, edfFile, subjInfo, ...);
% 弹出对话框，预填充 subjID，其他字段根据 CSV 历史智能填充
```

### 示例 4: 从主实验调用（空结构体）
```matlab
subjInfo = struct();
run_freeview_eyelink(wpnt, winRect, el, edfFile, subjInfo, ...);
% 弹出对话框，根据 CSV 建议下一个被试号
```

## 键盘监听管理流程

```
PTB 实验界面 (ListenChar(2) 激活)
    ↓
调用 run_freeview_eyelink()
    ↓
ListenChar(0) - 释放键盘监听
    ↓
InformationBox_FV() - 显示对话框，安全输入
    ↓
ListenChar(2) - 恢复键盘监听
    ↓
继续 FreeView 刺激呈现
```

## 数据文件
- `./Data/SubjInfo_FV.csv` - FreeView 被试信息汇总表
- `./Data/FV/` - FreeView EDF 数据输出目录

## 兼容性说明
此次更新修改了 `run_freeview_eyelink()` 的参数列表，需要同步更新所有调用点：
- ✅ `FreeViewExp_PTB_Eyelink.m` - 已更新
- ⚠️ 如有其他脚本调用此函数，需手动更新为新参数格式

## 测试建议
1. 测试 standalone 模式：直接运行 `run_freeview_eyelink()`
2. 测试完整信息传入：确认对话框正确预填充
3. 测试部分信息传入：确认智能补全逻辑
4. 测试空结构体传入：确认建议下一个被试号
5. 测试在 PTB 界面后调用：确认不受按键堆积影响
