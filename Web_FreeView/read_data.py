# %% [1] 导入库
import json
from pathlib import Path
import pandas as pd  # 如果没有安装，可注释掉并使用纯列表方式（见末尾附注）

# %% [2] 加载 JSON 文件（请将路径改为你的实际文件）
file_path = Path("E:\\AllDownloads\\DOWNLOAD\\WebFreeView_8_Ses1_20260430T111247.json")
with open(file_path, "r", encoding="utf-8") as f:
    data = json.load(f)

print(f"加载成功！共 {len(data['trials'])} 个试次。")

# %% [3] 定义时间提取函数
def get_trial_times(trial, time_format="epoch_ms"):
    """
    从单个 trial 字典中提取刺激呈现和按键响应的时间。
    
    time_format 可选值：
        'epoch_ms'      : Unix 毫秒时间戳（整数）
        'epoch_s'       : Unix 秒时间戳（浮点数）
        'iso'           : ISO 8601 字符串，如 '2026-04-30T03:12:20.647Z' (UTC)
        'local'         : 本地完整字符串，如 '2026/4/30 11:12:20'
        'local_compact' : 本地紧凑字符串，如 '2026-04-30 11:12:20.647'
        'utc'           : UTC 字符串，如 'Thu, 30 Apr 2026 03:12:20 GMT'
    """
    stim = trial.get("stimulus")
    resp = trial.get("response")

    def extract_time(time_obj, fmt):
        if time_obj is None:
            return None
        if fmt == "epoch_ms":
            return time_obj["epochMs"]
        elif fmt == "epoch_s":
            return time_obj["epochSeconds"]
        elif fmt == "iso":
            return time_obj["iso"]
        elif fmt == "local":
            return time_obj["localString"]
        elif fmt == "local_compact":
            return time_obj["localCompact"]
        elif fmt == "utc":
            return time_obj["utcString"]
        else:
            raise ValueError(f"不支持的时间格式: {fmt}")

    return {
        "trial": trial["trialIndex"],
        "stimulus_time": extract_time(stim, time_format),
        "response_time": extract_time(resp, time_format),
        "responded": trial.get("responded", False),
        "rt_ms": trial.get("responseRTMs")  # 反应时间（毫秒）
    }

# %% [4] 构建时间戳表格（使用 DataFrame，方便查看）
def build_timestamp_table(data, time_format="epoch_ms"):
    rows = [get_trial_times(t, time_format) for t in data["trials"]]
    return pd.DataFrame(rows)

# %% [5] 示例：以 Unix 毫秒时间戳显示
df_epoch = build_timestamp_table(data, "epoch_ms")
print("=== 刺激 & 响应时间（Unix 毫秒）===")
print(df_epoch)

# %% [6] 示例：以 ISO 8601 格式显示（人类可读，UTC）
df_iso = build_timestamp_table(data, "iso")
print("=== 刺激 & 响应时间（ISO 8601 UTC）===")
print(df_iso)

# %% [7] 示例：以本地紧凑格式显示
df_local = build_timestamp_table(data, "local_compact")
print("=== 刺激 & 响应时间（本地时间）===")
print(df_local)

# %% [8] （可选）直接访问某个 trial 的原始时间戳对象
trial_index = 1  # 第一个 trial 的 trialIndex = 1
trial = data["trials"][trial_index - 1]  # 转为 0‑based 索引

print(f"\nTrial {trial['trialIndex']} 详细信息：")
print(f"  刺激出现 (epochMs): {trial['stimulus']['epochMs']}")
print(f"  刺激出现 (ISO)  : {trial['stimulus']['iso']}")
print(f"  按键响应 (epochMs): {trial['response']['epochMs'] if trial['response'] else '无'}")
print(f"  按键响应 (ISO)  : {trial['response']['iso'] if trial['response'] else '无'}")
print(f"  反应时间 (ms)   : {trial['responseRTMs']:.1f}")

# %% [9] （可选）整理成纯 Python 列表，无需 pandas
def get_all_times_simple(data, fmt="epoch_ms"):
    """返回 [(trialIndex, stim_time, resp_time, rt_ms), ...] 的列表"""
    results = []
    for trial in data["trials"]:
        info = get_trial_times(trial, fmt)
        results.append((
            info["trial"],
            info["stimulus_time"],
            info["response_time"] if info["responded"] else None,
            info["rt_ms"]
        ))
    return results

times_list = get_all_times_simple(data, "iso")
print("\n=== 纯列表形式（ISO 格式）===")
for t in times_list:
    print(t)
# %%
