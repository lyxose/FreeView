# %% [1] 导入库
import json
from pathlib import Path
import pandas as pd  # 如果没有安装，可注释掉并使用纯列表方式（见末尾附注）

# %% [2] 加载 JSON 文件（请将路径改为你的实际文件）
file_path = Path("E:\\AllDownloads\\DOWNLOAD\\WebFreeView_5_Ses1_20260512T142151.json")
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

# %% [8] 根据第一次刺激的视频时间对齐，换算所有试次相对于视频 00:00:00.000 的时间点

# 8.1 定义时间字符串解析函数（输入格式：hh:mm:ss.fff）
def parse_videotime_to_ms(video_time_str):
    """
    将 'hh:mm:ss.fff' 字符串转换为毫秒整数。
    例如 '00:02:02.983' → 122983。
    支持带负号的时间（极少使用），返回负毫秒数。
    """
    s = video_time_str.strip()
    sign = 1
    if s.startswith('-'):
        sign = -1
        s = s[1:]

    # 分离小数部分（毫秒）
    if '.' in s:
        time_part, frac = s.split('.')
        ms = int(frac.ljust(3, '0')[:3])  # 取前三位，不足补零
    else:
        time_part = s
        ms = 0

    # 分离时、分、秒
    h, m, sec = map(int, time_part.split(':'))
    total_ms = ((h * 3600) + (m * 60) + sec) * 1000 + ms
    return sign * total_ms

# 8.2 输入第一次刺激出现的视频时间（可根据实际情况修改此字符串）
video_first_stim_str = "00:00:08.335"  # 示例值，请按需修改
# video_first_stim_str = "00:00:10.809"  # 示例值，请按需修改

video_first_stim_time_ms = parse_videotime_to_ms(video_first_stim_str)

# 8.3 获取第一个试次的刺激出现时间（Unix 毫秒）
first_stim_epoch_ms = df_epoch.loc[0, "stimulus_time"]  # trialIndex=0

# 8.4 计算偏移量（实验绝对时间 - 视频相对时间）
offset_ms = first_stim_epoch_ms - video_first_stim_time_ms

# 8.5 将所有试次的时间转换为相对于视频 00:00:00.000 的毫秒数
df_epoch["stimulus_video_ms"] = df_epoch["stimulus_time"] - offset_ms
df_epoch["response_video_ms"] = df_epoch["response_time"] - offset_ms

# 8.6 毫秒转视频时间字符串（强制 hh:mm:ss.fff）
def ms_to_videotime(ms):
    """将视频内毫秒数转换为 hh:mm:ss.fff 字符串，NaN 返回 None"""
    if pd.isna(ms):
        return None
    ms = int(ms)
    sign = "-" if ms < 0 else ""
    ms = abs(ms)
    total_sec = ms // 1000
    millis = ms % 1000
    minutes = total_sec // 60
    seconds = total_sec % 60
    hours = minutes // 60
    minutes = minutes % 60
    return f"{sign}{hours:02d}:{minutes:02d}:{seconds:02d}.{millis:03d}"

# 8.7 生成可读列
df_epoch["stimulus_video"] = df_epoch["stimulus_video_ms"].apply(ms_to_videotime)
df_epoch["response_video"] = df_epoch["response_video_ms"].apply(ms_to_videotime)

# 8.8 显示结果
print("=== 刺激 & 响应时间相对于视频 00:00:00.000 ===")
print(f"对齐锚点：第一个刺激的视频时间 = {video_first_stim_str}")
print(df_epoch[["trial", "stimulus_video", "response_video", "rt_ms", "responded"]])

# %%
