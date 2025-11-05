/FreeView_Paper
│  stat_Nfix_time.m
│  File_struct.md
│  
├─EXP_design   ← 实验设计（paradigms、scriptes、被试表）
│  │  paradigms.md
│  │  SubjInfo_v1.5.csv
│  │  SubjInfo_v1.csv
│  │  SubjInfo_v2.csv
│  │  
│  └─Scripts
│          FreeViewExp_PTB_TITTA.m ← 仅在FreeView_main_v2.m中会调用此脚本
│          FreeView_main_v1.5.m
│          FreeView_main_v1.m
│          FreeView_main_v2.m
│          getThreshold_v1.5.m
│          getThreshold_v1.m
│
├─EXP_notes
│      EXP1.5_notes.md
│      EXP1_notes.md
│      EXP2_notes.md
│
├─function_library_cus    ← 辅助函数
│  ├─ANA  ← 分析辅助函数
│  │      analyze_angle_curve.m
│  │      analyze_axis_obli_correlation.m
│  │      analyze_sector_fisherz.m
│  │      analyze_subject_spectrum.m
│  │      build_fixTable.m
│  │      cluster_based_permutation_test.m
│  │      compute_gaussian_window_series.m
│  │      compute_sliding_window_series.m
│  │      CrossV_collect_data.m
│  │      CrossV_plot_bar_anova.m
│  │      fdr_bh.m
│  │      filter_fixTable_for_analysis.m
│  │      gaussian_smooth_along_dim.m
│  │      get_eye_data_files.m
│  │      heatmap_stat_time.m
│  │      make_fixation_sector_video.m
│  │      normalize_by_dim.m
│  │      plot_angle_curve.m
│  │      plot_bar_multi.m
│  │      plot_comparison.m
│  │      plot_fixTable_distributions.m
│  │      plot_fixTable_heatmap.m
│  │      plot_sector_fisherz_clusters.m
│  │      plot_sector_stats_time_series.m
│  │      plot_sector_time_heatmap.m
│  │      plot_single.m
│  │      plot_single_prop.m
│  │      plot_sliding_window_analysis.m
│  │      plot_sliding_window_analysis_corr.m
│  │      plot_subject_spectrum.m
│  │      plot_trial_distributions.m
│  │      process_fixation_count_series.m
│  │      process_fix_density_series.m
│  │      process_time_count_series.m
│  │      process_time_density_series.m
│  │      save_subject_density_video.m
│  │      single_subject_fix_density.m
│  │      single_subject_time_density.m
│  │      UT.m
│  │
│  └─EXP   ← 实验辅助函数
│          checkend.m
│          demoCalCompletionFun.m
│          disp_rest.m
│          drawcenteredtext_dot.m
│          drawCentImg.m
│          expRun.m
│          genStim.m
│          getHeadDist.m
│          getLastFix.m
│          grating.m
│          InformationBox.m
│          overlapGabor.m
│          saveStiImg.m
│          showInstruc.m
│          show_fix.m
│          spatialJudge.m
│          tPinkNoise.m
│          UT.m
│          winOverlap.m
│
└─Results
    │  Axis_Effect.png
    │  Obli_Effect.png
    │
    ├─v1
    │      ALL_fixTable_26Subj.csv
    │      v1 (v1.1)--16扇区分别统计柱状图.png
    │      v1 (v1.1)--16扇区分别统计柱状图（比例）.png
    │      v1 (v1.1)--16扇区变异性时程.png
    │      v1 (v1.1)--16扇区热图.png
    │      v1 (v1.1)--45°fold扫描.png
    │      v1 (v1.1)--90°fold扫描.png
    │      v1 (v1.1)--AG总效应.png
    │      v1 (v1.1)--AG效应.png
    │      v1 (v1.1)--AG时程.png
    │      v1 (v1.1)--CO总效应.png
    │      v1 (v1.1)--CO效应.png
    │      v1 (v1.1)--CO时程.png
    │      v1 (v1.1)--trial滑窗全时程.png
    │      v1 (v1.1)--trial高斯滑窗.png
    │      v1 (v1.1)--全时窗热图.png
    │      v1 (v1.1)--全角度扫描.png
    │
    ├─v1.5
    │      ALL_fixTable_26Subj.csv
    │      v1.5 (v1.2)--16扇区分别统计柱状图.png
    │      v1.5 (v1.2)--16扇区分别统计柱状图（比例）.png
    │      v1.5 (v1.2)--16扇区变异性时程.png
    │      v1.5 (v1.2)--16扇区热图.png
    │      v1.5 (v1.2)--45°fold扫描.png
    │      v1.5 (v1.2)--90°fold扫描.png
    │      v1.5 (v1.2)--AG总效应.png
    │      v1.5 (v1.2)--AG效应.png
    │      v1.5 (v1.2)--AG时程.png
    │      v1.5 (v1.2)--CO总效应.png
    │      v1.5 (v1.2)--CO效应.png
    │      v1.5 (v1.2)--CO时程.png
    │      v1.5 (v1.2)--trial滑窗全时程.png
    │      v1.5 (v1.2)--trial高斯滑窗.png
    │      v1.5 (v1.2)--全时窗热图.png
    │      v1.5 (v1.2)--全角度扫描.png
    │
    └─v2
            ALL_fixTable_25Subj.csv
            v2 (v1.3)--16扇区分别统计柱状图.png
            v2 (v1.3)--16扇区分别统计柱状图（比例）.png
            v2 (v1.3)--16扇区变异性时程.png
            v2 (v1.3)--16扇区热图.png
            v2 (v1.3)--45°fold扫描.png
            v2 (v1.3)--90°fold扫描.png
            v2 (v1.3)--AG总效应.png
            v2 (v1.3)--AG效应.png
            v2 (v1.3)--AG时程.png
            v2 (v1.3)--CO总效应.png
            v2 (v1.3)--CO效应.png
            v2 (v1.3)--CO时程.png
            v2 (v1.3)--trial滑窗全时程.png
            v2 (v1.3)--trial高斯滑窗.png
            v2 (v1.3)--全时窗热图.png
            v2 (v1.3)--全角度扫描.png
            v2 (v1.3)--各扇区分别统计柱状图.png
