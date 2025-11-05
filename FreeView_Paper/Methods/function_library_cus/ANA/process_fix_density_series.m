function [sub_fix_bin_density, nfix_FT] = process_fix_density_series( ...
    start_FT, dur_FT, xpos_FT, ypos_FT, dnfix_FT, sub_FT, ses_FT, pairs_FT, ...
    video_width, video_height, ut, cfg)
% PROCESS_FIX_DENSITY_SERIES 计算每个被试按注视序号的空间密度分布
%
% 输入参数:
%   start_FT, dur_FT    - 注视开始时间和持续时间 (ms)（用于可选的视频生成）
%   xpos_FT, ypos_FT    - 注视点坐标 (pixel)
%   dnfix_FT            - 每个注视在trial内的剔除后序号（1..n）
%   sub_FT, ses_FT      - 被试ID和会话ID
%   pairs_FT            - [nsbj x 2] 被试-会话对
%   video_width, video_height - 视频尺寸 (pixel)
%   ut                  - 工具类（用于度数-像素转换）
%   cfg                 - 配置结构体，包含:
%       .keep_nFix      : 注视序号上限（默认11）
%       .sigma_deg      : 高斯核标准差 (deg，默认1)
%       .n_sectors      : 扇区数量 (默认16)
%       .sector_edges   : 扇区边界 (deg，默认-11.25:22.5:348.75)
%       .write_video    : 是否写视频文件 (默认false)
%       .useParallel    : 是否使用并行计算 (默认true)
%
% 输出参数:
%   sub_fix_bin_density - [nsbj x nfix x n_sectors] 去趋势后的密度
%   nfix_FT             - 注视序号上限

    % 参数解析与默认值
    if nargin < 12, cfg = struct(); end
    if ~isfield(cfg, 'keep_nFix'), cfg.keep_nFix = 11; end
    if ~isfield(cfg, 'sigma_deg'), cfg.sigma_deg = 1; end
    if ~isfield(cfg, 'n_sectors'), cfg.n_sectors = 16; end
    if ~isfield(cfg, 'write_video'), cfg.write_video = false; end
    if ~isfield(cfg, 'useParallel'), cfg.useParallel = true; end
    
    nfix_FT = max(1, round(cfg.keep_nFix));
    sigma_deg = cfg.sigma_deg;
    n_sectors = cfg.n_sectors;
    write_video = cfg.write_video;
    useParallel = cfg.useParallel;
    
    if ~isfield(cfg, 'sector_edges')
        cfg.sector_edges = linspace(-11.25, 360-11.25, n_sectors+1);
    end
    sector_edges = cfg.sector_edges;
    
    % 计算参数
    sigma_pix = max(0.5, ut.deg2pix(sigma_deg));
    nsbj = size(pairs_FT, 1);
    
    % 构建高斯核
    radius = max(1, ceil(3 * sigma_pix));
    [gx, gy] = meshgrid(-radius:radius, -radius:radius);
    gauss_template = exp(-(gx.^2 + gy.^2) / (2*sigma_pix^2));
    gauss_template = single(gauss_template);
    
    % 预计算扇区映射
    center_pix = [video_width, video_height] / 2;
    [xg, yg] = meshgrid(1:video_width, 1:video_height);
    angles_grid = mod(atan2d(yg - center_pix(2), xg - center_pix(1)), 360);
    sector_map = zeros(video_height, video_width, 'uint8');
    for si = 1:n_sectors
        l = mod(sector_edges(si), 360);
        r = mod(sector_edges(si+1), 360);
        if l < r
            sector_map(angles_grid >= l & angles_grid < r) = si;
        else
            sector_map(angles_grid >= l | angles_grid < r) = si;
        end
    end
    
    % 输出容器初始化
    fix_bin_sub_density = zeros(nfix_FT, n_sectors, nsbj, 'single');
    
    % 并行处理检查
    if useParallel
        try
            if isempty(gcp('nocreate'))
                parpool;
            end
        catch
            warning('Parallel pool not available. Falling back to serial execution.');
            useParallel = false;
        end
    end
    
    % 时间统计
    tAll = tic;
    timeStr = char(datetime("now","Format","HH:mm:ss"));
    fprintf('[%s] Subject-level fixation density: %d subjects, %d fixations, image %dx%d, sigma=%.2f deg (%.2f px)\n', ...
        timeStr, nsbj, nfix_FT, video_width, video_height, sigma_deg, sigma_pix);
    
    subj_dur = nan(nsbj, 1);
    
    % 主循环：逐被试处理
    if useParallel
        parfor si = 1:nsbj
            [fix_bin_sub_density(:,:,si), subj_dur(si)] = single_subject_fix_density(...
                si, pairs_FT, sub_FT, ses_FT, dnfix_FT, xpos_FT, ypos_FT, ...
                nfix_FT, video_width, video_height, ...
                radius, gauss_template, sector_map, n_sectors, ...
                write_video, start_FT, dur_FT);
        end
    else
        for si = 1:nsbj
            [fix_bin_sub_density(:,:,si), subj_dur(si)] = single_subject_fix_density(...
                si, pairs_FT, sub_FT, ses_FT, dnfix_FT, xpos_FT, ypos_FT, ...
                nfix_FT, video_width, video_height, ...
                radius, gauss_template, sector_map, n_sectors, ...
                write_video, start_FT, dur_FT);
        end
    end
    
    % 时间统计输出
    totTime = toc(tAll);
    timeStr = char(datetime("now","Format","HH:mm:ss"));
    fprintf('[%s] All subjects done in %.2fs | mean per subj: %.2fs | median: %.2fs\n', ...
        timeStr, totTime, mean(subj_dur,'omitnan'), median(subj_dur,'omitnan'));
    
    % 去趋势处理并转置为 [nsbj x nfix x n_sectors]
    sub_fix_bin_density = zeros(nsbj, nfix_FT, n_sectors, 'single');
    for si = 1:nsbj
        C = squeeze(fix_bin_sub_density(:,:,si));
        if isempty(C), continue; end
        C_dt = C - 0.5 * (C(:, [end, 1:n_sectors-1]) + C(:, [2:n_sectors, 1]));
        sub_fix_bin_density(si,:,:) = single(C_dt);
    end
end