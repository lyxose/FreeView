function [sector_density, dur] = single_subject_fix_density(...
    % 辅助函数：处理单个被试（按注视序号）
    si, pairs_FT, sub_FT, ses_FT, dnfix_FT, xpos_FT, ypos_FT, ...
    nfix_FT, video_width, video_height, ...
    radius, gauss_template, sector_map, n_sectors, ...
    write_video, start_FT, dur_FT)
    
    tSub = tic;
    timeStr = char(datetime("now","Format","HH:mm:ss"));
    fprintf('[%s] Start subj %d (Sub %d, Ses %d)\n', ...
        timeStr, si, pairs_FT(si,1), pairs_FT(si,2));
    
    sector_density = zeros(nfix_FT, n_sectors, 'single');
    
    % 筛选当前被试数据
    m_sub = (sub_FT == pairs_FT(si,1)) & (ses_FT == pairs_FT(si,2));
    if ~any(m_sub)
        dur = toc(tSub);
        fprintf('[%s] Skip subj %d (no data) in %.2fs\n', timeStr, si, dur);
        return;
    end
    
    xS = single(xpos_FT(m_sub));
    yS = single(ypos_FT(m_sub));
    dnS = dnfix_FT(m_sub);
    
    % 可选：视频帧缓冲
    if write_video && nargin >= 17
        sS = single(start_FT(m_sub));
        dS = single(dur_FT(m_sub));
        video_frames = zeros(video_height, video_width, nfix_FT, 'single');
    else
        write_video = false;
    end
    
    % 逐注视序号处理
    for fi = 1:nfix_FT
        active = (dnS == fi);
        if ~any(active)
            continue;
        end
        
        xa = round(xS(active));
        ya = round(yS(active));
        
        valid = xa >= 1 & xa <= video_width & ya >= 1 & ya <= video_height;
        xa = xa(valid); ya = ya(valid);
        if isempty(xa), continue; end
        
        % 构建密度图
        frame_density = zeros(video_height, video_width, 'single');
        for jj = 1:numel(xa)
            xc = xa(jj); yc = ya(jj);
            x_min = max(1, xc - radius); x_max = min(video_width, xc + radius);
            y_min = max(1, yc - radius); y_max = min(video_height, yc + radius);
            
            tx_min = 1 + (x_min - (xc - radius));
            ty_min = 1 + (y_min - (yc - radius));
            tx_max = tx_min + (x_max - x_min);
            ty_max = ty_min + (y_max - y_min);
            
            frame_density(y_min:y_max, x_min:x_max) = ...
                frame_density(y_min:y_max, x_min:x_max) + ...
                gauss_template(ty_min:ty_max, tx_min:tx_max);
        end
        
        % 归一化
        ssum = sum(frame_density(:), 'native');
        if ssum > 0
            frame_density = frame_density / single(ssum);
        end
        
        if write_video
            video_frames(:,:,fi) = frame_density;
        end
        
        % 计算扇区统计
        sc = accumarray(double(sector_map(:)), double(frame_density(:)), [n_sectors 1], @sum, 0);
        sector_density(fi, :) = single(sc);
    end
    
    % 保存视频（可选）
    if write_video
        save_subject_fixation_density_video(video_frames, pairs_FT(si,:), nfix_FT);
    end
    
    dur = toc(tSub);
    timeStr = char(datetime("now","Format","HH:mm:ss"));
    fprintf('[%s] Done subj %d in %.2fs\n', timeStr, si, dur);
end

