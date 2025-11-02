function save_subject_density_video(video_frames, pair, frame_time_res, n_frames)
% 辅助函数：保存视频
    try
        vname = sprintf('fixation_density_video_Sub%02d_Ses%02d.avi', pair(1), pair(2));
        vw = VideoWriter(vname);
        vw.FrameRate = 1000 / frame_time_res;
        open(vw);
        fh = figure('Visible','off'); ax = axes(fh);
        for fi = 1:n_frames
            imshow(video_frames(:,:,fi), [], 'Parent', ax);
            title(ax, sprintf('Sub %d Ses %d | Time: %.0f ms', pair(1), pair(2), (fi-1)*frame_time_res));
            frame = getframe(ax);
            writeVideo(vw, frame);
        end
        close(vw); close(fh);
    catch ME
        warning('Video writing failed for Sub %d Ses %d: %s', pair(1), pair(2), ME.message);
    end
end
