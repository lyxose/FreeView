function spec_result = analyze_subject_spectrum(seriesAxis)
% 对输入的 nsbj x n_trials 时序数据做频谱分析和置换检验
% 返回结构体 spec_result，包含所有分析结果和可直接用于绘图的字段

    [n_subj, n_trials] = size(seriesAxis);

    n_fft = 2^nextpow2(n_trials * 4);              % zero-pad for higher resolution
    f = (0:floor(n_fft/2)) / n_fft;                % frequency axis (cycles per trial)
    n_freqs = numel(f);
    hann_window = hann(n_trials)';

    % subject-level spectra
    subj_spectra = zeros(n_subj, n_freqs);
    for s = 1:n_subj
        subj_data = seriesAxis(s,:);
        detrended = detrend(subj_data);
        sd_val = std(detrended);
        if sd_val < eps, sd_val = 1; end
        normalized = detrended / sd_val;
        windowed = normalized .* hann_window;
        Y = fft(windowed, n_fft);
        P2 = abs(Y / n_trials).^2;
        P1 = P2(1:floor(n_fft/2)+1);
        if n_freqs > 2
            P1(2:end-1) = 2 * P1(2:end-1);
        end
        subj_spectra(s,:) = P1;
    end

    mean_spectrum = mean(subj_spectra, 1);
    se_spectrum = std(subj_spectra, 0, 1) / sqrt(n_subj);
    se_spectrum(se_spectrum == 0) = eps;

    smooth_win = 5;
    mean_spectrum_smooth = movmean(mean_spectrum, smooth_win);
    se_spectrum_smooth = movmean(se_spectrum, smooth_win);

    % observed t-statistic (one-sample)
    t_obs = mean_spectrum ./ se_spectrum;
    valid_mask = true(1, n_freqs);
    valid_mask(1) = false; % drop DC

    % mean-then-FFT spectrum
    mean_timeseries = mean(seriesAxis, 1);
    mt_detrended = detrend(mean_timeseries);
    mt_sd = std(mt_detrended);
    if mt_sd < eps, mt_sd = 1; end
    mt_windowed = (mt_detrended / mt_sd) .* hann_window;
    Y_mean = fft(mt_windowed, n_fft);
    P2_mean = abs(Y_mean / n_trials).^2;
    mean_fft_spectrum = P2_mean(1:floor(n_fft/2)+1);
    if n_freqs > 2
        mean_fft_spectrum(2:end-1) = 2 * mean_fft_spectrum(2:end-1);
    end
    mean_fft_spectrum_smooth = movmean(mean_fft_spectrum, smooth_win);

    % permutation test (time shuffle per subject)
    n_perms = 5000;
    alpha = 0.05;
    cluster_threshold = 0.05;
    t_threshold = tinv(1 - cluster_threshold, n_subj - 1);

    perm_max_cluster = zeros(n_perms,1);
    for p = 1:n_perms
        perm_spectra = zeros(n_subj, n_freqs);
        for s = 1:n_subj
            permuted = seriesAxis(s, randperm(n_trials));
            perm_detr = detrend(permuted);
            perm_sd = std(perm_detr);
            if perm_sd < eps, perm_sd = 1; end
            perm_windowed = (perm_detr / perm_sd) .* hann_window;
            Yp = fft(perm_windowed, n_fft);
            P2p = abs(Yp / n_trials).^2;
            P1p = P2p(1:floor(n_fft/2)+1);
            if n_freqs > 2
                P1p(2:end-1) = 2 * P1p(2:end-1);
            end
            perm_spectra(s,:) = P1p;
        end
        mean_perm = mean(perm_spectra, 1);
        se_perm = std(perm_spectra, 0, 1) / sqrt(n_subj);
        se_perm(se_perm == 0) = eps;
        t_perm = mean_perm ./ se_perm;
        t_perm(~isfinite(t_perm)) = 0;
        t_perm(~valid_mask) = 0;

        mask_perm = (t_perm > t_threshold) & valid_mask;
        if any(mask_perm)
            [p_starts, p_ends] = find_runs(mask_perm);
            masses = zeros(numel(p_starts),1);
            for c = 1:numel(p_starts)
                masses(c) = sum(t_perm(p_starts(c):p_ends(c)));
            end
            perm_max_cluster(p) = max(masses);
        else
            perm_max_cluster(p) = 0;
        end
    end

    % observed significant clusters
    obs_mask = (t_obs > t_threshold) & valid_mask;
    [obs_starts, obs_ends] = find_runs(obs_mask);
    cluster_threshold_value = prctile(perm_max_cluster, 100 * (1 - alpha));

    sig_clusters = [];
    if ~isempty(obs_starts)
        for c = 1:numel(obs_starts)
            mass = sum(t_obs(obs_starts(c):obs_ends(c)));
            p_cluster = (sum(perm_max_cluster >= mass) + 1) / (n_perms + 1);
            if mass > cluster_threshold_value
                sig_clusters = [sig_clusters; obs_starts(c), obs_ends(c), mass, p_cluster]; %#ok<AGROW>
            end
        end
        if ~isempty(sig_clusters)
            sig_clusters = sortrows(sig_clusters,1);
        end
    end

    % peak detection (subject-level mean)
    [f_pks, f_locs] = findpeaks(mean_spectrum_smooth(2:end), 'SortStr','descend', ...
        'NPeaks', 5, 'MinPeakDistance', 3);
    peak_freqs = f(f_locs + 1);
    peak_periods = 1 ./ peak_freqs;

    % peak detection (mean-then-FFT)
    [f_pks_mt, f_locs_mt] = findpeaks(mean_fft_spectrum_smooth(2:end), 'SortStr','descend', ...
        'NPeaks', 5, 'MinPeakDistance', 3);
    peak_freqs_mt = f(f_locs_mt + 1);
    peak_periods_mt = 1 ./ peak_freqs_mt;

    % 汇总结果
    spec_result = struct();
    spec_result.n_subj = n_subj;
    spec_result.n_trials = n_trials;
    spec_result.f = f;
    spec_result.mean_spectrum = mean_spectrum;
    spec_result.se_spectrum = se_spectrum;
    spec_result.mean_spectrum_smooth = mean_spectrum_smooth;
    spec_result.se_spectrum_smooth = se_spectrum_smooth;
    spec_result.mean_fft_spectrum = mean_fft_spectrum;
    spec_result.mean_fft_spectrum_smooth = mean_fft_spectrum_smooth;
    spec_result.t_obs = t_obs;
    spec_result.valid_mask = valid_mask;
    spec_result.sig_clusters = sig_clusters;
    spec_result.peak_freqs = peak_freqs;
    spec_result.peak_periods = peak_periods;
    spec_result.f_pks = f_pks;
    spec_result.peak_freqs_mt = peak_freqs_mt;
    spec_result.peak_periods_mt = peak_periods_mt;
    spec_result.f_pks_mt = f_pks_mt;
    spec_result.mean_timeseries_group = mean(seriesAxis, 1);
    spec_result.se_timeseries = std(seriesAxis, 0, 1) / sqrt(n_subj);
    spec_result.alpha = alpha;
end

