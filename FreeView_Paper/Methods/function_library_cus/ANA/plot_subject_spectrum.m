function plot_subject_spectrum(spec_result, AxisColor)
% 绘制频谱分析结果，支持自定义主色 AxisColor

    if nargin < 2 || isempty(AxisColor)
        light_yellow = [245, 229, 38] / 255;
    else
        light_yellow = AxisColor;
    end
    deep_yellow = min(1, light_yellow * 0.85 + 0.15);
    sig_color = [0.82, 0.24, 0.24];

    % 取结果
    f = spec_result.f;
    n_subj = spec_result.n_subj;
    n_trials = spec_result.n_trials;
    mean_spectrum_smooth = spec_result.mean_spectrum_smooth;
    se_spectrum_smooth = spec_result.se_spectrum_smooth;
    mean_fft_spectrum_smooth = spec_result.mean_fft_spectrum_smooth;
    peak_freqs = spec_result.peak_freqs;
    peak_periods = spec_result.peak_periods;
    f_pks = spec_result.f_pks;
    peak_freqs_mt = spec_result.peak_freqs_mt;
    peak_periods_mt = spec_result.peak_periods_mt;
    f_pks_mt = spec_result.f_pks_mt;
    sig_clusters = spec_result.sig_clusters;
    mean_timeseries_group = spec_result.mean_timeseries_group;
    se_timeseries = spec_result.se_timeseries;
    alpha = spec_result.alpha;

    figure('Position',[120, 80, 760, 810]);

    % 时域
    subplot(3,1,1);
    t_axis = 1:n_trials;
    h_fill = fill([t_axis, fliplr(t_axis)], ...
        [mean_timeseries_group + se_timeseries, fliplr(mean_timeseries_group - se_timeseries)], ...
        deep_yellow, 'EdgeColor','none', 'FaceAlpha',0.2);
    set(get(get(h_fill,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    hold on;
    plot(t_axis, mean_timeseries_group, 'Color', deep_yellow, 'LineWidth', 1.6);
    title('Time Series (Mean \pm SE Across Subjects)');
    xlabel('Trial Number');
    ylabel('Value');
    grid on;
    hold off;

    % 频域
    subplot(3,1,2);
    freq_axis = f(2:end);
    main_curve = mean_spectrum_smooth(2:end);
    main_se = se_spectrum_smooth(2:end);
    mt_curve = mean_fft_spectrum_smooth(2:end);

    h_fill = fill([freq_axis, fliplr(freq_axis)], ...
        [main_curve + main_se, fliplr(main_curve - main_se)], ...
        light_yellow, 'EdgeColor','none', 'FaceAlpha',0.2);
    set(get(get(h_fill,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    hold on;
    h_main = plot(freq_axis, main_curve, '-', 'Color', light_yellow, 'LineWidth', 1.6);
    h_mt = plot(freq_axis, mt_curve, '--', 'Color', deep_yellow, 'LineWidth', 1.6);
    set(gca,'XScale','log');
    xlabel('Frequency (cycles per trial)');
    ylabel('Power');
    title('Power Spectrum (Subject Mean vs. Mean-then-FFT)');
    grid on;
    xlim([1e-3, 1e-1]);

    for i = 1:min(3, numel(peak_freqs))
        plot(peak_freqs(i), f_pks(i), 'o', 'MarkerSize',5, ...
            'MarkerFaceColor', light_yellow, 'MarkerEdgeColor', [0.2 0.2 0.2]);
        text(peak_freqs(i)*1.05, f_pks(i)*1.03, sprintf('T=%.1f', peak_periods(i)), ...
            'Color', light_yellow*0.7, 'FontSize',8, 'HorizontalAlignment','left', 'VerticalAlignment','bottom');
    end
    for i = 1:min(3, numel(peak_freqs_mt))
        plot(peak_freqs_mt(i), f_pks_mt(i), 'd', 'MarkerSize',5, ...
            'MarkerFaceColor', deep_yellow, 'MarkerEdgeColor', [0.2 0.2 0.2]);
        text(peak_freqs_mt(i)*0.98, f_pks_mt(i)*0.94, sprintf('T=%.1f', peak_periods_mt(i)), ...
            'Color', deep_yellow*0.8, 'FontSize',8, 'HorizontalAlignment','right', 'VerticalAlignment','top');
    end

    if ~isempty(sig_clusters)
        yl = ylim;
        y_sig = yl(2) - 0.06*(yl(2)-yl(1));
        for c = 1:size(sig_clusters,1)
            fx0 = f(sig_clusters(c,1));
            fx1 = f(sig_clusters(c,2));
            plot([fx0, fx1], [y_sig, y_sig], '-', 'Color', sig_color, 'LineWidth', 3);
        end
    end
    legend([h_main, h_mt], {'Subject-level mean', 'Mean-then-FFT'}, 'Location','northeast');
    hold off;

    % 周期域
    subplot(3,1,3);
    periods = 1 ./ freq_axis;
    power_vals = main_curve;
    power_se = main_se;
    power_mt = mt_curve;

    min_period = 2;
    max_needed = max([peak_periods(:); peak_periods_mt(:); 96]);
    if isempty(max_needed) || ~isfinite(max_needed)
        max_needed = 200;
    end
    max_period = min(max(periods), max(ceil(max_needed*1.1), 120));
    period_mask = periods >= min_period & periods <= max_period;

    x_period = periods(period_mask);
    y_period = power_vals(period_mask);
    se_period = power_se(period_mask);
    y_period_mt = power_mt(period_mask);

    [x_period, idx_sort] = sort(x_period);
    y_period = y_period(idx_sort);
    se_period = se_period(idx_sort);
    y_period_mt = y_period_mt(idx_sort);

    h_fill = fill([x_period, fliplr(x_period)], ...
        [y_period + se_period, fliplr(y_period - se_period)], ...
        light_yellow, 'EdgeColor','none', 'FaceAlpha',0.2);
    set(get(get(h_fill,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    hold on;
    h_main = plot(x_period, y_period, '-', 'Color', light_yellow, 'LineWidth', 1.6);
    h_mt = plot(x_period, y_period_mt, '--', 'Color', deep_yellow, 'LineWidth', 1.6);
    xlabel('Period (trials per cycle)');
    ylabel('Power');
    title('Power Spectrum vs. Period');
    grid on;

    for i = 1:min(3, numel(peak_periods))
        if peak_periods(i) >= min_period && peak_periods(i) <= max_period
            [~, idx] = min(abs(x_period - peak_periods(i)));
            plot(peak_periods(i), y_period(idx), 'o', 'MarkerSize',5, ...
                'MarkerFaceColor', light_yellow, 'MarkerEdgeColor',[0.2 0.2 0.2]);
            text(peak_periods(i)*1.02, y_period(idx)*1.03, sprintf('T=%.1f', peak_periods(i)), ...
                'Color', light_yellow*0.7, 'FontSize',8, 'HorizontalAlignment','left','VerticalAlignment','bottom');
        end
    end
    for i = 1:min(3, numel(peak_periods_mt))
        if peak_periods_mt(i) >= min_period && peak_periods_mt(i) <= max_period
            [~, idx] = min(abs(x_period - peak_periods_mt(i)));
            plot(peak_periods_mt(i), y_period_mt(idx), 'd', 'MarkerSize',5, ...
                'MarkerFaceColor', deep_yellow, 'MarkerEdgeColor',[0.2 0.2 0.2]);
            text(peak_periods_mt(i)*0.98, y_period_mt(idx)*0.96, sprintf('T=%.1f', peak_periods_mt(i)), ...
                'Color', deep_yellow*0.8, 'FontSize',8, 'HorizontalAlignment','right','VerticalAlignment','top');
        end
    end

    xline(96, '--', 'Color', [0.1 0.3 0.8], 'LineWidth', 1.5);
    ylp = ylim;
    text(96 + 2, ylp(2) - 0.05*(ylp(2)-ylp(1)), '96 (Block Size)', 'Color',[0.1 0.3 0.8], ...
         'FontSize',8, 'HorizontalAlignment','left', 'VerticalAlignment','top');

    xlim([0 min(350, max(x_period))]);

    if ~isempty(sig_clusters)
        y_sig = ylp(2) - 0.06*(ylp(2)-ylp(1));
        for c = 1:size(sig_clusters,1)
            fx0 = f(sig_clusters(c,1));
            fx1 = f(sig_clusters(c,2));
            px = sort([1 / fx0, 1 / fx1]);
            px(px > 350) = 350;
            plot(px, [y_sig, y_sig], '-', 'Color', sig_color, 'LineWidth', 3);
        end
    end
    legend([h_main, h_mt], {'Subject-level mean', 'Mean-then-FFT'}, 'Location','northeast');
    hold off;

    % 控制台输出
    fprintf('Power spectrum analysis with subject-level statistics (n=%d):\n', n_subj);
    fprintf('Top periods found in the group-level data:\n');
    for i = 1:min(3, numel(peak_periods))
        fprintf('Peak %d: Frequency = %.4f cycles/trial, Period = %.1f trials/cycle\n', ...
            i, peak_freqs(i), peak_periods(i));
    end

    if ~isempty(peak_periods_mt)
        fprintf('\nTop periods from mean-then-FFT:\n');
        for i = 1:min(3, numel(peak_periods_mt))
            fprintf('Peak %d: Frequency = %.4f cycles/trial, Period = %.1f trials/cycle\n', ...
                i, peak_freqs_mt(i), peak_periods_mt(i));
        end
    end

    fprintf('\nSignificant spectral clusters (p < %.3f):\n', alpha);
    if isempty(sig_clusters)
        fprintf('No significant clusters found\n');
    else
        for c = 1:size(sig_clusters,1)
            freq_start = f(sig_clusters(c,1));
            freq_end   = f(sig_clusters(c,2));
            period_hi  = 1 / max(freq_start, realmin);
            period_lo  = 1 / max(freq_end,   realmin);
            fprintf('Cluster %d: Frequencies %.4f-%.4f cycles/trial (periods %.1f-%.1f trials), p=%.3f\n', ...
                c, freq_start, freq_end, period_lo, period_hi, sig_clusters(c,4));
        end
    end
end