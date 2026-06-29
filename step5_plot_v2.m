function fig = step5_plot_v2(rest, pert, cfg, pp, row, save_fig)
% STEP5_PLOT_V2  Two-panel figure per patient.
%
%   Panel 1 (BEFORE): full 40 s resting rho(t).
%     - First 10 s (burn-in) shaded pink and labelled.
%     - Remaining 30 s (steady state) shaded green and labelled.
%
%   Panel 2 (AFTER): post-perturbation recovery, trimmed to the shortest
%     settled realization (all traces trimmed to the same length).
%     - 5-95 percentile band shaded.
%     - Mean recovery trace plotted on top.
%     - Baseline mean and +/- 2 SD band marked.
%     - Trim point (end of common window) marked with a vertical line.
%     - Individual settling times shown as rug ticks on the x-axis.
%
%   All rho(t) line traces use LineJoin = chamfer.
%   Saved as .fig (MATLAB figure format).

    if nargin < 6; save_fig = false; end

    fig = figure('Name', sprintf('rho before/after v2 - Patient %d', row), ...
                 'Position', [80 80 1300 540]);

    % ===================== Panel 1: BEFORE (40 s resting) =====================
    ax1 = subplot(1, 2, 1);
    t_pre = rest.t_full;
    hold(ax1, 'on');

    yl_pre = [0, max(rest.rho_full) * 1.12 + eps];
    xb     = rest.burnin_s;

    % burn-in shading (pink)
    patch(ax1, [0 xb xb 0], [yl_pre(1) yl_pre(1) yl_pre(2) yl_pre(2)], ...
          [0.95 0.85 0.85], 'EdgeColor', 'none', 'FaceAlpha', 0.6, ...
          'HandleVisibility', 'off');

    % steady-state shading (pale green)
    patch(ax1, [xb t_pre(end) t_pre(end) xb], ...
          [yl_pre(1) yl_pre(1) yl_pre(2) yl_pre(2)], ...
          [0.85 0.92 0.85], 'EdgeColor', 'none', 'FaceAlpha', 0.5, ...
          'HandleVisibility', 'off');

    % rho(t) trace
    hPre = plot(ax1, t_pre, rest.rho_full, 'b', 'LineWidth', 0.9);
    set(hPre, 'LineJoin', 'chamfer');

    % baseline mean line
    hM = yline(ax1, rest.rho_baseline, 'k:', 'LineWidth', 1.0);
    try; set(hM, 'LineJoin', 'chamfer'); catch; end

    % burn-in boundary
    hBx = xline(ax1, xb, 'Color', [0.5 0 0], 'LineStyle', '--', 'LineWidth', 1.2);
    try; set(hBx, 'LineJoin', 'chamfer'); catch; end

    % region labels (placed near top)
    text(ax1, xb/2,                  yl_pre(2)*0.96, 'burn-in (10 s)', ...
         'HorizontalAlignment', 'center', 'Color', [0.5 0 0], 'FontWeight', 'bold', 'FontSize', 9);
    text(ax1, (xb + t_pre(end))/2,   yl_pre(2)*0.96, 'steady state (30 s)', ...
         'HorizontalAlignment', 'center', 'Color', [0 0.4 0], 'FontWeight', 'bold', 'FontSize', 9);

    xlabel(ax1, 'Time (s)');
    ylabel(ax1, '\rho(t) = |Z(t)|');
    title(ax1, 'Before perturbation  (40 s resting run)');
    legend(ax1, {'\rho(t)', 'baseline \rho_*', 'burn-in end'}, 'Location', 'best');
    grid(ax1, 'on');
    ylim(ax1, yl_pre);
    xlim(ax1, [0, t_pre(end)]);

    % ===================== Panel 2: AFTER (recovery, trimmed) =================
    ax2 = subplot(1, 2, 2);
    t_post = pert.t_trim;                          % trimmed time axis (s)
    hold(ax2, 'on');

    band = cfg.settle_nstd * rest.rho_std;
    yl_post = [0, max([max(pert.p95), rest.rho_baseline + band]) * 1.12 + eps];

    % 5-95 percentile shading
    fill(ax2, [t_post, fliplr(t_post)], [pert.p05, fliplr(pert.p95)], ...
         [0.95 0.80 0.70], 'EdgeColor', 'none', 'FaceAlpha', 0.55, ...
         'DisplayName', '5–95 pct band');

    % mean recovery trace
    hMean = plot(ax2, t_post, pert.mean_trace, ...
                 'Color', [0.85 0.33 0.10], 'LineWidth', 1.5, ...
                 'DisplayName', 'mean \rho(t)');
    set(hMean, 'LineJoin', 'chamfer');

    % baseline mean
    hBL = yline(ax2, rest.rho_baseline, 'k:', 'LineWidth', 1.0, ...
                'DisplayName', 'baseline \rho_*');
    try; set(hBL, 'LineJoin', 'chamfer'); catch; end

    % +/- 2 SD band lines
    % hUp = yline(ax2, rest.rho_baseline + band, '--', ...
    %             'Color', [0.3 0.3 0.3], 'LineWidth', 0.9, ...
    %             'DisplayName', ['\pm' num2str(cfg.settle_nstd) ' SD']);
    % hLo = yline(ax2, max(0, rest.rho_baseline - band), '--', ...
    %             'Color', [0.3 0.3 0.3], 'LineWidth', 0.9, ...
    %             'HandleVisibility', 'off');
    % try; set([hUp hLo], 'LineJoin', 'chamfer'); catch; end

    % %%% Can delete if ya want
    % % +/- niqr*IQR band lines (median-centered)
    % hUp = yline(ax2, rho_med + band, '--', ...
    %             'Color', [0.3 0.3 0.3], 'LineWidth', 0.9, ...
    %             'DisplayName', [num2str(niqr) '\timesIQR']);
    % hLo = yline(ax2, max(0, rho_med - band), '--', ...
    %             'Color', [0.3 0.3 0.3], 'LineWidth', 0.9, ...
    %             'HandleVisibility', 'off');
    % try; set([hUp hLo], 'LineJoin', 'chamfer'); catch; end
    % %%% Can delete ...

    % trim-point marker (end of common window = shortest settled realization)
    t_trim_end = t_post(end);
    hTrim = xline(ax2, t_trim_end, 'Color', [0.4 0 0.6], ...
                  'LineStyle', ':', 'LineWidth', 1.2, ...
                  'DisplayName', 'trim point (shortest settled)');
    try; set(hTrim, 'LineJoin', 'chamfer'); catch; end
    text(ax2, t_trim_end, yl_post(2)*0.88, ...
         sprintf('  trim\n  %.1f s', t_trim_end), ...
         'Color', [0.4 0 0.6], 'FontSize', 8);

    % rug ticks: individual settling times for each of the 24 realizations
    valid_st = pert.settle_times(~isnan(pert.settle_times));
    if ~isempty(valid_st)
        rug_y = yl_post(1) + 0.02 * diff(yl_post);
        plot(ax2, valid_st, repmat(rug_y, size(valid_st)), ...
             '|', 'Color', [0 0.45 0.74], 'MarkerSize', 6, ...
             'DisplayName', 'individual settle times');
    end

    xlabel(ax2, 'Time since kick (s)');
    ylabel(ax2, '\rho(t) = |Z(t)|');
    c_val = pert.target / rest.rho_baseline;
    title(ax2, sprintf('After perturbation  (n = %d realizations,  c = %.2f)', ...
          cfg.n_realizations, c_val));
    legend(ax2, 'Location', 'best');
    grid(ax2, 'on');
    xlim(ax2, [0, t_trim_end * 1.02]);
    ylim(ax2, yl_post);

    % ===================== Super-title ========================================
    sgtitle(sprintf( ...
        ['Patient %d   |   \\omega_0 = %.2f Hz,  K = %.2f,  ' ...
         '\\gamma = %.2f,  D = %.4f\n' ...
         '\\rho_* = %.4f \\pm %.4f (2SD band = \\pm%.4f),   ' ...
         'K_c = %.2f,   K - K_c = %.2f\n' ...
         'trim length = %.2f s   |   ' ...
         'settle times: min = %.2f s, max = %.2f s, median = %.2f s'], ...
        row, pp.omega0, pp.K, pp.gamma, pp.D, ...
        rest.rho_baseline, rest.rho_std, band, ...
        4*pi*pp.gamma + 2*pp.D, pp.K - (4*pi*pp.gamma + 2*pp.D), ...
        t_trim_end, ...
        min(pert.settle_times), max(pert.settle_times), median(pert.settle_times)), ...
        'FontWeight', 'bold', 'FontSize', 8);

    % ===================== Save as .fig =======================================
    if save_fig
        fname = sprintf('rho_before_after_v2_patient_%02d.fig', row);
        savefig(fig, fname);
        fprintf('Saved  %s\n', fname);
    end
end

% function fig = step5_plot_v2(rest, pert, cfg, pp, row, save_fig)
% % STEP5_PLOT_V2  Two-panel figure per patient.
% %
% %   Panel 1 (BEFORE): full 40 s resting rho(t).
% %     - First 10 s (burn-in) shaded pink and labelled.
% %     - Remaining 30 s (steady state) shaded green and labelled.
% %
% %   Panel 2 (AFTER): post-perturbation recovery, trimmed to the shortest
% %     settled realization (all traces trimmed to the same length).
% %     - 5-95 percentile band shaded.
% %     - Mean recovery trace plotted on top.
% %     - Baseline median and median-based IQR fence marked.
% %     - Trim point (end of common window) marked with a vertical line.
% %     - Individual settling times shown as rug ticks on the x-axis.
% %
% %   All rho(t) line traces use LineJoin = chamfer.
% %   Saved as .fig (MATLAB figure format).
% 
%     if nargin < 6; save_fig = false; end
% 
%     fig = figure('Name', sprintf('rho before/after v2 - Patient %d', row), ...
%                  'Position', [80 80 1300 540]);
% 
%     % ===================== Panel 1: BEFORE (40 s resting) =====================
%     ax1 = subplot(1, 2, 1);
%     t_pre = rest.t_full;
%     hold(ax1, 'on');
% 
%     yl_pre = [0, max(rest.rho_full) * 1.12 + eps];
%     xb     = rest.burnin_s;
% 
%     % burn-in shading (pink)
%     patch(ax1, [0 xb xb 0], [yl_pre(1) yl_pre(1) yl_pre(2) yl_pre(2)], ...
%           [0.95 0.85 0.85], 'EdgeColor', 'none', 'FaceAlpha', 0.6, ...
%           'HandleVisibility', 'off');
% 
%     % steady-state shading (pale green)
%     patch(ax1, [xb t_pre(end) t_pre(end) xb], ...
%           [yl_pre(1) yl_pre(1) yl_pre(2) yl_pre(2)], ...
%           [0.85 0.92 0.85], 'EdgeColor', 'none', 'FaceAlpha', 0.5, ...
%           'HandleVisibility', 'off');
% 
%     % rho(t) trace
%     hPre = plot(ax1, t_pre, rest.rho_full, 'b', 'LineWidth', 0.9);
%     set(hPre, 'LineJoin', 'chamfer');
% 
%     % baseline mean line
%     hM = yline(ax1, rest.rho_baseline, 'k:', 'LineWidth', 1.0);
%     try; set(hM, 'LineJoin', 'chamfer'); catch; end
% 
%     % burn-in boundary
%     hBx = xline(ax1, xb, 'Color', [0.5 0 0], 'LineStyle', '--', 'LineWidth', 1.2);
%     try; set(hBx, 'LineJoin', 'chamfer'); catch; end
% 
%     % region labels (placed near top)
%     text(ax1, xb/2,                  yl_pre(2)*0.96, 'burn-in (10 s)', ...
%          'HorizontalAlignment', 'center', 'Color', [0.5 0 0], 'FontWeight', 'bold', 'FontSize', 9);
%     text(ax1, (xb + t_pre(end))/2,   yl_pre(2)*0.96, 'steady state (30 s)', ...
%          'HorizontalAlignment', 'center', 'Color', [0 0.4 0], 'FontWeight', 'bold', 'FontSize', 9);
% 
%     xlabel(ax1, 'Time (s)');
%     ylabel(ax1, '\rho(t) = |Z(t)|');
%     title(ax1, 'Before perturbation  (40 s resting run)');
%     legend(ax1, {'\rho(t)', 'baseline \rho_*', 'burn-in end'}, 'Location', 'best');
%     grid(ax1, 'on');
%     ylim(ax1, yl_pre);
%     xlim(ax1, [0, t_pre(end)]);
% 
%     % ===================== Panel 2: AFTER (recovery, trimmed) =================
%     ax2 = subplot(1, 2, 2);
%     t_post = pert.t_trim;                          % trimmed time axis (s)
%     hold(ax2, 'on');
% 
%     % --- median-based settling band (pulled from what step3 actually used) ---
%     rho_med = pert.settle_center;     % = rest.rho_median
%     band    = pert.settle_band;       % = niqr * rest.rho_iqr
%     niqr    = 0.75;                    % for the legend label only
%     if isfield(cfg, 'settle_niqr'); niqr = cfg.settle_niqr; end
% 
%     yl_post = [0, max([max(pert.p95), rho_med + band]) * 1.12 + eps];
% 
%     % 5-95 percentile shading
%     fill(ax2, [t_post, fliplr(t_post)], [pert.p05, fliplr(pert.p95)], ...
%          [0.95 0.80 0.70], 'EdgeColor', 'none', 'FaceAlpha', 0.55, ...
%          'DisplayName', '5–95 pct band');
% 
%     % mean recovery trace
%     hMean = plot(ax2, t_post, pert.mean_trace, ...
%                  'Color', [0.85 0.33 0.10], 'LineWidth', 1.5, ...
%                  'DisplayName', 'mean \rho(t)');
%     set(hMean, 'LineJoin', 'chamfer');
% 
%     % baseline median (band center)
%     hBL = yline(ax2, rho_med, 'k:', 'LineWidth', 1.0, ...
%                 'DisplayName', 'baseline median \rho_*');
%     try; set(hBL, 'LineJoin', 'chamfer'); catch; end
% 
%     % +/- niqr*IQR band lines (median-centered)
%     hUp = yline(ax2, rho_med + band, '--', ...
%                 'Color', [0.3 0.3 0.3], 'LineWidth', 0.9, ...
%                 'DisplayName', [num2str(niqr) '\timesIQR fence']);
%     hLo = yline(ax2, max(0, rho_med - band), '--', ...
%                 'Color', [0.3 0.3 0.3], 'LineWidth', 0.9, ...
%                 'HandleVisibility', 'off');
%     try; set([hUp hLo], 'LineJoin', 'chamfer'); catch; end
% 
%     % trim-point marker (end of common window = shortest settled realization)
%     t_trim_end = t_post(end);
%     hTrim = xline(ax2, t_trim_end, 'Color', [0.4 0 0.6], ...
%                   'LineStyle', ':', 'LineWidth', 1.2, ...
%                   'DisplayName', 'trim point (shortest settled)');
%     try; set(hTrim, 'LineJoin', 'chamfer'); catch; end
%     text(ax2, t_trim_end, yl_post(2)*0.88, ...
%          sprintf('  trim\n  %.1f s', t_trim_end), ...
%          'Color', [0.4 0 0.6], 'FontSize', 8);
% 
%     % rug ticks: individual settling times for each of the 24 realizations
%     valid_st = pert.settle_times(~isnan(pert.settle_times));
%     if ~isempty(valid_st)
%         rug_y = yl_post(1) + 0.02 * diff(yl_post);
%         plot(ax2, valid_st, repmat(rug_y, size(valid_st)), ...
%              '|', 'Color', [0 0.45 0.74], 'MarkerSize', 6, ...
%              'DisplayName', 'individual settle times');
%     end
% 
%     xlabel(ax2, 'Time since kick (s)');
%     ylabel(ax2, '\rho(t) = |Z(t)|');
%     c_val = pert.target / rest.rho_baseline;
%     title(ax2, sprintf('After perturbation  (n = %d realizations,  c = %.2f)', ...
%           cfg.n_realizations, c_val));
%     legend(ax2, 'Location', 'best');
%     grid(ax2, 'on');
%     xlim(ax2, [0, t_trim_end * 1.02]);
%     ylim(ax2, yl_post);
% 
%     % ===================== Super-title ========================================
%     sgtitle(sprintf( ...
%         ['Patient %d   |   \\omega_0 = %.2f Hz,  K = %.2f,  ' ...
%          '\\gamma = %.2f,  D = %.4f\n' ...
%          '\\rho_{med} = %.4f   (%.2f\\timesIQR fence = \\pm%.4f),   ' ...
%          'K_c = %.2f,   K - K_c = %.2f\n' ...
%          'trim length = %.2f s   |   ' ...
%          'settle times: min = %.2f s, max = %.2f s, median = %.2f s'], ...
%         row, pp.omega0, pp.K, pp.gamma, pp.D, ...
%         rho_med, niqr, band, ...
%         4*pi*pp.gamma + 2*pp.D, pp.K - (4*pi*pp.gamma + 2*pp.D), ...
%         t_trim_end, ...
%         min(pert.settle_times), max(pert.settle_times), median(pert.settle_times)), ...
%         'FontWeight', 'bold', 'FontSize', 8);
% 
%     % ===================== Save as .fig =======================================
%     if save_fig
%         fname = sprintf('rho_before_after_v2_patient_%02d.fig', row);
%         savefig(fig, fname);
%         fprintf('Saved  %s\n', fname);
%     end
% end