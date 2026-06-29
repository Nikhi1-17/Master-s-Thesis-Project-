%% =========================================================
%  PARAMETER RECOVERY STUDY — Surrogate Optimisation
%  omega0 FIXED at template value; 
%  Optimise: gamma, K_val and D
%
%  CHANGES FROM ORIGINAL:
%    - run_v2_kuramoto inner loop VECTORISED (O(N) mean-field)
%    - Progress bars via DataQueue (per-run + per-optimiser step)
%    - parfor-safe progress tracking with atomic counter
%    - Misc: cleaner variable scoping, comments updated
%% =========================================================

% the same code has been used to produce results in foldernrep10Results. I
% am using it here withnum_rep = 20

clear; clc; close all;

set(groot, 'defaultFigureRenderer', 'painters')

%% =========================================================
%  GLOBAL SETTINGS
%% =========================================================

tot_time  = 30;
num_osc   = 50;
num_rep   = 20;
mfe       = 500;       % Max function evaluations for surrogate optimiser
Fs        = 2000;
bw_bp     = 4;         % Bandpass half-width (Hz): omega0_t +/- 4 Hz
num_runs  = 32;        % Number of repetitions of the recovery study


%% =========================================================
%  opt_params STRUCT
%% =========================================================

opt_params.dt     = 1 / Fs;
opt_params.t_exp  = tot_time;
opt_params.N      = num_osc;
opt_params.nRep   = num_rep;
opt_params.Fs_exp = Fs;

%% =========================================================
%  STORAGE
%  Columns: [omega0, gamma, K_val, D, rho_avg]
%% =========================================================

template_params  = zeros(num_runs, 5);
recovered_params = zeros(num_runs, 5);

%% =========================================================
%  RESULTS FOLDER
%% =========================================================

results_folder = fullfile(pwd, 'param_recovery_results_omega0fixed');
if ~exist(results_folder, 'dir')
    mkdir(results_folder);
end

%% =========================================================
%  PARALLEL POOL
%% =========================================================

if isempty(gcp('nocreate'))
    parpool;
end
pool = gcp;
nWorkers = pool.NumWorkers;
fprintf('Workers available: %d\n', nWorkers);

%% =========================================================
%  PROGRESS TRACKING via DataQueue
%% =========================================================

dq_log  = parallel.pool.DataQueue;
dq_prog = parallel.pool.DataQueue;

afterEach(dq_log, @(msg) fprintf('%s\n', msg));

prog_bar = waitbar(0, sprintf('Runs completed: 0 / %d', num_runs), ...
                   'Name', 'Parameter Recovery Progress');

afterEach(dq_prog, @(~) update_runbar(prog_bar, num_runs));

%% =========================================================
%  TIMING
%% =========================================================

t_start = tic;
send(dq_log, sprintf('[%s]  Starting %d runs on %d workers ...', ...
    datestr(now,'HH:MM:SS'), num_runs, nWorkers));

%% =========================================================
%  MAIN PARFOR LOOP
%% =========================================================

parfor run_idx = 1:num_runs

    task = getCurrentTask();
    t_run = tic;

    send(dq_log, sprintf('[%s]  Run %02d started on worker %d', ...
        datestr(now,'HH:MM:SS'), run_idx, task.ID));

    rng(1000 + run_idx);   % reproducibility

    %% --------------------------------------------------
    %  STEP 1 — Draw random template parameters
    %% --------------------------------------------------

    omega0_t = 15 + (40 - 15) * rand();   % U(15, 40)
    gamma_t  =  1 + ( 8 -  1) * rand();   % U(1, 8)
    K_val_t  =  5 + (60 -  5) * rand();   % U(5, 60)
    D_t      =  0.01 + (4   - 0.01) * rand();   % U(0.01,  4)

    %% --------------------------------------------------
    %  STEP 2 — Simulate "experimental" LFP
    %% --------------------------------------------------

    [lfp_template, ~, ~, rho_avg_t] = run_kuramoto_fast( ...
        omega0_t, gamma_t, K_val_t, D_t, opt_params);

    template_params(run_idx, :) = [omega0_t, gamma_t, K_val_t, D_t, rho_avg_t]; %#ok<PFOUS>

    %% --------------------------------------------------
    %  STEP 3 — Bandpass filter
    %% --------------------------------------------------

    bp_signal = matlab_bpass_filt(lfp_template, ...
        omega0_t - bw_bp, omega0_t + bw_bp, Fs);

    %% --------------------------------------------------
    %  STEP 4 — Feature extraction
    %% --------------------------------------------------

    exp_features = compute_features(bp_signal, Fs);

    %% --------------------------------------------------
    %  STEP 5 — Surrogate Optimisation (gamma, K_val, D)
    %% --------------------------------------------------

    lb = [1,  5, 0.01];
    ub = [8, 60, 4];

    eval_counter = 0;
    log_interval = 50;

    % x = [gamma, K_val, D]
    % obj_fn = @(x) obj_wrapper(x, omega0_t, exp_features, ...
    %                           opt_params, run_idx, task.ID, ...
    %                           eval_counter, log_interval, dq_log, mfe);

    obj_fn = @(x) obj_wrapper(x, omega0_t, exp_features, ...
                          opt_params, run_idx, task.ID, ...
                          eval_counter, log_interval, dq_log, mfe);

    surr_opts = optimoptions('surrogateopt', ...
        'MaxFunctionEvaluations', mfe, ...
        'BatchUpdateInterval',    1, ...
        'Display',                'off', ...
        'UseParallel',            false, ...
        'PlotFcn',                []);

    [x_opt, fval_opt] = surrogateopt(obj_fn, lb, ub, surr_opts);

    gamma_r = x_opt(1);
    K_val_r = x_opt(2);
    D_r     = x_opt(3);

    %% --------------------------------------------------
    %  STEP 6 — Simulate recovered model → rho_avg_r
    %% --------------------------------------------------

    [lfp_recovered, ~, ~, rho_avg_r] = run_kuramoto_fast( ...
        omega0_t, gamma_r, K_val_r, D_r, opt_params);

    recovered_params(run_idx, :) = [omega0_t, gamma_r, K_val_r, D_r, rho_avg_r]; %#ok<PFOUS>

    %% --------------------------------------------------
    %  STEP 7 — Per-run comparison plot
    %% --------------------------------------------------

    base_name = sprintf('run%02d_omega%.1f_g%.2f_K%.2f_D%.3f', ...
                        run_idx, omega0_t, gamma_r, K_val_r, D_r);

    model_features = compute_features(lfp_recovered, Fs);

    plot_comparison(exp_features, model_features, omega0_t, x_opt, fval_opt, ...
    omega0_t, bw_bp * 2, base_name, results_folder, Fs);

    elapsed = toc(t_run);
    send(dq_log, sprintf('[%s]  Run %02d DONE (worker %d)  |  gamma=%.3f  K=%.3f  D=%.4f  NSS=%.5f  |  %.1f min', ...
        datestr(now,'HH:MM:SS'), run_idx, task.ID, gamma_r, K_val_r, D_r, fval_opt, elapsed/60));

    send(dq_prog, run_idx);

end   % end parfor

%% =========================================================
%  STEP 8 — Close progress bar, report total time
%% =========================================================

if ishandle(prog_bar); close(prog_bar); end

total_elapsed = toc(t_start);
fprintf('\n============================================================\n');
fprintf('All %d runs complete in %.1f minutes (%.1f hours)\n', ...
    num_runs, total_elapsed/60, total_elapsed/3600);
fprintf('============================================================\n\n');

%% =========================================================
%  STEP 9 — Save results tables
%% =========================================================

param_names = {'omega0', 'gamma', 'K_val', 'D', 'rho_avg'};
T_template  = array2table(template_params,  'VariableNames', param_names);
T_recovered = array2table(recovered_params, 'VariableNames', param_names);

writetable(T_template,  fullfile(results_folder, 'template_params.csv'));
writetable(T_recovered, fullfile(results_folder, 'recovered_params.csv'));
fprintf('Results saved to: %s\n', results_folder);

%% =========================================================
%  STEP 10 — Recovery summary plot
%% =========================================================

param_labels = {'\gamma (Hz)', 'K', 'D', '\rho_{avg}'};
param_ranges = {[1 8], [5 60], [0.01 4], [0 1]};
plot_cols    = [2, 3, 4, 5];
colors       = [0.15 0.45 0.75];

fig = figure('Color', 'w', 'Units', 'normalized', ...
             'Position', [0.08 0.10 0.85 0.40]);

sgtitle('Parameter Recovery — Surrogate Optimisation (\omega_0 fixed)', ...
    'FontSize', 13, 'FontWeight', 'bold');

for p = 1:4

    ax  = subplot(1, 4, p);
    col = plot_cols(p);

    t_vals = template_params(:, col);
    r_vals = recovered_params(:, col);
    lo = param_ranges{p}(1);
    hi = param_ranges{p}(2);

    line([lo hi], [lo hi], 'Color', [0.6 0.6 0.6], ...
        'LineStyle', '--', 'LineWidth', 1.2); hold on;

    scatter(t_vals, r_vals, 60, colors, 'filled', ...
        'MarkerFaceAlpha', 0.85, 'MarkerEdgeColor', 'w', 'LineWidth', 0.8);

    cf    = polyfit(t_vals, r_vals, 1);
    x_fit = linspace(lo, hi, 100);
    plot(x_fit, polyval(cf, x_fit), '-', 'Color', colors, 'LineWidth', 1.8);

    ss_res = sum((r_vals - polyval(cf, t_vals)).^2);
    ss_tot = sum((r_vals - mean(r_vals)).^2);
    r2     = 1 - ss_res / max(ss_tot, eps);

    xlabel(['Template  ', param_labels{p}],  'FontSize', 9);
    ylabel(['Recovered ', param_labels{p}],  'FontSize', 9);
    title(sprintf('%s\nslope=%.2f  R^2=%.2f', param_labels{p}, cf(1), r2), ...
        'FontSize', 9);

    xlim([lo hi]); ylim([lo hi]);
    axis square; grid on; box on;
    set(ax, 'FontSize', 8.5);

end

fig_path = fullfile(results_folder, 'recovery_plot_omega0fixed_gammaK_D_rhoavg.png');
set(fig, 'PaperPositionMode', 'auto');
print(fig, fig_path, '-dpng', '-r150');
fprintf('Recovery plot saved: %s\n', fig_path);


%% =========================================================
%  LOCAL FUNCTIONS
%% =========================================================

% Nested function for waitbar update (needs to be accessible)
function update_runbar(wb, total)
    % Called from afterEach — runs on client
    persistent count;
    if isempty(count); count = 0; end
    count = count + 1;
    if ishandle(wb)
        waitbar(count / total, wb, ...
            sprintf('Runs completed: %d / %d', count, total));
    end
end

% ---------------------------------------------------------
%  KURAMOTO — VECTORISED  (O(N) mean-field, not O(N²))
%
%  KEY CHANGE: The original code computed
%      dtheta   = theta' - theta;                  % N×N matrix
%      coupling = (K/N) * sum(sin(dtheta), 2);     % still O(N²)
%
%  Using the mean-field identity:
%      (1/N) Σ_j sin(θ_j − θ_i) = Im[ Z · exp(−iθ_i) ]
%  where Z = (1/N) Σ_j exp(iθ_j)  (order parameter)
%
%  This reduces the coupling step from O(N²) → O(N) per time step,
%  giving ~100–500x speedup for N=200, nSteps=60001, nRep=20.
% ---------------------------------------------------------
function [lfp_z, rho_mean, psi_mean, rho_avg] = run_kuramoto_fast( ...
        omega0, gamma, K_val, D, opt_params)

    dt          = opt_params.dt;
    t_exp       = opt_params.t_exp;
    N           = opt_params.N;
    nRep        = opt_params.nRep;

    nSteps      = round(t_exp / dt) + 1;
    noise_scale = sqrt(2 * D * dt);

    lfp_z   = zeros(1, nSteps);
    rho_all = zeros(nRep, nSteps);
    psi_all = zeros(nRep, nSteps);

    for rep = 1:nRep

        % Draw natural frequencies from Cauchy distribution
        omegas = 2 * pi * cauchy_rnd(omega0, gamma, 1, N);   % 1×N
        theta  = 2 * pi * rand(1, N);                         % 1×N  (initial phases)

        rho_rep = zeros(1, nSteps);
        psi_rep = zeros(1, nSteps);
        lfp_rep = zeros(1, nSteps);

        % t = 0
        Z           = mean(exp(1i * theta));
        lfp_rep(1)  = real(Z);
        rho_rep(1)  = abs(Z);
        psi_rep(1)  = angle(Z);

        % -------------------------------------------------
        %  VECTORISED EULER STEP
        %  coupling_i = (K/N) Σ_j sin(θ_j − θ_i)
        %             = K · Im[ Z · exp(−i·θ_i) ]
        %  No N×N matrix needed.
        % -------------------------------------------------
        for i = 1:nSteps - 1

            Z         = mean(exp(1i * theta));           % order parameter (scalar)
            coupling  = K_val * imag(Z * exp(-1i * theta));  % 1×N, O(N)

            theta = theta + (omegas + coupling) * dt ...
                          + noise_scale * randn(1, N);

            Z           = mean(exp(1i * theta));
            lfp_rep(i+1) = real(Z);
            rho_rep(i+1) = abs(Z);
            psi_rep(i+1) = angle(Z);

        end

        rho_all(rep, :) = rho_rep;
        psi_all(rep, :) = psi_rep;

        s = std(lfp_rep);
        if s > 0
            lfp_z = lfp_z + (lfp_rep - mean(lfp_rep)) / s;
        end

    end

    lfp_z    = lfp_z / nRep;
    rho_mean = mean(rho_all, 1);
    psi_mean = mean(psi_all, 1);

    steady_idx = round(0.2 * nSteps) : nSteps;
    rho_avg    = mean(rho_mean(steady_idx));
end

% ---------------------------------------------------------
%  OBJECTIVE FUNCTION WRAPPER
%  Adds lightweight per-evaluation logging via dq_log.
%  eval_counter is passed by value so it resets each run —
%  this is intentional (parfor slice semantics).
% ---------------------------------------------------------
function f = obj_wrapper(x, omega0_fixed, exp_features, ...
                          opt_params, run_idx, worker_id, ...
                          eval_counter, log_interval, dq_log, mfe)

    eval_counter = eval_counter + 1;  % local to this call

    gamma_val = x(1);
    K_val     = x(2);
    D_val     = x(3);

    sim_out = run_kuramoto_fast(omega0_fixed, gamma_val, K_val, D_val, opt_params);
    sim_out = double(sim_out(:)');

    Fs_model = 1 / opt_params.dt;
    Fs_exp   = opt_params.Fs_exp;

    if abs(Fs_model - Fs_exp) > 0.5
        [p, q]  = rat(Fs_exp / Fs_model, 1e-6);
        sim_out = resample(sim_out, p, q);
    end

    model_feat = compute_features(sim_out, Fs_exp);
    f          = error_function(exp_features, model_feat);

    % Log progress every log_interval evaluations
    if mod(eval_counter, log_interval) == 0
        send(dq_log, sprintf([ ...
            '  Run %02d | Worker %d | eval %3d/%d | NSS=%.5f | gamma=%.3f K=%.2f D=%.4f'], ...
            run_idx, worker_id, eval_counter, mfe, f, gamma_val, K_val, D_val));
    end
end

% ---------------------------------------------------------
%  BANDPASS FILTER  (zero-phase FIR, pure MATLAB)
% ---------------------------------------------------------
function bp_signal = matlab_bpass_filt(signal, low_cut, high_cut, Fs)

    signal   = signal(:);
    Wn       = [low_cut, high_cut] / (Fs / 2);
    bp_order = round(3 * Fs / low_cut);
    if mod(bp_order, 2) ~= 0; bp_order = bp_order + 1; end

    b = fir1(bp_order, Wn, 'bandpass', hamming(bp_order + 1));

    pad_samples = min(round(2 * Fs), floor((length(signal) - 1) / 2));
    signal_pad  = [flipud(signal(1:pad_samples)); signal; ...
                   flipud(signal(end - pad_samples + 1 : end))];

    bp_pad    = filtfilt(b, 1, signal_pad);
    bp_signal = bp_pad(pad_samples + 1 : end - pad_samples);
end

% ---------------------------------------------------------
%  FEATURE EXTRACTION
% ---------------------------------------------------------
function features = compute_features(signal, Fs)

    signal = double(signal(:)');
    s = std(signal);
    if s > 0; signal = (signal - mean(signal)) / s; end

    window   = Fs * 1;
    noverlap = window / 2;

    [pxx_sig, f_sig] = pwelch(signal, window, noverlap, [], Fs);

    analytic = hilbert(signal);
    envelope = abs(analytic);

    [pdf_counts, bin_edges] = histcounts(envelope, 100, 'Normalization', 'pdf');
    pdf_smooth  = smoothdata(pdf_counts, 'gaussian', 5);
    bin_centres = bin_edges(1:end-1) + diff(bin_edges) / 2;

    [pxx_env, f_env] = pwelch(envelope, window, noverlap, [], Fs);

    features.signal_psd      = pxx_sig;
    features.signal_freq     = f_sig;
    features.env_pdf         = pdf_smooth;
    features.env_pdf_bins    = bin_centres;
    features.env_psd         = pxx_env;
    features.env_psd_freq    = f_env;
    features.filtered_signal = signal;
    features.envelope        = envelope;
end

% ---------------------------------------------------------
%  ERROR FUNCTION
% ---------------------------------------------------------
function f = error_function(exp_feat, model_feat)
    f = (1/3) * (nss(exp_feat.signal_psd, model_feat.signal_psd) + ...
                 nss(exp_feat.env_pdf,     model_feat.env_pdf   ) + ...
                 nss(exp_feat.env_psd,     model_feat.env_psd   ));
end

% ---------------------------------------------------------
%  NASH-SUTCLIFFE EFFICIENCY
% ---------------------------------------------------------
function val = nss(d, m)
    n   = min(length(d), length(m));
    d   = d(1:n);
    m   = m(1:n);
    num = sum((d - m).^2);
    den = sum((d - mean(d)).^2);
    if den == 0
        val = double(num ~= 0) * 1e6;
    else
        val = num / den;
    end
end

% ---------------------------------------------------------
%  CAUCHY RANDOM SAMPLES
% ---------------------------------------------------------
function x = cauchy_rnd(x0, gamma, m, n)
    x = x0 + gamma * tan(pi * (rand(m, n) - 0.5));
end

% ---------------------------------------------------------
%  INTERPOLATE VECTOR TO TARGET LENGTH
% ---------------------------------------------------------
function y_out = interp_to_grid(y_in, n_out)
    y_in = y_in(:)';
    n_in = length(y_in);
    if n_in == n_out; y_out = y_in; return; end
    y_out = interp1(linspace(0,1,n_in), y_in, linspace(0,1,n_out), 'linear', 0);
end

% ---------------------------------------------------------
%  COMPARISON PLOT  (7-panel, per-run)
% ---------------------------------------------------------
function plot_comparison(exp_feat, feat2, omega0_fixed, x2, err2, ...
                          f0, bw, base_name, results_folder, Fs)

    c_exp = [0.12 0.47 0.71];
    c_2   = [1.00 0.50 0.05];
    c_err = [0.20 0.20 0.20];
    lw    = 2.5;

    fig = figure('Color','w','Units','normalized', ...
                 'Position',[0.01 0.01 0.98 0.95],'Visible','off');

    sgtitle(sprintf('Kuramoto Fit  |  %s  |  Fs: %g Hz  |  \\omega_0 = %.3f Hz (fixed)', ...
        strrep(base_name,'_',' '), Fs, omega0_fixed), ...
        'FontSize', 11, 'FontWeight', 'bold');

    san = @(v) max(real(v(:)), 1e-12);

    f_sig  = exp_feat.signal_freq(:);
    f_env  = exp_feat.env_psd_freq(:);
    psd_e  = san(exp_feat.signal_psd);
    psd2   = san(interp_to_grid(feat2.signal_psd, length(f_sig)));
    epsd_e = san(exp_feat.env_psd);
    epsd2  = san(interp_to_grid(feat2.env_psd, length(f_env)));
    nb     = min(length(exp_feat.env_pdf), length(feat2.env_pdf));
    bins_e = exp_feat.env_pdf_bins(1:nb)';
    pdf_e  = san(exp_feat.env_pdf(1:nb));
    pdf2   = san(feat2.env_pdf(1:nb));

    fok  = f_sig >= 0;
    feok = f_env >= 0;

    % --- Panel 1: Signal PSD ---
    subplot(2,4,1);
    plot(f_sig(fok), psd_e(fok), 'Color', c_exp, 'LineWidth', lw); hold on;
    plot(f_sig(fok), psd2(fok),  'Color', c_2,   'LineWidth', lw);
    xlabel('Frequency (Hz)'); ylabel('PSD');
    title('1. Signal PSD', 'FontSize', 9);
    legend({'Experimental','Model'}, 'Location','southwest','FontSize',8);
    xlim([f0-8, f0+8]); grid on; box on;

    % --- Panel 2: Envelope PDF ---
    subplot(2,4,2);
    plot(bins_e, pdf_e, 'Color', c_exp, 'LineWidth', lw); hold on;
    plot(bins_e, pdf2,  'Color', c_2,   'LineWidth', lw);
    xlabel('Envelope Amplitude'); ylabel('PDF');
    title('2. Envelope PDF', 'FontSize', 9);
    legend({'Experimental','Model'}, 'Location','northeast','FontSize',8);
    grid on; box on;

    % --- Panel 3: Envelope PSD ---
    subplot(2,4,3);
    plot(f_env(feok), epsd_e(feok), 'Color', c_exp, 'LineWidth', lw); hold on;
    plot(f_env(feok), epsd2(feok),  'Color', c_2,   'LineWidth', lw);
    xlabel('Frequency (Hz)'); ylabel('PSD');
    title('3. Envelope PSD', 'FontSize', 9);
    legend({'Experimental','Model'}, 'Location','southwest','FontSize',8);
    xlim([0 5]); grid on; box on;

    % --- Panel 4: Parameter table ---
    subplot(2,4,4); axis off;
    tab_str = sprintf( ...
        ['Optimised parameters\n' ...
         '  omega_0 = %.4f Hz  (fixed)\n\n' ...
         '  gamma   = %.4f Hz  (optimised)\n\n' ...
         '  K       = %.3f     (optimised)\n\n' ...
         '  D       = %.4f     (optimised)\n\n' ...
         'Bandpass BW = %.1f Hz\n\n' ...
         'NSS error = %.5f'], ...
        omega0_fixed, x2(1), x2(2), x2(3), bw, err2);

    text(0.05, 0.97, tab_str, ...
        'Units','normalized','VerticalAlignment','top', ...
        'FontSize',8.5,'FontName','FixedWidth', ...
        'BackgroundColor',[0.96 0.96 0.96], ...
        'EdgeColor',[0.5 0.5 0.5],'LineWidth',1.2,'Interpreter','none');
    title('4. Best-Fit Parameters', 'FontSize', 9);

    % --- Panel 5: Signal PSD + error ---
    subplot(2,4,5);
    plot(f_sig(fok), psd_e(fok), 'Color', c_exp, 'LineWidth', lw); hold on;
    plot(f_sig(fok), psd2(fok),  'Color', c_2,   'LineWidth', lw);
    plot(f_sig(fok), san(abs(psd_e(fok)-psd2(fok))), '--', 'Color', c_err, 'LineWidth', 1.5);
    xlabel('Frequency (Hz)'); ylabel('PSD');
    title('5. Signal PSD Pointwise', 'FontSize', 9);
    legend({'Exp','Model','|Error|'}, 'Location','southwest','FontSize',8);
    xlim([f0-8, f0+8]); grid on; box on;

    % --- Panel 6: Envelope PDF + error ---
    subplot(2,4,6);
    plot(bins_e, pdf_e, 'Color', c_exp, 'LineWidth', lw); hold on;
    plot(bins_e, pdf2,  'Color', c_2,   'LineWidth', lw);
    plot(bins_e, san(abs(pdf_e-pdf2)), '--', 'Color', c_err, 'LineWidth', 1.5);
    xlabel('Envelope Amplitude'); ylabel('PDF');
    title('6. Envelope PDF Pointwise', 'FontSize', 9);
    legend({'Exp','Model','|Error|'}, 'Location','northeast','FontSize',8);
    grid on; box on;

    % --- Panel 7: Envelope PSD + error ---
    subplot(2,4,7);
    plot(f_env(feok), epsd_e(feok), 'Color', c_exp, 'LineWidth', lw); hold on;
    plot(f_env(feok), epsd2(feok),  'Color', c_2,   'LineWidth', lw);
    plot(f_env(feok), san(abs(epsd_e(feok)-epsd2(feok))), '--', 'Color', c_err, 'LineWidth', 1.5);
    xlabel('Frequency (Hz)'); ylabel('PSD');
    title('7. Envelope PSD Pointwise', 'FontSize', 9);
    legend({'Exp','Model','|Error|'}, 'Location','southwest','FontSize',8);
    xlim([0 5]); grid on; box on;

    % --- Save ---
    fig_path = fullfile(results_folder, [base_name '_comparison.fig']);
    set(fig, 'Visible', 'on');
    savefig(fig, fig_path);
    close(fig);
end
