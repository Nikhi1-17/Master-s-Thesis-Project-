%% =========================================================
%  KURAMOTO FIT PIPELINE — REAL LFP DATA
%
%  One row per STN side from table `alfo_opt_ready`.
%
%  Per row:
%    1. Welch PSD of best_epoch.signal -> FOOOF over BETA_BAND.
%       omega0 = CF of the highest-power peak (fixed thereafter).
%    2. Bandpass best_epoch.signal at omega0 +/- bw_bp, z-score
%       -> stored in column bp_zscore.
%    3. exp_features = compute_features(bp_zscore, Fs).
%    4. surrogateopt on x = [gamma, K, D], omega0 fixed,
%       repeated n_opt_runs times (different seeds).
%    5. Best fit = minimum-NSS run.
%    6. Forward-simulate best fit -> nRep independent LFP realisations.
%    7. Per-rep features -> mean + 5th/95th percentile bands.
%    8. Comparison plot: experimental (one colour) vs model mean +
%       shaded 5-95th percentile band (second colour).
%
%  CSV output: one row per (side x opt_run):
%       Name, side, omega0, gamma, K, D, nss
%
%  NOTE: stochastic spread on the plot is the MODEL's spread at the
%        single best parameter set (nRep realisations) — NOT optimiser
%        jitter and NOT parameter uncertainty. The n_opt_runs fits and
%        their NSS values (in the CSV) are the inference-spread record.
%% =========================================================

clc; clearvars; close all;

set(groot, 'defaultFigureRenderer', 'painters');

load ("alfo_opt_ready.mat");
alfo_opt_ready = Results;
clear ("Results")

%% =========================================================
%  GLOBAL SETTINGS
%% =========================================================

BETA_BAND   = [12 33];     % Hz — FOOOF fit range (matches table build)
bw_bp       = 3;           % Hz — bandpass half-width around CF
n_opt_runs  = 10;          % surrogateopt restarts per row
nRep        = 10;          % forward realisations of the best fit
mfe         = 500;         % MaxFunctionEvaluations per surrogateopt run
Fs_model    = 2000;        % Hz — model integration / sampling rate
num_osc     = 200;         % N oscillators
tot_time    = 30;          % s — simulated duration
pct_band    = [5 95];      % percentile band for the model envelope

% FOOOF settings (as specified)
fooof_settings.peak_width_limits = [2 12];
fooof_settings.max_n_peaks       = 6;
fooof_settings.min_peak_height   = 0;
fooof_settings.peak_threshold    = 2;
fooof_settings.aperiodic_mode    = 'fixed';
fooof_settings.verbose           = false;

%% =========================================================
%  opt_params STRUCT (model side; Fs_exp set per-row)
%% =========================================================

opt_params.dt     = 1 / Fs_model;
opt_params.t_exp  = tot_time;
opt_params.N      = num_osc;
opt_params.nRep   = nRep;
opt_params.Fs_exp = NaN;          % filled per row from best_epoch.time

%% =========================================================
%  RESULTS FOLDER
%% =========================================================

results_folder = fullfile(pwd, 'kuramoto_fit_results');
if ~exist(results_folder, 'dir'); mkdir(results_folder); end

%% =========================================================
%  PULL TABLE PIECES (avoid broadcasting whole table into parfor)
%% =========================================================

num_runs   = height(alfo_opt_ready);
be_all     = alfo_opt_ready.best_epoch;     % cell of structs (.signal, .time)
names_all  = string(alfo_opt_ready.Name);
sides_all  = string(alfo_opt_ready.side);

% Sliced outputs collected inside parfor, assembled afterwards
fooof_col  = cell(num_runs, 1);    % -> new column FOOOF
bpz_col    = cell(num_runs, 1);    % -> new column bp_zscore
fit_blocks = cell(num_runs, 1);    % per row: n_opt_runs x 5 [omega0 gamma K D nss]

%% =========================================================
%  PARALLEL POOL
%% =========================================================

if isempty(gcp('nocreate')); parpool; end
pool     = gcp;
nWorkers = pool.NumWorkers;
fprintf('Workers available: %d\n', nWorkers);

%% =========================================================
%  PROGRESS TRACKING via DataQueue
%% =========================================================

dq_log  = parallel.pool.DataQueue;
dq_prog = parallel.pool.DataQueue;

afterEach(dq_log, @(msg) fprintf('%s\n', msg));

prog_bar = waitbar(0, sprintf('Rows completed: 0 / %d', num_runs), ...
                   'Name', 'Kuramoto Fit Progress');
afterEach(dq_prog, @(~) update_runbar(prog_bar, num_runs));

%% =========================================================
%  TIMING
%% =========================================================

t_start = tic;
send(dq_log, sprintf('[%s]  Starting %d rows on %d workers ...', ...
    datestr(now,'HH:MM:SS'), num_runs, nWorkers));

%% =========================================================
%  MAIN PARFOR LOOP
%% =========================================================

parfor run_idx = 1:num_runs

    task  = getCurrentTask();
    t_run = tic;

    send(dq_log, sprintf('[%s]  Row %02d started on worker %d', ...
        datestr(now,'HH:MM:SS'), run_idx, task.ID));

    % Local copy of opt_params (so Fs_exp can be set per row inside parfor)
    op = opt_params;

    % -------- pull this row's epoch --------
    ep      = be_all{run_idx};
    sig_exp = double(ep.signal(:)');
    t_exp   = double(ep.time(:)');
    Fs      = 1 / mean(diff(t_exp));        % experimental Fs from time vector
    op.Fs_exp = Fs;

    %% --------------------------------------------------
    %  STEP 1 — FOOOF, select omega0 (highest-power peak)
    %% --------------------------------------------------

    win      = round(Fs * 1);
    noverlap = round(win / 2);
    [pxx, fxx] = pwelch(sig_exp, win, noverlap, [], Fs);   % linear power

    fm = fooof(fxx(:)', pxx(:)', BETA_BAND, fooof_settings, true);
    fooof_col{run_idx} = fm;

    pk = fm.peak_params;     % N x 3 : [CF, PW, BW]
    if isempty(pk)
        % --- no peak fallback: max raw power within BETA_BAND ---
        in_band   = fxx >= BETA_BAND(1) & fxx <= BETA_BAND(2);
        f_in      = fxx(in_band);
        p_in      = pxx(in_band);
        [~, imax] = max(p_in);
        omega0    = f_in(imax);
        send(dq_log, sprintf(['  Row %02d | FOOOF found NO peak — ' ...
            'falling back to max-power freq omega0=%.3f Hz'], run_idx, omega0));
    else
        [~, ipk] = max(pk(:, 2));    % strongest peak by power above aperiodic
        omega0   = pk(ipk, 1);       % its centre frequency
    end

    %% --------------------------------------------------
    %  STEP 2 — Bandpass omega0 +/- bw_bp, z-score
    %% --------------------------------------------------

    bp = matlab_bpass_filt(sig_exp, omega0 - bw_bp, omega0 + bw_bp, Fs);
    bp = bp(:)';
    sdv = std(bp);
    if sdv > 0; bp = (bp - mean(bp)) / sdv; end
    bpz_col{run_idx} = bp;

    %% --------------------------------------------------
    %  STEP 3 — Experimental features (single reference)
    %% --------------------------------------------------

    exp_features = compute_features(bp, Fs);

    %% --------------------------------------------------
    %  STEP 4 — surrogateopt x = [gamma, K, D], n_opt_runs times
    %% --------------------------------------------------

    lb = [1,  5, 0.01];
    ub = [8, 60, 4];

    fits = zeros(n_opt_runs, 4);   % [gamma K D nss]

    surr_opts = optimoptions('surrogateopt', ...
        'MaxFunctionEvaluations', mfe, ...
        'BatchUpdateInterval',    1, ...
        'Display',                'off', ...
        'UseParallel',            false, ...
        'PlotFcn',                []);

    for j = 1:n_opt_runs
        rng(1000 * run_idx + j);    % reproducible & distinct per restart

        obj_fn = @(x) obj_wrapper(x, omega0, exp_features, op, ...
                                  Fs_model, bw_bp);

        [x_opt, fval_opt] = surrogateopt(obj_fn, lb, ub, surr_opts);

        fits(j, :) = [x_opt(:)', fval_opt];

        send(dq_log, sprintf(['  Row %02d | opt %02d/%d | NSS=%.5f | ' ...
            'gamma=%.3f K=%.2f D=%.4f'], run_idx, j, n_opt_runs, ...
            fval_opt, x_opt(1), x_opt(2), x_opt(3)));
    end

    % store full block [omega0 gamma K D nss] for the CSV
    fit_blocks{run_idx} = [repmat(omega0, n_opt_runs, 1), fits];

    %% --------------------------------------------------
    %  STEP 5 — Best fit = minimum NSS
    %% --------------------------------------------------

    [fval_best, jbest] = min(fits(:, 4));
    x_best  = fits(jbest, 1:3);     % [gamma K D]
    gamma_r = x_best(1);
    K_r     = x_best(2);
    D_r     = x_best(3);

    %% --------------------------------------------------
    %  STEP 6 — Forward-simulate best fit, retain per-rep LFPs
    %% --------------------------------------------------

    lfp_all = run_kuramoto_fast(omega0, gamma_r, K_r, D_r, op);  % nRep x nSteps

    %% --------------------------------------------------
    %  STEP 7 — Per-rep features -> mean + percentile bands
    %% --------------------------------------------------

    model_bands = build_model_bands(lfp_all, Fs_model, Fs, ...
                                    omega0, bw_bp, exp_features, pct_band);

    %% --------------------------------------------------
    %  STEP 8 — Comparison plot
    %% --------------------------------------------------

    base_name = sprintf('row%02d_%s_%s_w%.1f_g%.2f_K%.2f_D%.3f', ...
        run_idx, sanitise(names_all(run_idx)), sides_all(run_idx), ...
        omega0, gamma_r, K_r, D_r);

    plot_comparison(exp_features, model_bands, omega0, x_best, ...
        fval_best, bw_bp, pct_band, base_name, results_folder);

    elapsed = toc(t_run);
    send(dq_log, sprintf(['[%s]  Row %02d DONE (worker %d) | ' ...
        'omega0=%.3f gamma=%.3f K=%.3f D=%.4f NSS=%.5f | %.1f min'], ...
        datestr(now,'HH:MM:SS'), run_idx, task.ID, ...
        omega0, gamma_r, K_r, D_r, fval_best, elapsed/60));

    send(dq_prog, run_idx);

end   % end parfor

%% =========================================================
%  STEP 9 — Close progress bar, report total time
%% =========================================================

if ishandle(prog_bar); close(prog_bar); end

total_elapsed = toc(t_start);
fprintf('\n============================================================\n');
fprintf('All %d rows complete in %.1f minutes (%.1f hours)\n', ...
    num_runs, total_elapsed/60, total_elapsed/3600);
fprintf('============================================================\n\n');

%% =========================================================
%  STEP 10 — Attach new columns, save table + CSV
%% =========================================================

alfo_opt_ready.FOOOF     = fooof_col;
alfo_opt_ready.bp_zscore = bpz_col;

save(fullfile(results_folder, 'alfo_opt_ready_fitted.mat'), 'alfo_opt_ready');

% Build CSV: one row per (side x opt_run)
Name_c = strings(0,1); side_c = strings(0,1);
omega_c = []; gamma_c = []; K_c = []; D_c = []; nss_c = [];

for r = 1:num_runs
    blk = fit_blocks{r};                 % n_opt_runs x 5 [omega0 gamma K D nss]
    if isempty(blk); continue; end
    nb      = size(blk, 1);
    Name_c  = [Name_c;  repmat(names_all(r), nb, 1)]; %#ok<AGROW>
    side_c  = [side_c;  repmat(sides_all(r), nb, 1)]; %#ok<AGROW>
    omega_c = [omega_c; blk(:,1)]; %#ok<AGROW>
    gamma_c = [gamma_c; blk(:,2)]; %#ok<AGROW>
    K_c     = [K_c;     blk(:,3)]; %#ok<AGROW>
    D_c     = [D_c;     blk(:,4)]; %#ok<AGROW>
    nss_c   = [nss_c;   blk(:,5)]; %#ok<AGROW>
end

T_fits = table(Name_c, side_c, omega_c, gamma_c, K_c, D_c, nss_c, ...
    'VariableNames', {'Name','side','omega0','gamma','K','D','nss'});

writetable(T_fits, fullfile(results_folder, 'all_fits.csv'));
fprintf('Results saved to: %s\n', results_folder);


%% =========================================================
%  LOCAL FUNCTIONS
%% =========================================================

% ---------------------------------------------------------
%  WAITBAR UPDATE (runs on client via afterEach)
% ---------------------------------------------------------
function update_runbar(wb, total)
    persistent count;
    if isempty(count); count = 0; end
    count = count + 1;
    if ishandle(wb)
        waitbar(count / total, wb, sprintf('Rows completed: %d / %d', count, total));
    end
end

% ---------------------------------------------------------
%  KURAMOTO — VECTORISED (O(N) mean-field)
%  Returns ALL nRep realisations (nRep x nSteps), each z-scored.
%  No rho_avg computed (not needed).
% ---------------------------------------------------------
function lfp_all = run_kuramoto_fast(omega0, gamma, K_val, D, opt_params)

    dt          = opt_params.dt;
    t_exp       = opt_params.t_exp;
    N           = opt_params.N;
    nRep        = opt_params.nRep;

    nSteps      = round(t_exp / dt) + 1;
    noise_scale = sqrt(2 * D * dt);

    lfp_all = zeros(nRep, nSteps);

    for rep = 1:nRep

        omegas = 2 * pi * cauchy_rnd(omega0, gamma, 1, N);   % 1xN
        theta  = 2 * pi * rand(1, N);                        % 1xN initial phases
        lfp_rep = zeros(1, nSteps);

        Z          = mean(exp(1i * theta));
        lfp_rep(1) = real(Z);

        for i = 1:nSteps - 1
            Z        = mean(exp(1i * theta));                % order parameter
            coupling = K_val * imag(Z * exp(-1i * theta));   % 1xN, O(N)
            theta    = theta + (omegas + coupling) * dt ...
                             + noise_scale * randn(1, N);
            Z            = mean(exp(1i * theta));
            lfp_rep(i+1) = real(Z);
        end

        s = std(lfp_rep);
        if s > 0; lfp_rep = (lfp_rep - mean(lfp_rep)) / s; end
        lfp_all(rep, :) = lfp_rep;
    end
end

% ---------------------------------------------------------
%  FEATURISE A MODEL SIGNAL
%  resample to Fs_exp -> bandpass omega0 +/- bw_bp -> z-score -> features.
%  Matches the experimental signal's treatment exactly.
% ---------------------------------------------------------
function feat = featurize_model(sig, Fs_model, Fs_exp, omega0, bw_bp)
    sig = double(sig(:)');
    if abs(Fs_model - Fs_exp) > 0.5
        [p, q] = rat(Fs_exp / Fs_model, 1e-6);
        sig    = resample(sig, p, q);
    end
    bp  = matlab_bpass_filt(sig, omega0 - bw_bp, omega0 + bw_bp, Fs_exp);
    bp  = bp(:)';
    s   = std(bp);
    if s > 0; bp = (bp - mean(bp)) / s; end
    feat = compute_features(bp, Fs_exp);
end

% ---------------------------------------------------------
%  OBJECTIVE WRAPPER
%  Simulate -> mean across reps -> featurise -> NSS error.
% ---------------------------------------------------------
function f = obj_wrapper(x, omega0, exp_features, opt_params, Fs_model, bw_bp)

    lfp_all = run_kuramoto_fast(omega0, x(1), x(2), x(3), opt_params);
    sig_avg = mean(lfp_all, 1);    % mean of z-scored realisations

    model_feat = featurize_model(sig_avg, Fs_model, opt_params.Fs_exp, ...
                                 omega0, bw_bp);
    f = error_function(exp_features, model_feat);
end

% ---------------------------------------------------------
%  BUILD MODEL BANDS (mean + percentile envelope over nRep reps)
%  All curves aligned to the experimental grids.
% ---------------------------------------------------------
function mb = build_model_bands(lfp_all, Fs_model, Fs_exp, ...
                                omega0, bw_bp, exp_features, pct_band)

    nRep    = size(lfp_all, 1);
    f_sig   = exp_features.signal_freq(:)';
    f_env   = exp_features.env_psd_freq(:)';
    bins    = exp_features.env_pdf_bins(:)';

    SIG = zeros(nRep, numel(f_sig));
    ENV = zeros(nRep, numel(f_env));
    PDF = zeros(nRep, numel(bins));

    for rep = 1:nRep
        feat = featurize_model(lfp_all(rep, :), Fs_model, Fs_exp, omega0, bw_bp);

        % PSDs share the grid (same Fs & window) -> length-align defensively
        SIG(rep, :) = interp_to_grid(feat.signal_psd, numel(f_sig));
        ENV(rep, :) = interp_to_grid(feat.env_psd,    numel(f_env));

        % env PDF bins differ per rep -> interpolate onto experimental bins
        PDF(rep, :) = interp1(feat.env_pdf_bins(:)', feat.env_pdf(:)', ...
                              bins, 'linear', 0);
    end

    lo = pct_band(1); hi = pct_band(2);

    mb.signal_freq    = f_sig;
    mb.env_psd_freq   = f_env;
    mb.env_pdf_bins   = bins;

    mb.signal_psd_mean = mean(SIG, 1);
    mb.signal_psd_lo   = prctile(SIG, lo, 1);
    mb.signal_psd_hi   = prctile(SIG, hi, 1);

    mb.env_psd_mean    = mean(ENV, 1);
    mb.env_psd_lo      = prctile(ENV, lo, 1);
    mb.env_psd_hi      = prctile(ENV, hi, 1);

    mb.env_pdf_mean    = mean(PDF, 1);
    mb.env_pdf_lo      = prctile(PDF, lo, 1);
    mb.env_pdf_hi      = prctile(PDF, hi, 1);
end

% ---------------------------------------------------------
%  BANDPASS FILTER (zero-phase FIR, pure MATLAB)
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

    window   = round(Fs * 1);
    noverlap = round(window / 2);

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
%  ERROR FUNCTION (mean NSS across the three features)
% ---------------------------------------------------------
function f = error_function(exp_feat, model_feat)
    f = (1/3) * (nss(exp_feat.signal_psd, model_feat.signal_psd) + ...
                 nss(exp_feat.env_pdf,     model_feat.env_pdf   ) + ...
                 nss(exp_feat.env_psd,     model_feat.env_psd   ));
end

% ---------------------------------------------------------
%  NASH-SUTCLIFFE (here returned as SS ratio; lower = better)
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
%  FILENAME SANITISER
% ---------------------------------------------------------
function out = sanitise(str)
    out = regexprep(char(str), '[^a-zA-Z0-9]', '');
    if isempty(out); out = 'x'; end
end

% ---------------------------------------------------------
%  COMPARISON PLOT (3 feature panels + parameter table)
%  Experimental = one colour; Model = mean line + shaded
%  5-95th percentile band. All lines use 'LineJoin','chamfer'.
% ---------------------------------------------------------
function plot_comparison(exp_feat, mb, omega0, x_best, nss_best, ...
                         bw, pct_band, base_name, results_folder)

    c_exp = [0.12 0.47 0.71];     % experimental
    c_mod = [1.00 0.50 0.05];     % model
    lw    = 2.5;
    LJ    = {'LineJoin','chamfer'};

    fig = figure('Color','w','Units','normalized', ...
                 'Position',[0.05 0.15 0.92 0.55],'Visible','off');

    sgtitle(sprintf(['Kuramoto Fit  |  %s  |  \\omega_0 = %.3f Hz (fixed)  |  ' ...
        'BP \\pm %.1f Hz'], strrep(base_name,'_',' '), omega0, bw), ...
        'FontSize', 11, 'FontWeight', 'bold');

    san = @(v) max(real(v(:)'), 1e-12);

    f_sig = exp_feat.signal_freq(:)';
    f_env = exp_feat.env_psd_freq(:)';
    bins  = exp_feat.env_pdf_bins(:)';

    psd_e  = san(exp_feat.signal_psd);
    epsd_e = san(exp_feat.env_psd);
    pdf_e  = san(exp_feat.env_pdf);

    band_lbl = sprintf('Model %d–%dth pct', pct_band(1), pct_band(2));

    % ---- Panel 1: Signal PSD ----
    subplot(1,4,1);
    shaded_band(f_sig, san(mb.signal_psd_lo), san(mb.signal_psd_hi), c_mod);
    hold on;
    plot(f_sig, san(mb.signal_psd_mean), 'Color', c_mod, 'LineWidth', lw, LJ{:});
    plot(f_sig, psd_e, 'Color', c_exp, 'LineWidth', lw, LJ{:});
    xlabel('Frequency (Hz)'); ylabel('PSD');
    title('Signal PSD', 'FontSize', 9);
    legend({band_lbl,'Model mean','Experimental'}, 'Location','northeast','FontSize',7);
    xlim([omega0-8, omega0+8]); grid on; box on;

    % ---- Panel 2: Envelope PDF ----
    subplot(1,4,2);
    shaded_band(bins, san(mb.env_pdf_lo), san(mb.env_pdf_hi), c_mod);
    hold on;
    plot(bins, san(mb.env_pdf_mean), 'Color', c_mod, 'LineWidth', lw, LJ{:});
    plot(bins, pdf_e, 'Color', c_exp, 'LineWidth', lw, LJ{:});
    xlabel('Envelope Amplitude'); ylabel('PDF');
    title('Envelope PDF', 'FontSize', 9);
    legend({band_lbl,'Model mean','Experimental'}, 'Location','northeast','FontSize',7);
    grid on; box on;

    % ---- Panel 3: Envelope PSD ----
    subplot(1,4,3);
    shaded_band(f_env, san(mb.env_psd_lo), san(mb.env_psd_hi), c_mod);
    hold on;
    plot(f_env, san(mb.env_psd_mean), 'Color', c_mod, 'LineWidth', lw, LJ{:});
    plot(f_env, epsd_e, 'Color', c_exp, 'LineWidth', lw, LJ{:});
    xlabel('Frequency (Hz)'); ylabel('PSD');
    title('Envelope PSD', 'FontSize', 9);
    legend({band_lbl,'Model mean','Experimental'}, 'Location','northeast','FontSize',7);
    xlim([0 5]); grid on; box on;

    % ---- Panel 4: Parameter table ----
    subplot(1,4,4); axis off;
    tab_str = sprintf( ...
        ['Best-fit parameters\n\n' ...
         '  omega_0 = %.4f Hz  (fixed)\n\n' ...
         '  gamma   = %.4f Hz\n\n' ...
         '  K       = %.3f\n\n' ...
         '  D       = %.4f\n\n' ...
         'Bandpass BW = %.1f Hz\n\n' ...
         'NSS error = %.5f'], ...
        omega0, x_best(1), x_best(2), x_best(3), bw*2, nss_best);
    text(0.02, 0.97, tab_str, 'Units','normalized', ...
        'VerticalAlignment','top','FontSize',9,'FontName','FixedWidth', ...
        'BackgroundColor',[0.96 0.96 0.96], ...
        'EdgeColor',[0.5 0.5 0.5],'LineWidth',1.2,'Interpreter','none');
    title('Best-Fit Parameters', 'FontSize', 9);

    % ---- Save ----
    png_path = fullfile(results_folder, [base_name '_comparison.png']);
    fig_path = fullfile(results_folder, [base_name '_comparison.fig']);
    set(fig, 'PaperPositionMode', 'auto');
    print(fig, png_path, '-dpng', '-r150');
    set(fig, 'Visible', 'on');
    savefig(fig, fig_path);
    close(fig);
end

% ---------------------------------------------------------
%  SHADED PERCENTILE BAND (fill between lo and hi)
% ---------------------------------------------------------
function shaded_band(x, lo, hi, col)
    x  = x(:)'; lo = lo(:)'; hi = hi(:)';
    fill([x, fliplr(x)], [lo, fliplr(hi)], col, ...
        'FaceAlpha', 0.20, 'EdgeColor', 'none', 'LineJoin', 'chamfer');
end