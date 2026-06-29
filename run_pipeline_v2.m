
function results = run_pipeline_v2(csvfile, outfile, rows, save_figs)
% RUN_PIPELINE_V2  v2 driver implementing the revised protocol.
%
%   For each parameter set in the table:
%     Step 1  40 s resting run, discard 10 s burn-in, measure rho_baseline
%             and SD over the remaining 30 s, store the oscillator END-STATE.
%     Step 2  choose c (noise-floor rule) and kick width sigma_kick.
%     Step 3  24 perturbation realizations (parfor), all starting from the
%             stored end-state; compute mean + 5-95 percentile band; compute
%             settling time (first sustained entry into +/-2 SD of baseline).
%     Step 5  two-panel .fig figure (before: 40 s with burn-in marked;
%             after: mean + 5-95 band, +/-2 SD band, settling time marker).
%
% Usage:
%   results = run_pipeline_v2();                         % all rows, defaults
%   results = run_pipeline_v2(csv, out, [1 3 5], true);  % subset + save .fig
%
% Inputs (all optional):
%   csvfile   default 'less_nss_all_THE_fit_with_KANG.csv'
%   outfile   default 'THE_fit_v2_results.mat'
%   rows      row indices to process (default: all)
%   save_figs save each figure as .fig  (default false)

    if nargin < 1 || isempty(csvfile);   csvfile   = 'less_nss_all_THE_fit_with_KANG.csv'; end
    if nargin < 2 || isempty(outfile);   outfile   = 'THE_fit_v2_results.mat'; end
    if nargin < 4 || isempty(save_figs); save_figs = false; end

    T   = readtable(csvfile);
    cfg = default_config_v2();
    if nargin < 3 || isempty(rows); rows = 1:height(T); end

    % Start a parallel pool if one is not already running (used by parfor in Step 3)
    if isempty(gcp('nocreate'))
        try
            parpool;
        catch ME
            warning('Could not start parpool (%s). parfor will run serially.', ME.message);
        end
    end

    nR = numel(rows);

    % pre-allocate result columns
    rho_baseline  = nan(nR, 1);
    rho_std_col   = nan(nR, 1);
    c_used        = nan(nR, 1);
    meets_floor   = false(nR, 1);
    settle_med    = nan(nR, 1);
    settle_min    = nan(nR, 1);
    settle_max    = nan(nR, 1);
    n_settled     = nan(nR, 1);
    trim_len_s    = nan(nR, 1);
    Kc_col        = nan(nR, 1);
    is_synced_col = false(nR, 1);
    patient_row   = rows(:);

    for idx = 1:nR
        row = rows(idx);

        pp.omega0 = T.omega0(row);
        pp.gamma  = T.gamma(row);
        pp.K      = T.K(row);
        pp.D      = T.D(row);

        Kc            = 4*pi*pp.gamma + 2*pp.D;
        Kc_col(idx)   = Kc;
        is_synced_col(idx) = pp.K > Kc;

        % ---- chained steps ----
        rest = step1_resting_sim_v2(pp, cfg);
        kick = step2_choose_c_and_kick_v2(rest, cfg);
        pert = step3_run_perturbation_v2(rest, kick, pp, cfg);
        step5_plot_v2(rest, pert, cfg, pp, row, save_figs);

        % ---- collect metrics ----
        % step3 returns a VECTOR of 24 per-realization settle times
        % (pert.settle_times); summarize it to one number per patient.
        rho_baseline(idx)  = rest.rho_baseline;
        rho_std_col(idx)   = rest.rho_std;
        c_used(idx)        = kick.c;
        meets_floor(idx)   = kick.meets_floor;

        st = pert.settle_times;                  % 1x24, NaN where never settled
        settle_med(idx)    = median(st, 'omitnan');
        settle_min(idx)    = min(st, [], 'omitnan');
        settle_max(idx)    = max(st, [], 'omitnan');
        n_settled(idx)     = sum(~isnan(st));    % how many of the 24 settled
        trim_len_s(idx)    = pert.t_trim(end);   % common trimmed window length (s)

        fprintf(['Patient %3d/%d | rho* = %.4f +/- %.4f | c = %.2f | ' ...
                 'meets_floor = %d | t_settle(med) = %.3f s | %d/%d settled | synced = %d\n'], ...
                idx, nR, rho_baseline(idx), rho_std_col(idx), c_used(idx), ...
                meets_floor(idx), settle_med(idx), n_settled(idx), ...
                cfg.n_realizations, is_synced_col(idx));
    end

    results = table(patient_row, rho_baseline, rho_std_col, c_used, meets_floor, ...
                    settle_med, settle_min, settle_max, n_settled, trim_len_s, ...
                    Kc_col, is_synced_col, ...
        'VariableNames', {'row', 'rho_baseline', 'rho_std', 'c_used', ...
                          'meets_floor', 'settle_med_s', 'settle_min_s', ...
                          'settle_max_s', 'n_settled', 'trim_len_s', ...
                          'Kc', 'is_synced'});

    save(outfile, 'results', 'cfg');
    fprintf('\nSaved results to  %s\n', outfile);
end