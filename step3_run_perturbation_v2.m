function pert = step3_run_perturbation_v2(rest, kick, pp, cfg)
% STEP3_RUN_PERTURBATION_V2  24 perturbation realizations per patient.
%
% Changes from previous version:
%   - Kick tolerance is now relative to the finite-N noise floor 1/sqrt(2N)
%     rather than to rho_baseline. This ensures the tolerance is achievable.
%   - Settling criterion uses a smoothed mean trace (cfg.smooth_s window)
%     before applying the band test, removing single-sample excursions.
%   - settle_hold defaults to 2 s (was 30 s).
%
% Per realization:
%   1. Start from the stored oscillator end-state (rest.theta_end).
%   2. Apply dispersing phase kick of width sigma_kick.
%   3. Simulate forward in chunks.
%   4. Continue until SMOOTHED mean stays within rho_baseline +/- 2 SD
%      for settle_hold seconds continuously.
%   5. Store rho(t) for that realization.

    N  = cfg.N;
    dt = cfg.dt;
    R  = cfg.n_realizations;

    % broadcast scalars for parfor
    theta_end  = rest.theta_end;
    omegas     = rest.omegas;
    K          = pp.K;
    D          = pp.D;
    sigma_kick = kick.sigma_kick;
    target     = kick.target;
    maxAtt     = cfg.max_kick_attempts;
    rho_base   = rest.rho_baseline;
    rho_std    = rest.rho_std;
    band       = cfg.settle_nstd * rho_std;
    hold_n     = round(cfg.settle_hold / dt);
    chunk_n    = round(cfg.chunk_s / dt);
    t_max_n    = round(cfg.t_max_recover / dt);
    smooth_n   = max(1, round(cfg.smooth_s / dt));

    % FIX 1: tolerance relative to finite-N noise floor
    % For N oscillators the std of |Z| at any instant is ~1/sqrt(2N)
    tol = cfg.kick_tol_factor / sqrt(2 * N);

    traces_cell  = cell(1, R);
    settle_steps = nan(1, R);

    parfor r = 1:R

        % ---- kick from stored end-state ----
        theta_kicked = theta_end;
        best_xi      = sigma_kick * randn(1, N);  % fallback if never within tol
        best_err     = Inf;

        for attempt = 1:maxAtt
            xi           = sigma_kick * randn(1, N);
            theta_try    = theta_end + xi;
            rho_realized = abs(mean(exp(1i * theta_try)));
            err          = abs(rho_realized - target);

            if err < best_err
                best_err  = err;
                best_xi   = xi;
            end
            if err <= tol
                break;
            end
        end
        theta_kicked = theta_end + best_xi;   % best draw found in maxAtt attempts

        % ---- simulate in chunks until settled ----
        theta_now  = theta_kicked;
        rho_trace  = [];
        settled_at = NaN;

        while numel(rho_trace) < t_max_n

            n_this = min(chunk_n, t_max_n - numel(rho_trace));
            rho_chunk = kuramoto_core(omegas, theta_now, K, D, dt, n_this + 1, []);

            if isempty(rho_trace)
                rho_trace = rho_chunk;
            else
                rho_trace = [rho_trace, rho_chunk(2:end)]; %#ok<AGROW>
            end

            % advance theta_now
            [~, theta_now] = kuramoto_core(omegas, theta_now, K, D, dt, n_this + 1, n_this + 1);

            % FIX 2: smooth rho_trace before checking the band
            % (removes single-sample excursions that block settling detection)
            if numel(rho_trace) >= smooth_n
                rho_smooth = conv(rho_trace, ones(1, smooth_n)/smooth_n, 'same');
            else
                rho_smooth = rho_trace;
            end

            in_band = abs(rho_smooth - rho_base) <= band;
            n_tr    = numel(in_band);

            for i = max(1, n_tr - chunk_n - hold_n) : max(1, n_tr - hold_n)
                if i + hold_n <= n_tr && all(in_band(i : i + hold_n))
                    settled_at = i + hold_n;
                    break;
                end
            end

            if ~isnan(settled_at)
                rho_trace = rho_trace(1 : min(settled_at, numel(rho_trace)));
                break;
            end
        end

        traces_cell{r}  = rho_trace;
        settle_steps(r) = settled_at;
    end

    % ---- trim all traces to the shortest settled length ----
    lengths  = cellfun(@numel, traces_cell);
    trim_len = min(lengths);

    traces = zeros(R, trim_len);
    for r = 1:R
        traces(r, :) = traces_cell{r}(1 : trim_len);
    end

    t_trim     = (0 : trim_len - 1) * dt;
    mean_trace = mean(traces, 1);
    p05        = prctile_rows(traces, 5);
    p95        = prctile_rows(traces, 95);

    pert.traces_cell  = traces_cell;
    pert.traces       = traces;
    pert.mean_trace   = mean_trace;
    pert.p05          = p05;
    pert.p95          = p95;
    pert.t_trim       = t_trim;
    pert.settle_steps = settle_steps;
    pert.settle_times = (settle_steps - 1) * dt;
    pert.target       = target;
    pert.trim_len     = trim_len;
    pert.kick_tol     = tol;              % log for diagnostics
    pert.kick_target  = target;
end

% -----------------------------------------------------------------------
function v = prctile_rows(M, p)
    [R, C] = size(M);
    v = zeros(1, C);
    pos_target = p / 100;
    for c = 1:C
        x    = sort(M(:, c));
        pos  = ((1:R) - 0.5) / R;
        v(c) = interp1(pos, x, pos_target, 'linear', 'extrap');
    end
end

% function pert = step3_run_perturbation_v2(rest, kick, pp, cfg)
% % STEP3_RUN_PERTURBATION_V2  24 perturbation realizations per patient.
% %
% % Changes from previous version:
% %   - Kick tolerance is now relative to the finite-N noise floor 1/sqrt(2N)
% %     rather than to rho_baseline. This ensures the tolerance is achievable.
% %   - Settling criterion uses a smoothed mean trace (cfg.smooth_s window)
% %     before applying the band test, removing single-sample excursions.
% %   - Settling band is a MEDIAN-based outlier rule centered on rest.rho_median
% %     with half-width niqr * rest.rho_iqr (Tukey-style fence), replacing the
% %     previous mean +/- settle_nstd * SD band.
% %   - settle_hold defaults to 2 s (was 30 s).
% %
% % Per realization:
% %   1. Start from the stored oscillator end-state (rest.theta_end).
% %   2. Apply dispersing phase kick of width sigma_kick.
% %   3. Simulate forward in chunks.
% %   4. Continue until SMOOTHED mean stays within rho_median +/- niqr*IQR
% %      (median-based outlier band) for settle_hold seconds continuously.
% %   5. Store rho(t) for that realization.
% 
%     N  = cfg.N;
%     dt = cfg.dt;
%     R  = cfg.n_realizations;
% 
%     % broadcast scalars for parfor
%     theta_end  = rest.theta_end;
%     omegas     = rest.omegas;
%     K          = pp.K;
%     D          = pp.D;
%     sigma_kick = kick.sigma_kick;
%     target     = kick.target;
%     maxAtt     = cfg.max_kick_attempts;
% 
%     % --- median-based settling band ---
%     % Center on the baseline median, half-width = niqr * baseline IQR.
%     % To keep the band centered on the mean baseline instead, swap rho_med
%     % for rest.rho_baseline in the in_band test below.
%     rho_med    = rest.rho_median;
%     rho_iqr    = rest.rho_iqr;
%     niqr       = 0.75;                       % Tukey outlier factor
%     if isfield(cfg, 'settle_niqr')
%         niqr = cfg.settle_niqr;
%     end
%     band       = niqr * rho_iqr;
% 
%     hold_n     = round(cfg.settle_hold / dt);
%     chunk_n    = round(cfg.chunk_s / dt);
%     t_max_n    = round(cfg.t_max_recover / dt);
%     smooth_n   = max(1, round(cfg.smooth_s / dt));
% 
%     % FIX 1: tolerance relative to finite-N noise floor
%     % For N oscillators the std of |Z| at any instant is ~1/sqrt(2N)
%     tol = cfg.kick_tol_factor / sqrt(2 * N);
% 
%     traces_cell  = cell(1, R);
%     settle_steps = nan(1, R);
% 
%     parfor r = 1:R
% 
%         % ---- kick from stored end-state ----
%         theta_kicked = theta_end;
%         best_xi      = sigma_kick * randn(1, N);  % fallback if never within tol
%         best_err     = Inf;
% 
%         for attempt = 1:maxAtt
%             xi           = sigma_kick * randn(1, N);
%             theta_try    = theta_end + xi;
%             rho_realized = abs(mean(exp(1i * theta_try)));
%             err          = abs(rho_realized - target);
% 
%             if err < best_err
%                 best_err  = err;
%                 best_xi   = xi;
%             end
%             if err <= tol
%                 break;
%             end
%         end
%         theta_kicked = theta_end + best_xi;   % best draw found in maxAtt attempts
% 
%         % ---- simulate in chunks until settled ----
%         theta_now  = theta_kicked;
%         rho_trace  = [];
%         settled_at = NaN;
% 
%         while numel(rho_trace) < t_max_n
% 
%             n_this = min(chunk_n, t_max_n - numel(rho_trace));
%             rho_chunk = kuramoto_core(omegas, theta_now, K, D, dt, n_this + 1, []);
% 
%             if isempty(rho_trace)
%                 rho_trace = rho_chunk;
%             else
%                 rho_trace = [rho_trace, rho_chunk(2:end)]; %#ok<AGROW>
%             end
% 
%             % advance theta_now
%             [~, theta_now] = kuramoto_core(omegas, theta_now, K, D, dt, n_this + 1, n_this + 1);
% 
%             % FIX 2: smooth rho_trace before checking the band
%             % (removes single-sample excursions that block settling detection)
%             if numel(rho_trace) >= smooth_n
%                 rho_smooth = conv(rho_trace, ones(1, smooth_n)/smooth_n, 'same');
%             else
%                 rho_smooth = rho_trace;
%             end
% 
%             % median-based outlier band, centered on rho_med
%             in_band = abs(rho_smooth - rho_med) <= band;
%             n_tr    = numel(in_band);
% 
%             for i = max(1, n_tr - chunk_n - hold_n) : max(1, n_tr - hold_n)
%                 if i + hold_n <= n_tr && all(in_band(i : i + hold_n))
%                     settled_at = i + hold_n;
%                     break;
%                 end
%             end
% 
%             if ~isnan(settled_at)
%                 rho_trace = rho_trace(1 : min(settled_at, numel(rho_trace)));
%                 break;
%             end
%         end
% 
%         traces_cell{r}  = rho_trace;
%         settle_steps(r) = settled_at;
%     end
% 
%     % ---- trim all traces to the shortest settled length ----
%     lengths  = cellfun(@numel, traces_cell);
%     trim_len = min(lengths);
% 
%     traces = zeros(R, trim_len);
%     for r = 1:R
%         traces(r, :) = traces_cell{r}(1 : trim_len);
%     end
% 
%     t_trim     = (0 : trim_len - 1) * dt;
%     mean_trace = mean(traces, 1);
%     p05        = prctile_rows(traces, 5);
%     p95        = prctile_rows(traces, 95);
% 
%     pert.traces_cell  = traces_cell;
%     pert.traces       = traces;
%     pert.mean_trace   = mean_trace;
%     pert.p05          = p05;
%     pert.p95          = p95;
%     pert.t_trim       = t_trim;
%     pert.settle_steps = settle_steps;
%     pert.settle_times = (settle_steps - 1) * dt;
%     pert.target       = target;
%     pert.trim_len     = trim_len;
%     pert.kick_tol     = tol;              % log for diagnostics
%     pert.kick_target  = target;
%     pert.settle_band  = band;             % log the median-based half-width
%     pert.settle_center = rho_med;         % log the band center
% end
% 
% % -----------------------------------------------------------------------
% function v = prctile_rows(M, p)
%     [R, C] = size(M);
%     v = zeros(1, C);
%     pos_target = p / 100;
%     for c = 1:C
%         x    = sort(M(:, c));
%         pos  = ((1:R) - 0.5) / R;
%         v(c) = interp1(pos, x, pos_target, 'linear', 'extrap');
%     end
% end