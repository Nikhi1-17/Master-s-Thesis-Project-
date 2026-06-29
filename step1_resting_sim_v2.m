function rest = step1_resting_sim_v2(pp, cfg)
% STEP1_RESTING_SIM_V2  One 40 s resting run; discard first 10 s; measure
%   rho_baseline and SD over the remaining 30 s; STORE the final oscillator
%   state (theta_end) so every perturbation starts from this same instance.
%
% Input
%   pp  : struct with omega0, gamma, K, D
%   cfg : default_config_v2()
% Output (struct 'rest')
%   omegas        quenched natural frequencies (rad/s)
%   theta_end     1xN final phase state (seeds all perturbations)
%   rho_full      full 40 s rho(t) (for the "before" plot)
%   rho_stat      stationary 30 s portion
%   t_full        time vector for rho_full (s)
%   rho_baseline  mean of rho over stationary 30 s
%   rho_std       standard deviation of rho over stationary 30 s
%   varRho        variance (= rho_std^2)
%   sigma_floor   = rho_std
%   burnin_s      burn-in length (s)
%   dt

    N  = cfg.N;
    dt = cfg.dt;
    nSteps    = round(cfg.t_rest_total / dt) + 1;
    burnSteps = round(cfg.t_burnin / dt);

    % quenched frequency draw + random initial phases
    omegas = 2*pi * cauchy_rnd(pp.omega0, pp.gamma, 1, N);
    theta0 = 2*pi * rand(1, N);

    % integrate the full 40 s, capturing the final state as a "snapshot"
    [rho_full, theta_snaps, theta_end] = ...
        kuramoto_core(omegas, theta0, pp.K, pp.D, dt, nSteps, nSteps); %#ok<ASGLU>

    rho_full = movmean(rho_full, cfg.smooth_s); %%% Smoothening step added

    rho_stat = rho_full(burnSteps + 1 : end);          % last 30 s

    rest.omegas       = omegas;
    rest.theta_end    = theta_end;                     % <-- stored end-state
    rest.rho_full     = rho_full;
    rest.rho_stat     = rho_stat;
    rest.t_full       = (0:nSteps-1) * dt;
    rest.rho_baseline = mean(rho_stat);
    rest.rho_median   = median(rho_stat); %%% NEW ADDITION 
    rest.rho_iqr      = iqr (rho_stat);
    rest.rho_std      = std(rho_stat);
    rest.varRho       = var(rho_stat);
    rest.sigma_floor  = rest.rho_std;
    rest.burnin_s     = cfg.t_burnin;
    rest.dt           = dt;
end

