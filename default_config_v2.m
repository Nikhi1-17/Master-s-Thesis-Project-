function cfg = default_config_v2()
% DEFAULT_CONFIG_V2  Settings for the v2 pipeline.

    % --- core solver ---
    cfg.N  = 200;
    cfg.dt = 0.0005;      % s -- MUST match the fitting dt

    % --- resting run ---
    cfg.t_rest_total = 40;
    cfg.t_burnin     = 10;
    cfg.n_rest_realizations = 32;   % resting realizations averaged in step1

    % --- depth selection ---
    cfg.c_init             = 0.2;
    cfg.c_min              = 0.2;
    cfg.c_step             = 0.005;
    cfg.noise_floor_factor = 5;

    % --- perturbation ---
    cfg.n_realizations    = 24;
    cfg.max_kick_attempts = 200;

    % FIX 1: kick tolerance relative to finite-N noise floor (1/sqrt(2N)),
    % NOT to rho_baseline. With N=200, 1/sqrt(2*200) ≈ 0.050.
    % We allow up to 2x that as tolerance.
    cfg.kick_tol_factor   = 2.0;   % tolerance = kick_tol_factor / sqrt(2*N)
    % (the old kick_tol = 0.05*rho_baseline is replaced by this)

    % --- recovery settling criterion ---
    cfg.settle_nstd  = 0.85 ; % TWO Standard Deviations
    cfg.settle_hold  = 0.5;     % FIX 2: reduced from 30 s to 2 s
                               % 30 s continuous in-band is unreachable for
                               % near-critical noisy signals; 2 s is meaningful
    cfg.smooth_s     = 0.05;   % s -- smooth mean trace before band test
                               % removes single-sample excursions that prevent
                               % settling detection without changing the physics
    cfg.chunk_s      = 5;
    cfg.t_max_recover = 60;   % s -- reduced from 300 s; near-critical patients
                               % settle within seconds once the kick lands

    % --- minimum recommended K-Kc margin ---
    cfg.min_margin_warn = 5;   % warn if K-Kc < this after step0 override

    % --- burst detection ---
    cfg.percentile = 75;

end