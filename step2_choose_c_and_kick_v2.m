function kick = step2_choose_c_and_kick_v2(rest, cfg)
% STEP2_CHOOSE_C_AND_KICK_V2  Pick retained fraction c by the noise-floor
%   rule, then convert to the per-oscillator kick width.
%
%   (a) require (1-c)*rho_baseline > factor*sigma_floor, reduce c until met.
%   (b) sigma_kick = sqrt(-2*ln c)   (rho_baseline cancels -> standardized deficit)
%
% Input : rest (step1_v2), cfg
% Output: kick struct with c, sigma_kick, deficit, target, meets_floor

    c      = cfg.c_init;
    factor = cfg.noise_floor_factor;

    while ((1 - c) * rest.rho_baseline <= factor * rest.sigma_floor) && (c > cfg.c_min)
        c = c - cfg.c_step;
    end
    c = max(c, cfg.c_min);

    kick.c           = c;
    kick.sigma_kick  = sqrt(-2 * log(c));
    kick.deficit     = (1 - c) * rest.rho_baseline;
    kick.target      = c * rest.rho_baseline;
    kick.meets_floor = kick.deficit > factor * rest.sigma_floor;
end
