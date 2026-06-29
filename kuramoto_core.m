function [rho_series, theta_snaps, theta_final] = kuramoto_core(omegas, theta0, K, D, dt, nSteps, snap_idx)
% KURAMOTO_CORE  Euler-Maruyama integration of the noisy Kuramoto model.
%   Mirrors the user's run_kuramoto_fast stepping scheme, but:
%     - tracks rho = |Z| (the order-parameter MAGNITUDE), not real(Z)
%     - does NOT z-score (absolute rho is needed for the OU / burst analysis)
%     - can be initialised from an arbitrary theta0 (for perturbation recovery)
%
% Inputs:
%   omegas   1xN natural frequencies (rad/s), already multiplied by 2*pi
%   theta0   1xN initial phases
%   K, D     coupling and noise (same units/conventions as the fit)
%   dt       timestep (s)
%   nSteps   number of recorded samples (including the initial one)
%   snap_idx vector of step indices at which to store theta (may be empty)
%
% Outputs:
%   rho_series  1 x nSteps,  rho(t) = |mean(exp(1i*theta))|
%   theta_snaps numel(snap_idx) x N,  theta at requested indices
%   theta_final 1 x N,  final phase configuration

    N = numel(theta0);
    noise_scale = sqrt(2 * D * dt);
    theta = theta0(:).';                 % ensure row

    rho_series = zeros(1, nSteps);
    if nargin < 7 || isempty(snap_idx); snap_idx = []; end
    nSnap = numel(snap_idx);
    theta_snaps = zeros(nSnap, N);
    snapPtr = 1;

    Z = mean(exp(1i * theta));
    rho_series(1) = abs(Z);
    if nSnap > 0 && snap_idx(1) == 1
        theta_snaps(1, :) = theta; snapPtr = 2;
    end

    for i = 1:nSteps - 1
        Z        = mean(exp(1i * theta));
        coupling = K * imag(Z * exp(-1i * theta));
        theta    = theta + (omegas + coupling) * dt + noise_scale * randn(1, N);
        Z        = mean(exp(1i * theta));
        rho_series(i + 1) = abs(Z);
        if snapPtr <= nSnap && (i + 1) == snap_idx(snapPtr)
            theta_snaps(snapPtr, :) = theta;
            snapPtr = snapPtr + 1;
        end
    end
    theta_final = theta;
end
