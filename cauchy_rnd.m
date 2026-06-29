function x = cauchy_rnd(x0, gamma, m, n)
% CAUCHY_RND  Draw m-by-n Cauchy variates with location x0 and scale gamma.
% gamma is the SCALE parameter (= HWHM). Matches the user's existing code.
    x = x0 + gamma * tan(pi * (rand(m, n) - 0.5));
end
