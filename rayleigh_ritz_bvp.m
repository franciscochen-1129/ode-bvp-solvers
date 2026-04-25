function [N_values, h_values, solutions] = rayleigh_ritz_bvp(P, Q, f, a, b, Nmin, Nmax)
%RAYLEIGH_RITZ_BVP Public entry point for the Rayleigh-Ritz solver.

[N_values, h_values, solutions] = odebvp.rayleigh_ritz_bvp(P, Q, f, a, b, Nmin, Nmax);
end
