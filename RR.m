function [N_values, h_values, solutions] = RR(P, Q, f, a, b, Nmin, Nmax)
%RR Legacy wrapper for rayleigh_ritz_bvp.

[N_values, h_values, solutions] = rayleigh_ritz_bvp(P, Q, f, a, b, Nmin, Nmax);
end
