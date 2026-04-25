function [N_values, h_values, solutions] = Finite_difference(P, Q, R, a, b, alpha, beta, Nmin, Nmax)
%FINITE_DIFFERENCE Legacy wrapper for finite_difference_bvp.

[N_values, h_values, solutions] = finite_difference_bvp( ...
    P, Q, R, a, b, alpha, beta, Nmin, Nmax);
end
