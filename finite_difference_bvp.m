function [N_values, h_values, solutions] = finite_difference_bvp(P, Q, R, a, b, alpha, beta, Nmin, Nmax)
%FINITE_DIFFERENCE_BVP Public entry point for the finite-difference solver.

[N_values, h_values, solutions] = odebvp.finite_difference_bvp( ...
    P, Q, R, a, b, alpha, beta, Nmin, Nmax);
end
