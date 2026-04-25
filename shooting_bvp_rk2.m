function [x, y, s_final, it] = shooting_bvp_rk2(F, a, b, alpha, beta, n, s0, s1, tol, maxit)
%SHOOTING_BVP_RK2 Public entry point for the shooting solver.

[x, y, s_final, it] = odebvp.shooting_bvp_rk2(F, a, b, alpha, beta, n, s0, s1, tol, maxit);
end
