function [x, y, s_final, it] = shooting_rk2(F, a, b, alpha, beta, n, s0, s1, tol, maxit)
%SHOOTING_RK2 Legacy wrapper for shooting_bvp_rk2.

[x, y, s_final, it] = shooting_bvp_rk2(F, a, b, alpha, beta, n, s0, s1, tol, maxit);
end
