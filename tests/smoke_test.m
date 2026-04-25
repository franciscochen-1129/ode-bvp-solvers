function smoke_test()
%SMOKE_TEST Minimal checks that the public solvers run and return valid output.

addpath(genpath(fileparts(fileparts(mfilename('fullpath')))));

P = @(x) 0;
Q = @(x) -1;
R = @(x) -x;
[N_fd, h_fd, sols_fd] = finite_difference_bvp(P, Q, R, 0, 1, 0, 0, 8, 10); %#ok<ASGLU>
assert(numel(sols_fd) == numel(N_fd));
assert(all(isfinite(sols_fd{end}.u)));

P2 = @(x) 1 + 0*x;
Q2 = @(x) 0*x;
f2 = @(x) pi^2 * sin(pi*x);
[N_rr, h_rr, sols_rr] = rayleigh_ritz_bvp(P2, Q2, f2, 0, 1, 8, 10); %#ok<ASGLU>
assert(numel(sols_rr) == numel(N_rr));
assert(all(isfinite(sols_rr{end}.v)));

F = @(x, y, yp) -y; %#ok<INUSD>
[x, y, s_final, it] = shooting_bvp_rk2(F, 0, pi/2, 0, 1, 100, 0.5, 1.5, 1e-8, 25); %#ok<ASGLU>
assert(all(isfinite(y)));
assert(isfinite(s_final));
assert(it >= 1);

disp('smoke_test passed');
end
