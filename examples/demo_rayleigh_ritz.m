P = @(x) 1 + 0*x;
Q = @(x) 0*x;
f = @(x) pi^2 * sin(pi*x);

[~, ~, solutions] = rayleigh_ritz_bvp(P, Q, f, 0, 1, 8, 32);
sol = solutions{end};

plot(sol.x, sol.v, '-o', 'LineWidth', 1.2);
xlabel('x');
ylabel('v(x)');
title('Rayleigh-Ritz Approximation');
grid on;
