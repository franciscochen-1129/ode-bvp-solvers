P = @(x) 0;
Q = @(x) -1;
R = @(x) -x;

[~, ~, solutions] = finite_difference_bvp(P, Q, R, 0, 1, 0, 0, 8, 32);
sol = solutions{end};

plot(sol.x, sol.u, '-o', 'LineWidth', 1.2);
xlabel('x');
ylabel('u(x)');
title('Finite Difference Solution');
grid on;
