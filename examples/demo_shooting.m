F = @(x, y, yp) -y;
[x, y, s_final, it] = shooting_bvp_rk2(F, 0, pi/2, 0, 1, 100, 0.5, 1.5, 1e-8, 25); %#ok<ASGLU>

plot(x, y, 'LineWidth', 1.5);
xlabel('x');
ylabel('y(x)');
title(sprintf('Shooting RK2 Solution (iterations = %d)', it));
grid on;
