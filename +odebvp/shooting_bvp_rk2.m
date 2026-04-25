function [x, y, s_final, it] = shooting_bvp_rk2(F, a, b, alpha, beta, n, s0, s1, tol, maxit)
%SHOOTING_BVP_RK2 Solve a scalar second-order BVP by shooting + RK2 + secant.
%   Solves
%       y'' = F(x, y, y')
%       y(a) = alpha,  y(b) = beta
%   by converting to an IVP in the unknown initial slope y'(a).

    narginchk(8, 10);
    if nargin < 10
        maxit = 20;
    end
    if nargin < 9
        tol = 1e-6;
    end

    validate_inputs_local(F, a, b, n, tol, maxit);

    N = n - 1;
    h = (b - a) / N;
    x = linspace(a, b, n);

    [y_end0, ~] = solve_with_s_local(s0);
    [y_end1, y_all1] = solve_with_s_local(s1);

    it = 1;
    while abs(y_end1 - beta) > tol && it < maxit
        denom = y_end1 - y_end0;
        if abs(denom) < eps
            error('Secant update failed because consecutive shooting residuals are identical.');
        end

        s_new = s1 - (y_end1 - beta) * (s1 - s0) / denom;

        s0 = s1;
        y_end0 = y_end1;
        s1 = s_new;
        [y_end1, y_all1] = solve_with_s_local(s1);
        it = it + 1;
    end

    y = y_all1;
    s_final = s1;

    function [y_end, y_all] = solve_with_s_local(s)
        Y = zeros(2, n);
        Y(:, 1) = [alpha; s];

        for k = 1:N
            xk = a + (k - 1) * h;
            yk = Y(1, k);
            ypk = Y(2, k);

            K1 = [ypk; F(xk, yk, ypk)];
            y_temp = yk + h * K1(1);
            yp_temp = ypk + h * K1(2);
            K2 = [yp_temp; F(xk + h, y_temp, yp_temp)];

            Y(:, k + 1) = Y(:, k) + (h / 2) * (K1 + K2);
        end

        y_all = Y(1, :);
        y_end = y_all(end);
    end
end

function validate_inputs_local(F, a, b, n, tol, maxit)
if ~isa(F, 'function_handle')
    error('F must be a function handle.');
end
validateattributes(a, {'numeric'}, {'scalar', 'real', 'finite'});
validateattributes(b, {'numeric'}, {'scalar', 'real', 'finite', '>', a});
validateattributes(n, {'numeric'}, {'scalar', 'integer', '>=', 2});
validateattributes(tol, {'numeric'}, {'scalar', 'real', 'positive', 'finite'});
validateattributes(maxit, {'numeric'}, {'scalar', 'integer', '>=', 1});
end
