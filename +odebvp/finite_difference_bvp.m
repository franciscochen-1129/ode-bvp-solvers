function [N_values, h_values, solutions] = finite_difference_bvp(P, Q, R, a, b, alpha, beta, Nmin, Nmax)
%FINITE_DIFFERENCE_BVP Solve a linear two-point BVP with central differences.
%   Solves
%       u''(x) + P(x) u'(x) + Q(x) u(x) = R(x),  a < x < b
%       u(a) = alpha,  u(b) = beta
%   on a sequence of uniform grids using a second-order finite-difference
%   discretization.
%
%   The returned cell array stores one struct per grid size with fields:
%       .x  grid nodes
%       .u  nodal solution values

    narginchk(9, 9);
    validate_inputs_local(P, Q, R, a, b, Nmin, Nmax);

    N_values = Nmin:Nmax;
    h_values = (b - a) ./ N_values;
    solutions = cell(numel(N_values), 1);

    for idx = 1:numel(N_values)
        N = N_values(idx);
        h = (b - a) / N;
        x_nodes = linspace(a, b, N + 1).';

        A = zeros(N - 1, N - 1);
        b_vec = zeros(N - 1, 1);

        for i = 1:N-1
            x_i = a + i * h;
            Pi = P(x_i);
            Qi = Q(x_i);
            Ri = R(x_i);

            if i == 1
                A(i, i) = 2 + h^2 * Qi;
                A(i, i + 1) = -1 + (h / 2) * Pi;
                b_vec(i) = -h^2 * Ri + alpha * (1 + (h / 2) * Pi);
            elseif i == N - 1
                A(i, i - 1) = -1 - (h / 2) * Pi;
                A(i, i) = 2 + h^2 * Qi;
                b_vec(i) = -h^2 * Ri + beta * (1 - (h / 2) * Pi);
            else
                A(i, i - 1) = -1 - (h / 2) * Pi;
                A(i, i) = 2 + h^2 * Qi;
                A(i, i + 1) = -1 + (h / 2) * Pi;
                b_vec(i) = -h^2 * Ri;
            end
        end

        u_inner = A \ b_vec;

        sol = struct();
        sol.x = x_nodes;
        sol.u = [alpha; u_inner; beta];
        sol.N = N;
        sol.h = h;
        solutions{idx} = sol;
    end
end

function validate_inputs_local(P, Q, R, a, b, Nmin, Nmax)
if ~isa(P, 'function_handle') || ~isa(Q, 'function_handle') || ~isa(R, 'function_handle')
    error('P, Q, and R must be function handles.');
end
validateattributes(a, {'numeric'}, {'scalar', 'real', 'finite'});
validateattributes(b, {'numeric'}, {'scalar', 'real', 'finite', '>', a});
validateattributes(Nmin, {'numeric'}, {'scalar', 'integer', '>=', 2});
validateattributes(Nmax, {'numeric'}, {'scalar', 'integer', '>=', Nmin});
end
