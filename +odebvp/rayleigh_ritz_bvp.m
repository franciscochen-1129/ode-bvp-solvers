function [N_values, h_values, solutions] = rayleigh_ritz_bvp(P, Q, f, a, b, Nmin, Nmax)
%RAYLEIGH_RITZ_BVP Solve a homogeneous-endpoint BVP with piecewise linear FEM.
%   Solves
%       -(P(x) v'(x))' + Q(x) v(x) = f(x),  a < x < b
%       v(a) = v(b) = 0
%   with piecewise linear basis functions on a uniform mesh.

    narginchk(7, 7);
    validate_inputs_local(P, Q, f, a, b, Nmin, Nmax);

    N_values = Nmin:Nmax;
    h_values = (b - a) ./ N_values;
    solutions = cell(numel(N_values), 1);

    for idx = 1:numel(N_values)
        N = N_values(idx);
        h = (b - a) / N;
        x_nodes = linspace(a, b, N + 1);

        A = zeros(N - 1, N - 1);
        b_vec = zeros(N - 1, 1);

        for i = 1:N-1
            x_left = x_nodes(i);
            x_right = x_nodes(i + 2);
            b_vec(i) = integral(@(x) f(x) .* phi_local(x, i, h, a), x_left, x_right);
        end

        for i = 1:N-1
            x_left = x_nodes(i);
            x_mid = x_nodes(i + 1);
            x_right = x_nodes(i + 2);

            A(i, i) = integral(@(x) ...
                P(x) .* phi_prime_local(x, i, h, a) .* phi_prime_local(x, i, h, a) + ...
                Q(x) .* phi_local(x, i, h, a) .* phi_local(x, i, h, a), ...
                x_left, x_right);

            if i > 1
                val = integral(@(x) ...
                    P(x) .* phi_prime_local(x, i, h, a) .* phi_prime_local(x, i - 1, h, a) + ...
                    Q(x) .* phi_local(x, i, h, a) .* phi_local(x, i - 1, h, a), ...
                    x_left, x_mid);
                A(i, i - 1) = val;
                A(i - 1, i) = val;
            end

            if i < N - 1
                val = integral(@(x) ...
                    P(x) .* phi_prime_local(x, i, h, a) .* phi_prime_local(x, i + 1, h, a) + ...
                    Q(x) .* phi_local(x, i, h, a) .* phi_local(x, i + 1, h, a), ...
                    x_mid, x_right);
                A(i, i + 1) = val;
                A(i + 1, i) = val;
            end
        end

        c = A \ b_vec;
        v_nodes = zeros(N + 1, 1);
        v_nodes(2:N) = c;

        sol = struct();
        sol.x = x_nodes(:);
        sol.v = v_nodes;
        sol.N = N;
        sol.h = h;
        solutions{idx} = sol;
    end
end

function val = phi_local(x, i, h, a)
xi_minus = a + (i - 1) * h;
xi = a + i * h;
xi_plus = a + (i + 1) * h;

val = zeros(size(x));
left = (x >= xi_minus) & (x <= xi);
right = (x >= xi) & (x <= xi_plus);

val(left) = (x(left) - xi_minus) / h;
val(right) = (xi_plus - x(right)) / h;
end

function val = phi_prime_local(x, i, h, a)
xi_minus = a + (i - 1) * h;
xi = a + i * h;
xi_plus = a + (i + 1) * h;

val = zeros(size(x));
left = (x >= xi_minus) & (x <= xi);
right = (x >= xi) & (x <= xi_plus);

val(left) = 1 / h;
val(right) = -1 / h;
end

function validate_inputs_local(P, Q, f, a, b, Nmin, Nmax)
if ~isa(P, 'function_handle') || ~isa(Q, 'function_handle') || ~isa(f, 'function_handle')
    error('P, Q, and f must be function handles.');
end
validateattributes(a, {'numeric'}, {'scalar', 'real', 'finite'});
validateattributes(b, {'numeric'}, {'scalar', 'real', 'finite', '>', a});
validateattributes(Nmin, {'numeric'}, {'scalar', 'integer', '>=', 2});
validateattributes(Nmax, {'numeric'}, {'scalar', 'integer', '>=', Nmin});
end
