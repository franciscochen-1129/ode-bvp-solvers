function [N_values, h_values, solutions] = RR(P, Q, f, a, b, Nmin, Nmax)
% Rayleigh–Ritz method with piecewise linear basis functions on [a,b]
% for the problem:
%   -(P(x) v'(x))' + Q(x) v(x) = f(x),   a < x < b,
%   v(a) = v(b) = 0.
%
% Inputs:
%   P, Q, f   : function handles, e.g. @(x) ...
%   a, b      : interval endpoints
%   Nmin      : minimum number of subintervals N
%   Nmax      : maximum number of subintervals N
%
% Outputs:
%   N_values  : vector of N used
%   h_values  : corresponding mesh sizes h = (b-a)/N
%   solutions : cell array; solutions{k} has fields:
%                 .x  = node coordinates (N_k+1 x 1)
%                 .v  = nodal values of v (N_k+1 x 1), with v(a)=v(b)=0

    if nargin < 7
        error('Usage: RR(P, Q, f, a, b, Nmin, Nmax)');
    end

    N_values = Nmin:Nmax;
    h_values = (b - a) ./ N_values;
    solutions = cell(length(N_values), 1);

    for idx = 1:length(N_values)
        N = N_values(idx);
        h = (b - a) / N;

        % Nodes: x_0,...,x_N
        x_nodes = linspace(a, b, N+1);

        % Initialize system matrix and right-hand side
        A     = zeros(N-1, N-1);
        b_vec = zeros(N-1, 1);

        % Assemble right-hand side: b_i = ∫ f(x) φ_i(x) dx
        for i = 1:N-1
            x_left  = x_nodes(i);     % x_{i-1}
            x_right = x_nodes(i+2);   % x_{i+1}
            b_vec(i) = integral(@(x) f(x) .* phi(x, i, h, a), x_left, x_right);
        end

        % Assemble stiffness matrix A
        for i = 1:N-1
            x_left = x_nodes(i);      % x_{i-1}
            x_mid  = x_nodes(i+1);    % x_i
            x_right = x_nodes(i+2);   % x_{i+1}

            % Diagonal entry A(i,i): ∫_{x_{i-1}}^{x_{i+1}} [...]
            A(i,i) = integral(@(x) ...
                P(x) .* phi_prime(x, i, h, a) .* phi_prime(x, i, h, a) + ...
                Q(x) .* phi      (x, i, h, a) .* phi      (x, i, h, a), ...
                x_left, x_right);

            % Lower diagonal A(i,i-1) and symmetric A(i-1,i)
            if i > 1
                val = integral(@(x) ...
                    P(x) .* phi_prime(x, i,   h, a) .* phi_prime(x, i-1, h, a) + ...
                    Q(x) .* phi      (x, i,   h, a) .* phi      (x, i-1, h, a), ...
                    x_left, x_mid);
                A(i, i-1) = val;
                A(i-1, i) = val;
            end

            % Upper diagonal A(i,i+1) and symmetric A(i+1,i)
            if i < N-1
                val = integral(@(x) ...
                    P(x) .* phi_prime(x, i,   h, a) .* phi_prime(x, i+1, h, a) + ...
                    Q(x) .* phi      (x, i,   h, a) .* phi      (x, i+1, h, a), ...
                    x_mid, x_right);
                A(i, i+1) = val;
                A(i+1, i) = val;
            end
        end

        % Solve linear system A * c = b
        c = A \ b_vec;          % interior nodal values v(x_i), i=1,...,N-1

        % Build full nodal vector including boundary zeros
        v_nodes = zeros(N+1, 1);
        v_nodes(2:N) = c;

        % Store solution for this N
        sol.x = x_nodes(:);
        sol.v = v_nodes;
        solutions{idx} = sol;
    end
end


function val = phi(x, i, h, a)
    % Piecewise linear hat function φ_i with support on [x_{i-1}, x_{i+1}]
    xi_minus = a + (i-1)*h;   % x_{i-1}
    xi       = a + i*h;       % x_i
    xi_plus  = a + (i+1)*h;   % x_{i+1}

    val = zeros(size(x));

    left  = (x >= xi_minus) & (x <= xi);
    right = (x >= xi)       & (x <= xi_plus);

    val(left)  = (x(left)  - xi_minus) / h;
    val(right) = (xi_plus  - x(right)) / h;
end


function val = phi_prime(x, i, h, a)
    % Derivative of φ_i: piecewise constant on [x_{i-1},x_i] and [x_i,x_{i+1}]
    xi_minus = a + (i-1)*h;
    xi       = a + i*h;
    xi_plus  = a + (i+1)*h;

    val = zeros(size(x));

    left  = (x >= xi_minus) & (x <= xi);
    right = (x >= xi)       & (x <= xi_plus);

    val(left)  =  1/h;
    val(right) = -1/h;
end
