function [N_values, h_values, solutions] = Finite_difference(P, Q, R, a, b, alpha, beta, Nmin, Nmax)
% Finite difference solver for the BVP:
%   u''(x) + P(x) u'(x) + Q(x) u(x) = R(x),   a < x < b
%   u(a) = alpha,  u(b) = beta
%
% using second-order central differences and the scheme:
%   -1 - (h/2)P_i,   2 + h^2 Q_i,   -1 + (h/2)P_i.
%
% Inputs:
%   P, Q, R   : function handles @(x) ...
%   a, b      : interval endpoints
%   alpha     : boundary value u(a)
%   beta      : boundary value u(b)
%   Nmin      : minimum number of subintervals (integer >= 2)
%   Nmax      : maximum number of subintervals
%
% Outputs:
%   N_values  : vector of N used
%   h_values  : corresponding mesh sizes h = (b-a)/N
%   solutions : cell array; solutions{k} has fields:
%                 .x = node coordinates (N_k+1 x 1)
%                 .u = nodal values of u (N_k+1 x 1)

    if nargin < 9
        error('Usage: bvp_fd_no_exact(P,Q,R,a,b,alpha,beta,Nmin,Nmax)');
    end

    N_values = Nmin:Nmax;
    h_values = (b - a) ./ N_values;
    solutions = cell(length(N_values), 1);

    for idx = 1:length(N_values)
        N = N_values(idx);
        h = (b - a) / N;

        % nodes x_0,...,x_N
        x_nodes = linspace(a, b, N+1)';

        % System size: N-1 interior points
        A     = zeros(N-1, N-1);
        b_vec = zeros(N-1, 1);

        % Build the tridiagonal system
        for i = 1:N-1
            x_i = a + i*h;
            Pi  = P(x_i);
            Qi  = Q(x_i);
            Ri  = R(x_i);

            if i == 1
                % left neighbor involves u_0 = alpha
                A(i,i)   = 2 + h^2 * Qi;
                A(i,i+1) = -1 + (h/2) * Pi;
                b_vec(i) = -h^2 * Ri + alpha * (1 + (h/2) * Pi);

            elseif i == N-1
                % right neighbor involves u_N = beta
                A(i,i-1) = -1 - (h/2) * Pi;
                A(i,i)   = 2 + h^2 * Qi;
                b_vec(i) = -h^2 * Ri + beta * (1 - (h/2) * Pi);

            else
                % interior point
                A(i,i-1) = -1 - (h/2) * Pi;
                A(i,i)   = 2 + h^2 * Qi;
                A(i,i+1) = -1 + (h/2) * Pi;
                b_vec(i) = -h^2 * Ri;
            end
        end

        % Solve for interior values
        u_inner = A \ b_vec;

        % Assemble full solution including boundary values
        u_approx = [alpha; u_inner; beta];

        % Store solution
        sol.x = x_nodes;
        sol.u = u_approx;
        solutions{idx} = sol;
    end
end
