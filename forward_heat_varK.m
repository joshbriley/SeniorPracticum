function [x, t, U] = forward_heat_varK(k0, k1, N, dt, T, u0fun, T_infy, h)
% forward_heat_varK
% Solves u_t = (k(x) u_x)_x on x in (0,1), u(0,t)=u(1,t)=0, u(x,0)=u0(x) using
% flux-form finite differences + Backward Euler.
%
% Inputs:
%   k0,k1  : diffusion parameters, k(x)=k0+k1*sin(pi*x), require k0>|k1|
%   N      : number of spatial subintervals (grid has N+1 nodes)
%   dt     : time step
%   T      : final time
%   u0fun  : function handle for initial condition u0(x)
%   T_infy : background temp
%   h      : material parameter
%
% Outputs:
%   x      : grid (N+1 x 1)
%   t      : time vector (1 x Nt+1)
%   U      : solution snapshots at all grid nodes (N+1 x Nt+1)

    if k0 <= abs(k1)
        error('Need k0 > |k1| so that k(x) stays positive.');
    end

    % Grid discretization
    delta_x = 1 / N;
    x = (0:N)' * delta_x; % points (N+1)x1
    xi = x(2:N);    % interior points (N-1)x1

    % Time
    Nt = T / dt;
    t = (0:Nt) * dt;

    % Initial condition (include boundaries as zeros)
    U = zeros(N+1, Nt+1);
    U(2:N,1) = u0fun(xi); 

    % Midpoint evaluation
    xface = (x(1:N) + x(2:N+1)) / 2;   % x_{i+1/2}, i=0..N-1  (N x 1)
    kface = k0 + k1 * sin(pi * xface); % kappa_{i+1/2}        (N x 1)

    % Build tridiagonal matrix A = I - dt*L on interior unknowns
    % For interior i=1..N-1:
    % lower diag: -dt/delta_x^2 * k_{i-1/2}
    % main diag : 1 + dt/delta_x^2*(k_{i-1/2}+k_{i+1/2})
    % upper diag: -dt/delta_x^2 * k_{i+1/2}

    alpha = (dt / delta_x^2);

    k_imhalf = kface(1:N-1); % k_{i-1/2} for interior i=1..N-1  (N-1 x 1)
    k_iphalf = kface(2:N); % k_{i+1/2} for interior i=1..N-1  (N-1 x 1)

    main = 1 + alpha * (k_imhalf + k_iphalf); 
    lower = -alpha * k_imhalf;   
    upper = -alpha * k_iphalf; 

    A = spdiags([lower, main, upper], [-1, 0, 1], N-1, N-1);

    % Define additional row for Robin BCs
    robin_row = zeros(1, N);
    robin_row(1) = k0 - h * delta_x; 
    robin_row(2) = -k0;  

    % Define additional column for Robin BCs
    robin_col = zeros(N-1, 1); 
    robin_col(1) = -alpha;

    % Add new row/column to A matrix
    A = [robin_col, A];
    A = [robin_row; A];

    % Time stepping: A * u^{n+1} = u^n
    for n = 1:Nt
        rhs = U(1:N, n);
        rhs(1) = -h * delta_x * T_infy; % Additional parameter for Robin
        U(1:N, n+1) = A \ rhs; % solve tridiagonal sparse system
    end
end
