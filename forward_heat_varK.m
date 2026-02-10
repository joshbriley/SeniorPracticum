function [x, t, U, full_A, rhs, biot] = forward_heat_varK(k0, k1, N, dt, T, u0fun, T_infy, h, big_U, L)
    if k0 <= abs(k1) 
	    error('Need k0 > |k1| so that k(x) stays positive.'); 
    end
    rho = 8.9; heat_capacity = 385; % Grid discretization 
    delta_x = 1 / N; x = (0:N)' * delta_x; % points (N+1)x1 
    xi = x(2:N); % interior points (N-1)x1 
    
    % Time 
    Nt = round(T / dt); 
    t = (0:Nt) * dt; 
    
    % Initial condition
    U = zeros(N+1, Nt+1); 
    U(2:N,1) = u0fun(xi); 
    U = (U - T_infy)./big_U; 
    
    % Midpoint evaluation 
    xface = (x(1:N) + x(2:N+1)) / 2; % x_{i+1/2}, i=0..N-1 (N x 1) 
    alphaface = k0 + k1 .* sin(pi .* xface); % kappa_{i+1/2} (N x 1) 
    kface = rho.*heat_capacity.*alphaface; 
    biot = (h.*L)./kface; 
    
    % Build tridiagonal matrix A = I - dt*L on interior unknowns 
    % For interior i=1..N-1: % lower diag: -dt/delta_x^2 * k_{i-1/2} 
    % main diag : 1 + dt/delta_x^2*(k_{i-1/2}+k_{i+1/2}) 
    % upper diag: -dt/delta_x^2 * k_{i+1/2} 
    
    lambda = (dt / delta_x^2); 
    a_imhalf = alphaface(1:N-1); 
    a_iphalf = alphaface(2:N); 
    main = 1 + lambda .* (a_imhalf + a_iphalf); 
    lower = -lambda .* a_imhalf; 
    upper = -lambda .* a_iphalf; 
    A = spdiags([lower, main, upper], [-1, 0, 1], N-1, N-1); 
    
    % Define additional row for Robin BCs 
    robin_row = zeros(1, N); 
    robin_row(1) = biot(1)*delta_x - alphaface(1); 
    robin_row(2) = -alphaface(1); 
    
    % Define additional column for Robin BCs 
    robin_col = zeros(N-1, 1); 
    robin_col(1) = lower(1); 
    
    % Add new row/column to A matrix 
    A = [robin_col, A]; 
    A = [robin_row; A]; 
    A(end+1, end+1) = 1; 
    A(end, end-1) = -1; 
    A(end-1, end) = upper(1); 
    full_A = full(A); 

    % rhs = zeros(N+1, 1);
    
    % Time stepping: A * u^{n+1} = u^n 
    for n = 1:Nt 
        rhs = U(1:N, n); 
        rhs(1) = biot(1); % Additional parameter for Robin 
        rhs(end+1) = 0;
        U(1:N+1, n+1) = A \ rhs; % solve tridiagonal sparse system 
    end 
    U = big_U.*U + T_infy; 
end
