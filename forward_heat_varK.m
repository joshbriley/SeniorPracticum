function [x, t, U, full_A, rhs, biot] = forward_heat_varK( ...
    alpha, N, dt, T, u0fun, T_infy, h, big_U, L)

    % --- constants / grid ---
    rho = 8.9;
    Cp  = 385;

    dx = 1/N;
    x  = ((0:N)' * dx)./L;          % (N+1)x1 grid nodes
    xi = x(2:N)./L;               % interior nodes

    % --- time ---
    Nt = round(T/dt);
    t  = (0:Nt) * dt;

    % --- initial condition (dimensionless) ---
    U        = zeros(N+1, Nt+1);
    U(2:N,1) = u0fun(xi);

    % --- Biot ---
    biot = (h * L) ./ alpha;

    % --- Define Lambda and matrix diags ---
    lam   = alpha*(dt / dx^2);
    main  = 1 + 2*lam;
    low   = -lam;
    up    = -lam;

    Aint = spdiags([low, main, up], [-1 0 1], N-1, N-1);

    % --- build full system: (u0, interior..., uN) ---
    % Size: (N+1)x(N+1)
    A = spalloc(N+1, N+1, 3*(N-1) + 10);

    % Put interior block into rows/cols 2..N
    A(2:N, 2:N) = Aint;

    % Robin boundary row at x=0 (your existing formula)
    % robin_row(1) = biot(1)*dx - alphaFace(1);
    % robin_row(2) = -alphaFace(1);
    A(1,1) = biot*dx - 1;
    A(1,2) = 1;

    % Couple Robin equation into first interior equation
    A(2,1) = low(1);

    % Neumann at x=1
    A(end,end)   = 1;
    A(end,end-1) = -1;

    % Couple last interior equation to boundary uN
    A(end-1,end) = up(1);

    full_A = full(A);  % debuging

    % --- time stepping ---
    rhs = zeros(N+1,1);
    for n = 1:Nt
        rhs(:) = U(:,n);

        rhs(1)   = 0;
        rhs(end) = 0;

        U(:,n+1) = A \ rhs;
    end
end

