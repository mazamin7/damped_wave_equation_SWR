function u = run_fdtd_1D(u0, v0, Lx, T, c, dh, dt, gamma, nu)
    % RUN_FDTD_REFERENCE - Compute reference solution using FDTD
    % Input: explicit parameters
    % Output: u - reference solution matrix
    
    % Derived parameters
    Nx = round(Lx / dh) + 1;
    Nt = floor(T / dt);
    CFL = c * dt / dh;

    x_axis = linspace(0, Lx, Nx);
    t_axis = linspace(0, T, Nt);

    %% Assembling FDTD matrices
    % scheme B
    a1 = -0.5*nu * dt/dh^2;
    a3 = 1 + 0.5*gamma*dt + nu * dt/dh^2;
    b1 = CFL^2;
    b3 = 2 * (1 - CFL^2);
    c1 = -0.5*nu * dt/dh^2;
    c3 = -(1 - 0.5*gamma*dt - nu * dt/dh^2);

    A = sparse(Nx);
    B = sparse(Nx);
    C = sparse(Nx);

    i = 1;
    A(i,i:i+1) = [a3, a1];
    B(i,i:i+1) = [b3, b1];
    C(i,i:i+1) = [c3, c1];

    for i = 2:Nx-1
        A(i,i-1:i+1) = [a1, a3, a1];
        B(i,i-1:i+1) = [b1, b3, b1];
        C(i,i-1:i+1) = [c1, c3, c1];
    end

    i = Nx;
    A(i,i-1:i) = [a1, a3];
    B(i,i-1:i) = [b1, b3];
    C(i,i-1:i) = [c1, c3];

    % Dirichlet boundary conditions
    A_mod = A;
    A_mod(1,:) = 0; A_mod(1,1) = 1;
    A_mod(Nx,:) = 0; A_mod(Nx,Nx) = 1;

    %% Pure FDTD reference solution
    u = zeros(Nx, Nt);
    u(:,1) = u0(x_axis);
    u(:,2) = u0(x_axis) + dt * v0(x_axis);

    for n = 2:Nt-1
        RHS = B*u(:,n) + C*u(:,n-1);
        
        % Dirichlet boundary conditions
        RHS(1) = 0;
        RHS(Nx) = 0;
        
        u(:,n+1) = A_mod \ RHS;
    end
end