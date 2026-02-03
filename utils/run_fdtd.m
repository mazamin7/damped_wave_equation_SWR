function u = run_fdtd(u0, v0, Lx, T, c, dh, dt, gamma, nu)
    % --- grid-exact assertion for dh, dt ---
    Nx = round(Lx / dh) + 1;
    Nt = round(T  / dt) + 1;
    dh_eff = Lx / (Nx - 1);
    dt_eff = T  / (Nt - 1);

    tol = 1e-12;
    assert(abs(dh_eff - dh) < tol, ...
      'run_fdtd_1D:dh_mismatch', 'Requested dh=%.16g not grid-exact. Use dh=%.16g.', dh, dh_eff);
    assert(abs(dt_eff - dt) < tol, ...
      'run_fdtd_1D:dt_mismatch', 'Requested dt=%.16g not grid-exact. Use dt=%.16g.', dt, dt_eff);
    assert(Nt >= 2, 'run_fdtd_1D:Nt_too_small', 'Nt must be at least 2.');
    CFL = c*dt/dh;
    assert(CFL <= 1+1e-14, 'run_fdtd_1D:CFL', 'CFL=c*dt/dh=%.3g > 1 may be unstable.', CFL);

    % grid
    x  = linspace(0, Lx, Nx).';

    % scheme-B coefficients
    a1 = -0.5*nu * dt/dh^2;
    a3 = 1 + 0.5*gamma*dt + nu * dt/dh^2;
    b1 = CFL^2;
    b3 = 2 * (1 - CFL^2);
    c1 = -0.5*nu * dt/dh^2;
    c3 = -(1 - 0.5*gamma*dt - nu * dt/dh^2);

    A = spalloc(Nx,Nx,3*Nx);
    B = spalloc(Nx,Nx,3*Nx);
    C = spalloc(Nx,Nx,3*Nx);

    for i = 2:Nx-1
        A(i,i-1:i+1) = [a1, a3, a1];
        B(i,i-1:i+1) = [b1, b3, b1];
        C(i,i-1:i+1) = [c1, c3, c1];
    end
    A(1,1:2)        = [a3, a1];      A(Nx,Nx-1:Nx) = [a1, a3];
    B(1,1:2)        = [b3, b1];      B(Nx,Nx-1:Nx) = [b1, b3];
    C(1,1:2)        = [c3, c1];      C(Nx,Nx-1:Nx) = [c1, c3];

    % Dirichlet via row replacement
    A_mod = A;  A_mod(1,:)=0; A_mod(1,1)=1;  A_mod(Nx,:)=0; A_mod(Nx,Nx)=1;

    % second-order startup: u_tt = c^2 u_xx + nu v_xx - gamma v
    u   = zeros(Nx,Nt);

    % FORCE COLUMN SHAPE here
    u0v = u0(x); u0v = u0v(:);
    v0v = v0(x); v0v = v0v(:);

    Dxx = spdiags([ones(Nx,1) -2*ones(Nx,1) ones(Nx,1)], [-1 0 1], Nx, Nx)/dh^2;
    Dxx(1,:)=0; Dxx(Nx,:)=0;

    utt0 = c^2*(Dxx*u0v) + nu*(Dxx*v0v) - gamma*v0v;

    u(:,1)   = u0v;
    if Nt >= 2
        u(:,2) = u0v + dt*v0v + 0.5*dt^2*utt0;
        u(1,1:2)=0; u(end,1:2)=0;
    end

    for n = 2:Nt-1
        RHS = B*u(:,n) + C*u(:,n-1);
        RHS(1)=0; RHS(Nx)=0;
        u(:,n+1) = A_mod \ RHS;
    end
end
