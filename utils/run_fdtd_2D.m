function u = run_fdtd_2D(u0, v0, Lx, Ly, T, c, dh, dt, gamma, nu)
%RUN_FDTD_2D - Run 2D FDTD simulation with damping
% Inputs:
%   u0, v0 - function handles for initial conditions u0(x,y), v0(x,y)
%   Lx, Ly - domain dimensions
%   T - total simulation time
%   c - wave speed
%   dh - spatial step size
%   gamma - damping coefficient for ∂_t term
%   nu - damping coefficient for ∂_t∂_xx term
% Output:
%   u - 3D array (Nx, Ny, Nt) of solution

% Derived parameters
Nx = round(Lx / dh) + 1;
Ny = round(Ly / dh) + 1;
Nt = floor(T / dt);

x_axis = linspace(0, Lx, Nx);
y_axis = linspace(0, Ly, Ny);

%% Assembling FDTD matrices
% 2D Laplacian
ex = ones(Nx,1); ey = ones(Ny,1);
Tx = spdiags([ex -2*ex ex], -1:1, Nx, Nx) / dh^2;
Ty = spdiags([ey -2*ey ey], -1:1, Ny, Ny) / dh^2;
Ld = kron(speye(Ny), Tx) + kron(Ty, speye(Nx));

% centered 3-level with damping
dt2 = dt^2;
mu0 = 1/dt2 + gamma/(2*dt);
mu1 = 1/dt2 - gamma/(2*dt);
beta = nu/(2*dt);

A = mu0*speye(Nx*Ny) - beta*Ld;
B = (2/dt2)*speye(Nx*Ny) + c^2*Ld;
C = -mu1*speye(Nx*Ny) - beta*Ld;

% Precompute boundary indices
x0_indices = 1:Nx:Nx*Ny;
xL_indices = Nx:Nx:Nx*Ny;
y0_indices = 1:Nx;
yL_indices = (Ny-1)*Nx+1:Nx*Ny;
boundary_indices = unique([x0_indices, xL_indices, y0_indices, yL_indices]);
boundary_pattern = speye(Nx*Ny);
boundary_pattern = boundary_pattern(boundary_indices, :);

%% Initialize and run simulation
u = zeros(Nx, Ny, Nt);

% Impose initial conditions
for jj = 1:Ny
    u(:,jj,1) = u0(x_axis, y_axis(jj));
    u(:,jj,2) = u0(x_axis, y_axis(jj)) + dt * v0(x_axis, y_axis(jj));
end

% Apply boundary conditions
A_mod = A;
A_mod(boundary_indices, :) = boundary_pattern;
A_mod = decomposition(A_mod);

% Simulation loop
for n = 2:Nt-1
    % Assembling RHS
    RHS = B*reshape(u(:,:,n), Nx*Ny, 1) + C*reshape(u(:,:,n-1), Nx*Ny, 1);
    
    % Apply boundary conditions
    RHS_mod = RHS;
    RHS_mod(boundary_indices) = 0;
    
    % Solve
    u(:,:,n+1) = reshape(A_mod \ RHS_mod, Nx, Ny);
end

end