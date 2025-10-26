function [ud, final_res, res_history] = run_swr_2D(u0, v0, N, a, M, T, c, dh, dt, gamma, nu, theta1, theta2, k, u_init, u_ref)
%RUN_SWR_2D - Run 2D Schwarz Waveform Relaxation
% Inputs:
%   u0, v0 - function handles for initial conditions u0(x,y), v0(x,y)
%   N - number of subdomains
%   a - length of each subdomain (except last has length a + M)
%   M - overlap size
%   T - total simulation time
%   c - wave speed
%   dh - spatial step size
%   gamma - damping coefficient for ∂_t term
%   nu - damping coefficient for ∂_t∂_xx term
%   theta1, theta2 - transmission condition parameters (not used in current scheme)
%   k - number of SWR iterations
%   u_init - initial guess for SWR
%   u_ref - reference solution for residuals (Nx, Ny, Nt)
% Outputs:
%   ud - final SWR solution (Nx, Ny, Nt)
%   final_res - final residual achieved
%   res_history - residual history array

% Domain parameters
aj = @(j) a*(j-1);
bj = @(j) aj(j+1) + M;
Lx = bj(N);
Ly = 5;  % Fixed Ly

ajd = @(j) round(aj(j)/dh) + 1;
bjd = @(j) round(bj(j)/dh) + 1;
Nxj = bjd(1);

% Derived parameters
Nx = round(Lx / dh) + 1;
Ny = round(Ly / dh) + 1;
Nt = floor(T / dt);
CFL = c * dt / dh;

x_axis = linspace(0, Lx, Nx);
y_axis = linspace(0, Ly, Ny);

%% Assembling FDTD matrix for subdomains
% 2D Laplacian on subdomain
exj = ones(Nxj,1); ey = ones(Ny,1);
Txj = spdiags([exj -2*exj exj], -1:1, Nxj, Nxj) / dh^2;
Ty = spdiags([ey -2*ey ey], -1:1, Ny, Ny) / dh^2;
Lj = kron(speye(Ny), Txj) + kron(Ty, speye(Nxj));

% centered 3-level with damping
dt2 = dt^2;
mu0 = 1/dt2 + gamma/(2*dt);
mu1 = 1/dt2 - gamma/(2*dt);
beta = nu/(2*dt);

Aj = mu0*speye(Nxj*Ny) - beta*Lj;
Bj = (2/dt2)*speye(Nxj*Ny) + c^2*Lj;
Cj = -mu1*speye(Nxj*Ny) - beta*Lj;

% Precompute boundary indices for subdomains
xa_indices_j = 1:Nxj:Nxj*Ny;           % x = aj boundary  
xb_indices_j = Nxj:Nxj:Nxj*Ny;         % x = bj boundary
y0_indices_j = 1:Nxj;                  % y = 0 boundary
yL_indices_j = (Ny-1)*Nxj+1:Nxj*Ny;    % y = Ly boundary

% All boundary indices for subdomains
boundary_indices_j = unique([xa_indices_j, xb_indices_j, y0_indices_j, yL_indices_j]);
boundary_pattern_j = speye(Nxj*Ny);
boundary_pattern_j = boundary_pattern_j(boundary_indices_j, :);

%% Initialize arrays
ud = u_init;
% Initialize ujnew from u_init for each subdomain
ujnew = zeros(N, Nx, Ny, Nt);
for j = 1:N
    ujnew(j, :, :, :) = u_init;
end

% Set initial conditions from reference solution
for j = 1:N
    ujnew(j, ajd(j):bjd(j), :, 1:2) = u_ref(ajd(j):bjd(j), :, 1:2);
end

% Sweeping orders
red_order = 1:2:N;
black_order = 2:2:N;
order = [red_order, black_order];

% Pre-allocate residual history array
ckN = ceil(k/N);
res_history = zeros(1, N * ckN);

iter_k = 0;

% Apply boundary conditions
for j = 1:N
    Aj_mod = Aj;
    Aj_mod(boundary_indices_j, :) = boundary_pattern_j;
    Aj_mod = decomposition(Aj_mod);
end

%% SWR iterations
for iter = 1:ckN
    for idx = 1:length(order)
        iter_k = iter_k + 1;
        
        j = order(idx);

        % Get neighbor solutions
        if j > 1
            ujleft = squeeze(ujnew(j-1, :, :, :));
        else
            ujleft = [];
        end

        if j < N
            ujright = squeeze(ujnew(j+1, :, :, :));
        else
            ujright = [];
        end

        % Simulation loop for subdomain j
        for n = 2:Nt-1
            % Assembling RHS
            RHS = Bj * reshape(ujnew(j, ajd(j):bjd(j), :, n), Nxj*Ny, 1) + ...
                  Cj * reshape(ujnew(j, ajd(j):bjd(j), :, n-1), Nxj*Ny, 1);
        
            % Apply boundary conditions
            RHS_mod = RHS;
            RHS_mod(boundary_indices_j) = 0;
            
            % Apply interface conditions for left and right boundaries
            % Left boundary (x = aj)
            if j > 1
                for jj = 0:Ny-1
                    idx_left = 1 + jj*Nxj;
                    RHS_mod(idx_left) = ujleft(ajd(j), jj+1, n+1);
                end
            end
            
            % Right boundary (x = bj)
            if j < N
                for jj = 0:Ny-1
                    idx_right = Nxj + jj*Nxj;
                    RHS_mod(idx_right) = ujright(bjd(j), jj+1, n+1);
                end
            end
        
            % Solve the wave equation for subdomain
            ujnew(j, ajd(j):bjd(j), :, n+1) = reshape(Aj_mod \ RHS_mod, Nxj, Ny);
        end

        % Update global solution
        ud(ajd(j):bjd(j), :, :) = ujnew(j, ajd(j):bjd(j), :, :);

        % Compute and store residual after each complete iteration
        res = max(abs(ud - u_ref), [], 'all') / max(abs(u_ref), [], 'all');
        res_history(iter_k) = res;
    end
end

% Final residual is the last computed residual
final_res = res_history(end);

end