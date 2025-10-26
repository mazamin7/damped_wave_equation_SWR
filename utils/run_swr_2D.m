function [ud, final_res, res_history] = run_swr_2D(u0, v0, N, a_val, M, Ly, T, c, dh, dt, gamma, nu, theta1, theta2, k, u_init, u_ref)
% 2D SWR with Robin transmission (Λ = p ∂t + q).
% Optimized version with precomputed boundary patterns and efficient interface handling.

% ---------------- geometry and grids ----------------
aj  = @(j) a_val*(j-1);
bj  = @(j) aj(j+1) + M;
Lx  = bj(N);

Nx  = round(Lx/dh) + 1;
Ny  = round(Ly/dh) + 1;

ajd = @(j) round(aj(j)/dh) + 1;   % global x-index for x=a_j
bjd = @(j) round(bj(j)/dh) + 1;   % global x-index for x=b_j
Nxj = bjd(1);                     % local x-size

Nt  = floor(T/dt);
CFL = c*dt/dh;

x_axis = linspace(0,Lx,Nx);
y_axis = linspace(0,Ly,Ny);

% Local indexing function
ix  = @(ii,jj) (jj-1)*Nxj + ii;

% --------------- scheme-B coefficients ---------------
% 2D Laplacian
ex = ones(Nxj,1); ey = ones(Ny,1);
Tx = spdiags([ex -2*ex ex], -1:1, Nxj, Nxj) / dh^2;
Ty = spdiags([ey -2*ey ey], -1:1, Ny, Ny) / dh^2;
Ld = kron(speye(Ny), Tx) + kron(Ty, speye(Nxj));

% centered 3-level with damping
dt2 = dt^2;
mu0 = 1/dt2 + gamma/(2*dt);
mu1 = 1/dt2 - gamma/(2*dt);
beta = nu/(2*dt);

Aj = mu0*speye(Nxj*Ny) - beta*Ld;
Bj = (2/dt2)*speye(Nxj*Ny) + c^2*Ld;
Cj = -mu1*speye(Nxj*Ny) - beta*Ld;

% % -------------- SWR arrays (ud and ujnew) --------------
% ud    = u_init;                             % global SWR solution
% ujnew = zeros(N, Nx, Ny, Nt);               % per-subdomain fields
% for j = 1:N, ujnew(j,:,:,:) = u_init; end

% -------------- SWR arrays (ud and ujnew) --------------
ud = u_init;                                  % Nx x Ny x Nt
ujnew = repmat(reshape(u_init,[1,size(u_init)]), [N,1,1,1]);  % N x Nx x Ny x Nt

% initial conditions from u0,v0
% U0 = zeros(Nx,Ny); U1 = zeros(Nx,Ny);
% for jj = 1:Ny
%     U0(:,jj) = u0(x_axis, y_axis(jj));
%     U1(:,jj) = u0(x_axis, y_axis(jj)) + dt * v0(x_axis, y_axis(jj));
% end

[Xg,Yg] = ndgrid(x_axis, y_axis);           % Nx-by-Ny
U0 = u0(Xg, Yg);                             % Nx-by-Ny
U1 = u0(Xg, Yg) + dt * v0(Xg, Yg);           % Nx-by-Ny

ud(:,:,1) = U0; ud(:,:,2) = U1;
for j = 1:N
    ujnew(j,ajd(j):bjd(j),:,1) = U0(ajd(j):bjd(j),:);
    ujnew(j,ajd(j):bjd(j),:,2) = U1(ajd(j):bjd(j),:);
end

% -------------- Robin coefficients (same on both sides) --------------
p = theta1; q = theta2;
A0_0 = 3/dh + 3*p/dt + 2*q;    % left rows coef at i=1
A0_1 = -4/dh;
A0_2 =  1/dh;
AL_N = A0_0;  AL_N1 = A0_1;  AL_N2 = A0_2; % right rows at i=Nxj

% ------- PRE-COMPUTE ALL BOUNDARY AND INTERFACE INDICES -------
boundary_data = cell(N,1);  % Store all boundary/interface data per subdomain

for j = 1:N
    data = struct();
    
    % Precompute all boundary indices
    % Dirichlet boundaries (y=0 and y=Ly)
    data.y0_rows = 1:Nxj;
    data.yL_rows = (Ny-1)*Nxj+1 : Nxj*Ny;
    data.all_dirichlet_rows = [data.y0_rows, data.yL_rows];
    
    % Interface boundaries
    data.xa_rows = 1:Nxj:Nxj*Ny;  % x = a_j (left boundary)
    data.xb_rows = Nxj:Nxj:Nxj*Ny; % x = b_j (right boundary)
    
    % Remove corners from interface rows (already handled by Dirichlet)
    data.xa_rows_no_corners = setdiff(data.xa_rows, [1, (Ny-1)*Nxj+1]);
    data.xb_rows_no_corners = setdiff(data.xb_rows, [Nxj, Nxj*Ny]);
    
    % Create boundary modification patterns
    n_total = Nxj*Ny;
    
    % Dirichlet pattern (identity rows for y boundaries)
    data.dirichlet_pattern = speye(n_total);
    data.dirichlet_pattern = data.dirichlet_pattern(data.all_dirichlet_rows, :);
    
    % Left interface pattern (for Robin conditions)
    if j > 1
        data.left_pattern = sparse(length(data.xa_rows_no_corners), n_total);
        for idx = 1:length(data.xa_rows_no_corners)
            row_idx = data.xa_rows_no_corners(idx);
            jj = floor((row_idx-1)/Nxj) + 1;  % Recover jj from row index
            data.left_pattern(idx, [ix(1,jj), ix(2,jj), ix(3,jj)]) = [A0_0, A0_1, A0_2];
        end
    else
        data.left_pattern = speye(length(data.xa_rows)); % Dirichlet at global boundary
    end
    
    % Right interface pattern (for Robin conditions)
    if j < N
        data.right_pattern = sparse(length(data.xb_rows_no_corners), n_total);
        for idx = 1:length(data.xb_rows_no_corners)
            row_idx = data.xb_rows_no_corners(idx);
            jj = floor((row_idx-1)/Nxj) + 1;  % Recover jj from row index
            data.right_pattern(idx, [ix(Nxj-2,jj), ix(Nxj-1,jj), ix(Nxj,jj)]) = [AL_N2, AL_N1, AL_N];
        end
    else
        data.right_pattern = speye(length(data.xb_rows)); % Dirichlet at global boundary
    end
    
    % Precompute global indices for neighbor access
    if j > 1
        data.left_neighbor_global = ajd(j);  % Same physical point as right boundary of j-1
    end
    if j < N
        data.right_neighbor_global = bjd(j); % Same physical point as left boundary of j+1
    end
    
    boundary_data{j} = data;
end

% ------- PRE-COMPUTE MODIFIED MATRICES WITH FACTORIZATION -------
Aj_modified = cell(N,1);

for j = 1:N
    data = boundary_data{j};
    A_mod = Aj;  % Start with original matrix
    
    % Apply Dirichlet conditions (y boundaries)
    A_mod(data.all_dirichlet_rows, :) = 0;
    for idx = 1:length(data.all_dirichlet_rows)
        row = data.all_dirichlet_rows(idx);
        A_mod(row, row) = 1;
    end
    
    % Apply left boundary conditions
    if j == 1
        % Global left boundary - Dirichlet
        for idx = 1:length(data.xa_rows)
            row = data.xa_rows(idx);
            A_mod(row, :) = 0;
            A_mod(row, row) = 1;
        end
    else
        % Interface - Robin conditions
        for idx = 1:length(data.xa_rows_no_corners)
            row = data.xa_rows_no_corners(idx);
            jj = floor((row-1)/Nxj) + 1;
            A_mod(row, :) = 0;
            A_mod(row, [ix(1,jj), ix(2,jj), ix(3,jj)]) = [A0_0, A0_1, A0_2];
        end
    end
    
    % Apply right boundary conditions
    if j == N
        % Global right boundary - Dirichlet
        for idx = 1:length(data.xb_rows)
            row = data.xb_rows(idx);
            A_mod(row, :) = 0;
            A_mod(row, row) = 1;
        end
    else
        % Interface - Robin conditions
        for idx = 1:length(data.xb_rows_no_corners)
            row = data.xb_rows_no_corners(idx);
            jj = floor((row-1)/Nxj) + 1;
            A_mod(row, :) = 0;
            A_mod(row, [ix(Nxj-2,jj), ix(Nxj-1,jj), ix(Nxj,jj)]) = [AL_N2, AL_N1, AL_N];
        end
    end
    
    % Factorize for efficient solving
    Aj_modified{j} = decomposition(A_mod);
end

% ---------------- SWR iteration (multiplicative) ----------------
red_order   = 1:2:N;
black_order = 2:2:N;
order       = [red_order, black_order];

ckN = ceil(k/N);
res_history = zeros(1, N*ckN);
iter_k = 0;

% Precompute time coefficients for Robin conditions
time_coeff1 = (4*p/dt);
time_coeff2 = (p/dt);

for iter = 1:ckN
    for idx = 1:length(order)
        iter_k = iter_k + 1;
        j = order(idx);

        % Get neighbor solutions
        if j > 1, ujleft  = squeeze(ujnew(j-1,:,:,:)); end
        if j < N, ujright = squeeze(ujnew(j+1,:,:,:)); end

        % Get boundary data for current subdomain
        data = boundary_data{j};
        
        % Time stepping on subdomain j
        for n = 2:Nt-1
            % Local slices at time n, n-1
            un   = reshape(squeeze(ujnew(j, ajd(j):bjd(j), :, n)), Nxj*Ny, 1);
            unm1 = reshape(squeeze(ujnew(j, ajd(j):bjd(j), :, n-1)), Nxj*Ny, 1);

            RHS = Bj*un + Cj*unm1;

            % Apply Dirichlet conditions (y boundaries)
            RHS(data.all_dirichlet_rows) = 0;
            
            % Apply left boundary conditions
            if j == 1
                % Global left boundary - Dirichlet
                RHS(data.xa_rows) = 0;
            else
                % Interface - Robin conditions
                for idx = 1:length(data.xa_rows_no_corners)
                    row = data.xa_rows_no_corners(idx);
                    jj = floor((row-1)/Nxj) + 1;
                    
                    % Local history term
                    rloc = time_coeff1 * un(row) - time_coeff2 * unm1(row);
                    
                    % Neighbor trace using forward stencil
                    ia = data.left_neighbor_global;
                    neighbor_term = A0_0 * ujleft(ia, jj, n+1) ...
                                  + A0_1 * ujleft(ia+1, jj, n+1) ...
                                  + A0_2 * ujleft(ia+2, jj, n+1) ...
                                  - time_coeff1 * ujleft(ia, jj, n) ...
                                  + time_coeff2 * ujleft(ia, jj, n-1);
                    
                    RHS(row) = rloc + neighbor_term;
                end
            end
            
            % Apply right boundary conditions
            if j == N
                % Global right boundary - Dirichlet
                RHS(data.xb_rows) = 0;
            else
                % Interface - Robin conditions
                for idx = 1:length(data.xb_rows_no_corners)
                    row = data.xb_rows_no_corners(idx);
                    jj = floor((row-1)/Nxj) + 1;
                    
                    % Local history term
                    rloc = time_coeff1 * un(row) - time_coeff2 * unm1(row);
                    
                    % Neighbor trace using backward stencil
                    ib = data.right_neighbor_global;
                    neighbor_term = AL_N2 * ujright(ib-2, jj, n+1) ...
                                  + AL_N1 * ujright(ib-1, jj, n+1) ...
                                  + AL_N  * ujright(ib, jj, n+1) ...
                                  - time_coeff1 * ujright(ib, jj, n) ...
                                  + time_coeff2 * ujright(ib, jj, n-1);
                    
                    RHS(row) = rloc + neighbor_term;
                end
            end

            % Solve for n+1
            unp1 = Aj_modified{j} \ RHS;
            ujnew(j, ajd(j):bjd(j), :, n+1) = reshape(unp1, Nxj, Ny);
        end

        % Update global solution
        ud(ajd(j):bjd(j),:,:) = ujnew(j, ajd(j):bjd(j), :, :);

        % Compute residual
        res = max(abs(ud - u_ref), [], 'all') / max(abs(u_ref), [], 'all');
        res_history(iter_k) = res;
        
        % fprintf('Iter %d, Subdomain %d, Residual: %.2e\n', iter_k, j, res);
    end
end

final_res = res_history(iter_k);
end