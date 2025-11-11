function [ud, final_res, res_history] = run_swr_1D(u0, v0, N, a_val, M, T, c, dh, dt, gamma, nu, theta1, theta2, k, u_init, u_ref)
    % RUN_SWR - Execute Waveform Relaxation with fixed iterations
    % Input: explicit parameters
    % Output: 
    %   ud - final SWR solution
    %   final_res - final residual achieved
    %   res_history - residual history array
    
    % Derived parameters
    aj = @(j) a_val*(j-1);
    bj = @(j) aj(j+1) + M;
    Lx = bj(N);
    dx = dh;
    
    ajd = @(j) round(aj(j)/dh) + 1;
    bjd = @(j) round(bj(j)/dh) + 1;
    Nxj = bjd(1);
    
    % Derived parameters
    Nx = round(Lx/dh) + 1;
    Nt = round(T/dt) + 1;
    CFL = c * dt / dh;

    x_axis = linspace(0, Lx, Nx);
    t_axis = linspace(0, T, Nt);

    %% Assembling FDTD matrices for subdomains
    % scheme B coefficients
    a1 = -0.5*nu * dt/dh^2;
    a3 = 1 + 0.5*gamma*dt + nu * dt/dh^2;
    b1 = CFL^2;
    b3 = 2 * (1 - CFL^2);
    c1 = -0.5*nu * dt/dh^2;
    c3 = -(1 - 0.5*gamma*dt - nu * dt/dh^2);

    Aj = sparse(Nxj);
    Bj = sparse(Nxj);
    Cj = sparse(Nxj);

    % First row
    i = 1;
    Aj(i,i:i+1) = [a3, a1];
    Bj(i,i:i+1) = [b3, b1];
    Cj(i,i:i+1) = [c3, c1];

    % Middle rows
    for i = 2:Nxj-1
        Aj(i,i-1:i+1) = [a1, a3, a1];
        Bj(i,i-1:i+1) = [b1, b3, b1];
        Cj(i,i-1:i+1) = [c1, c3, c1];
    end

    % Last row
    i = Nxj;
    Aj(i,i-1:i) = [a1, a3];
    Bj(i,i-1:i) = [b1, b3];
    Cj(i,i-1:i) = [c1, c3];

    %% Initialize solution arrays
    ud = u_init;
    % Initialize ujnew from u_init for each subdomain
    ujnew = zeros(N, Nx, Nt);
    for j = 1:N
        ujnew(j, :, :) = u_init;
    end

    u0_val = u_ref(:,1)';
    u1_val = u_ref(:,2)';

    % Set initial conditions from reference solution
    for j = 1:N
        ujnew(j, ajd(j):bjd(j), 1:2) = [u0_val(ajd(j):bjd(j));u1_val(ajd(j):bjd(j))]';
    end

    % Robin coefficients
    A0_0 = 3/dx + 3*theta1/dt + 2*theta2;
    A0_1 = -4/dx;
    A0_2 = 1/dx;
    AL_N = A0_0;
    AL_N1 = A0_1;
    AL_N2 = A0_2;

    % Create modified matrices for each subdomain with Robin conditions
    Aj_robin = cell(N, 1);
    for j = 1:N
        Aj_robin{j} = Aj;
        
        % Left boundary
        if j == 1
            % Global left boundary - Dirichlet
            Aj_robin{j}(1,:) = 0;
            Aj_robin{j}(1,1) = 1;
        else
            % Interface - Robin
            Aj_robin{j}(1,:) = 0;
            Aj_robin{j}(1,1:3) = [A0_0, A0_1, A0_2];
        end
        
        % Right boundary
        if j == N
            % Global right boundary - Dirichlet
            Aj_robin{j}(end,:) = 0;
            Aj_robin{j}(end,end) = 1;
        else
            % Interface - Robin
            Aj_robin{j}(end,:) = 0;
            Aj_robin{j}(end,end-2:end) = [AL_N2, AL_N1, AL_N];
        end
    end

    %% SWR iteration with fixed number of iterations
    red_order = 1:2:N;
    black_order = 2:2:N;
    order = [red_order, black_order];
    
    res_history = zeros(1, N*ceil(k/N));  % Pre-allocate residual history array

    iter_k = 0;
    ckN = ceil(k/N);
    
    for iter = 1:ckN
        for idx = 1:length(order)
            iter_k = iter_k + 1;

            j = order(idx);
            
            % Get neighbor solutions
            if j > 1
                ujleft = squeeze(ujnew(j-1, :, :));
            end
            if j < N
                ujright = squeeze(ujnew(j+1, :, :));
            end
            
            % Time stepping for subdomain j
            for n = 2:Nt-1
                % Extract subdomain solution
                subdomain_indices = ajd(j):bjd(j);
                uj_sub = squeeze(ujnew(j, subdomain_indices, n))';
                uj_sub_prev = squeeze(ujnew(j, subdomain_indices, n-1))';
                
                % Compute RHS for interior points
                RHS = Bj * uj_sub + Cj * uj_sub_prev;
                
                % Apply boundary conditions
                if j == 1
                    % Left boundary: Dirichlet
                    RHS(1) = 0;
                else
                    % Left interface: Robin
                    var_aj = ajd(j);
                    RHS(1) = (4*theta1/dt) * ujnew(j, var_aj, n) - (theta1/dt) * ujnew(j, var_aj, n-1) ...
                                + [A0_0, A0_1, A0_2] * ujleft(var_aj:var_aj+2,n+1) ...
                                - (4*theta1/dt) * ujleft(var_aj, n) + (theta1/dt) * ujleft(var_aj, n-1);
                end
                
                if j == N
                    % Right boundary: Dirichlet
                    RHS(end) = 0;
                else
                    % Right interface: Robin
                    var_bj = bjd(j);
                    RHS(end) = (4*theta1/dt) * ujnew(j, var_bj, n) - (theta1/dt) * ujnew(j, var_bj, n-1) ...
                                + [AL_N2, AL_N1, AL_N] * ujright(var_bj-2:var_bj,n+1) ...
                                - (4*theta1/dt) * ujright(var_bj, n) + (theta1/dt) * ujright(var_bj, n-1);
                end
                
                % Solve for next time step
                ujnew(j, subdomain_indices, n+1) = Aj_robin{j} \ RHS;
            end
        
            % Update global solution
            ud(ajd(j):bjd(j), :) = ujnew(j, ajd(j):bjd(j), :);

            % figure()
            % mesh(x_axis, t_axis, abs(ud - u_ref)' / max(abs(u_ref), [], 'all'))
            % xlabel('space')
            % ylabel('time')
            % set(gca, 'YDir', 'normal')
            % title_str = sprintf('Iteration %d', iter_k);
            % title(title_str)

            % Compute and store residual after each complete iteration
            res = max(abs(ud - u_ref), [], 'all') / max(abs(u_ref), [], 'all');
            % res = max(abs(ud(:,2:end-1) - u_ref(:,2:end-1)), [], 'all') / max(abs(u_ref(:,2:end-1)), [], 'all');
            res_history(iter_k) = res;
        end
    end
    
    % Final residual is the last computed residual
    final_res = res_history(end);
end