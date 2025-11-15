function res_surface = swr_error_surface_versusgt_1D(N, a, M, T, c, dh, dt, gamma, nu, k, THETA1, THETA2)
    % SWR_RESIDUAL_SURFACE - Parameter sweep over theta1, theta2
    % Input: explicit parameters + THETA1, THETA2 grids
    % Output: final residuals for each parameter combination

    Lx = N*a + M;  % Calculate domain length

    % --- Modal initial condition
    m_mode  = 1;
    k0      = m_mode*pi/Lx;
    A0      = 1.0;      % Initial displacement amplitude
    v0amp   = 0.0;      % Initial velocity amplitude
    
    u0  = @(x) A0 * sin(k0*x);
    v0  = @(x) v0amp * sin(k0*x);


    % Ground-truth analytical solution
    Nx = round(Lx/dh) + 1;
    Nt = round(T/dt) + 1;
    x_grid = linspace(0,Lx,Nx);
    t_grid = linspace(0,T,Nt);
    u_gt = analytic_solution_single_mode(x_grid, t_grid, c, gamma, nu, ...
                                         k0, A0, v0amp);

    
    % Get dimensions from input matrices
    [num_rows, num_cols] = size(THETA1);
    res_surface = zeros(num_rows, num_cols);
    
    total_iterations = num_rows * num_cols;
    current_iteration = 0;

    u_init = rand(Nx,Nt);
    % u_init = zeros(Nx,Nt);
    
    % Main parameter sweep
    for ii = 1:num_cols
        for jj = 1:num_rows
            current_iteration = current_iteration + 1;
            theta1 = THETA1(jj, ii);
            theta2 = THETA2(jj, ii);
            
            if mod(current_iteration, 10) == 1 || current_iteration == total_iterations
                fprintf('Progress: %.1f%% (%d/%d) - theta1=%.5g, theta2=%.5g\n', ...
                    current_iteration/total_iterations*100, current_iteration, total_iterations, theta1, theta2);
            end
            
            % Run SWR for this parameter combination
            [~, final_res, ~] = run_swr_1D(u0, v0, N, a, M, T, c, dh, dt, gamma, nu, theta1, theta2, k, u_init, u_gt);
            
            % Store final residual
            res_surface(jj, ii) = final_res;
            
            % fprintf('  Final residual after %d iterations: %.3e\n', k, final_res);
        end
    end
end