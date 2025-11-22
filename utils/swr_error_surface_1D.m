function res_surface = swr_error_surface_1D(N, a, M, T, c, dh, dt, gamma, nu, k, THETA1, THETA2, u0, v0, u_init, u_ref)
    % SWR_ERROR_SURFACE_1D - Parameter sweep over theta1, theta2
    % Input: 
    %   Standard params (N...k)
    %   THETA1, THETA2 : Grids of parameters
    %   u0, v0         : Function handles for initial conditions
    %   u_init         : The specific random initial guess matrix (Nx x Nt)
    %   u_ref          : The pre-computed FDTD reference solution
    %
    % Output: final residuals for each parameter combination

    % Get dimensions from input matrices
    [num_rows, num_cols] = size(THETA1);
    res_surface = zeros(num_rows, num_cols);
    
    total_iterations = num_rows * num_cols;
    current_iteration = 0;
    
    % Error at iteration 0 (same for all parameter pairs since u_init and u_ref are fixed)
    err0_test = max(abs(u_init(:) - u_ref(:))) / max(abs(u_ref(:)));
    
    % Number of iterations for the test run (to calculate amplification factor)
    k_test = 1;

    fprintf('Starting Grid Search (%d points)...\n', total_iterations);

    % Main parameter sweep
    for ii = 1:num_cols
        for jj = 1:num_rows
            current_iteration = current_iteration + 1;
            theta1 = THETA1(jj, ii);
            theta2 = THETA2(jj, ii);
            
            % Progress print
            if mod(current_iteration, 50) == 1 || current_iteration == total_iterations
                fprintf('  Grid Progress: %.1f%% (%d/%d)\n', ...
                    current_iteration/total_iterations*100, current_iteration, total_iterations);
            end
            
            % ------------------------------------------------------------
            % TEST RUN: 1 SWR iteration to get amplification factor F
            % ------------------------------------------------------------
            [~, ~, res_history_test] = run_swr_1D( ...
                u0, v0, N, a, M, T, c, dh, dt, gamma, nu, ...
                theta1, theta2, k_test, u_init, u_ref);
            
            % Error after first SWR iteration in the test run
            E1_test = res_history_test(1);
            
            % Amplification factor F = E1 / E0
            F = E1_test / err0_test;
            
            % Guard against degenerate case
            if F == 0
                F = 1;
            end
            
            % ------------------------------------------------------------
            % REAL RUN: rescaled initial condition u_init / F
            % ------------------------------------------------------------
            u_init_scaled = u_init / F;
            
            [~, final_res, ~] = run_swr_1D( ...
                u0, v0, N, a, M, T, c, dh, dt, gamma, nu, ...
                theta1, theta2, k, u_init_scaled, u_ref);
            
            % Store final residual
            res_surface(jj, ii) = final_res;
        end
    end
end