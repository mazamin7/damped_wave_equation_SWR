function res_surface = swr_residual_surface_2D(N, a, M, Ly, T, c, dh, dt, gamma, nu, k, THETA1, THETA2)
    % SWR_RESIDUAL_SURFACE - Parameter sweep over theta1, theta2
    % Input: explicit parameters + THETA1, THETA2 grids
    % Output: final residuals for each parameter combination
    
    % Compute reference solution once
    fprintf('Computing reference FDTD solution...\n');
    Lx = N*a + M;  % Calculate domain length


    % Initial conditions
    u0 = @(x,y) 0.*x.*y;
    v0 = @(x,y) exp(-sqrt((x-Lx/3).^2+(y-Ly/3).^2));


    u_ref = run_fdtd_2D(u0, v0, Lx, Ly, T, c, dh, dt, gamma, nu);
    
    % Get dimensions from input matrices
    [num_rows, num_cols] = size(THETA1);
    res_surface = zeros(num_rows, num_cols);
    
    total_iterations = num_rows * num_cols;
    current_iteration = 0;
    
    Nx = round(Lx / dh) + 1;
    Ny = round(Ly / dh) + 1;
    Nt = floor(T / dt);

    u_init = rand(Nx,Ny,Nt);
    % u_init = zeros(Nx,Ny,Nt);
    
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
            [~, final_res, ~] = run_swr_2D(u0, v0, N, a, M, Ly, T, c, dh, dt, gamma, nu, theta1, theta2, k, u_init, u_ref);
            
            % Store final residual
            res_surface(jj, ii) = final_res;
            
            % fprintf('  Final residual after %d iterations: %.3e\n', k, final_res);
        end
    end
end