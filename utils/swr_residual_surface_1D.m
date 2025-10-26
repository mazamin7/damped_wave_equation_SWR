function res_surface = swr_residual_surface_1D(N, a, M, T, c, dh, dt, gamma, nu, k, THETA1, THETA2)
    % SWR_RESIDUAL_SURFACE - Parameter sweep over theta1, theta2
    % Input: explicit parameters + THETA1, THETA2 grids
    % Output: final residuals for each parameter combination
    
    % Compute reference solution once
    fprintf('Computing reference FDTD solution...\n');
    Lx = N*a + M;  % Calculate domain length


    % Initial conditions components
    gaussian = @(r,mu,sigma) 1/(2*pi*sigma^2) * exp(-(r-mu).^2/(2*sigma^2));
    
    % Define n_max - you can adjust this value as needed
    n_max = 100;
    
    % Compute max values for normalization
    x_test = linspace(0,Lx,1000);
    gaussian_max = max(gaussian(x_test, Lx/4, Lx/20));
    
    % Compute max of sine sum
    sine_sum_vals = zeros(size(x_test));
    for n = 1:n_max
        sine_sum_vals = sine_sum_vals + sin(n*pi*x_test/Lx);
    end
    sine_sum_max = max(abs(sine_sum_vals));
    
    % Create normalized components as function handles
    gaussian_normalized = @(x) gaussian(x, Lx/4, Lx/20) / gaussian_max;
    
    % Create normalized sine sum function handle
    sine_sum_normalized = @(x) arrayfun(@(xi) sum(sin((1:n_max)' * pi * xi / Lx)) / sine_sum_max, x);
    
    % Create the final initial condition function
    u0 = @(x) gaussian_normalized(x) + sine_sum_normalized(x);
    v0 = @(x) 0;


    dt = dh/c;     % Calculate time step
    u_ref = run_fdtd_1D(u0, v0, Lx, T, c, dh, dt, gamma, nu);
    
    % Get dimensions from input matrices
    [num_rows, num_cols] = size(THETA1);
    res_surface = zeros(num_rows, num_cols);
    
    total_iterations = num_rows * num_cols;
    current_iteration = 0;
    
    Nx = round(Lx / dh) + 1;
    Nt = floor(T / dt);

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
            [~, final_res, ~] = run_swr_1D(u0, v0, N, a, M, T, c, dh, dt, gamma, nu, theta1, theta2, k, u_init, u_ref);
            
            % Store final residual
            res_surface(jj, ii) = final_res;
            
            % fprintf('  Final residual after %d iterations: %.3e\n', k, final_res);
        end
    end
end