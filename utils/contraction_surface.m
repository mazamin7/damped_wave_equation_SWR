function [Z_p2, Z_inf] = contraction_surface(N, a, M, Ly, y_mode, T, c, dh, dt, gamma, nu, J, THETA1, THETA2)
    % CONTRACTION_SURFACE - Compute contraction factors over theta grid
    % Input: explicit parameters + THETA grids
    % Output: contraction factor surfaces for p=2 and p=inf norms
    
    % Calculate derived parameters
    b = a + M;  % b = a + M
    Lx = N*a + M;

	if Ly ~= 0
		ky = y_mode*pi/Ly;
	else
		ky = 0;
	end

    % T = 0;
    
    % Compute surfaces for both p=2 and p=inf
    Z_p2 = zeros(size(THETA1));
    Z_inf = zeros(size(THETA1));
    
    % Display progress
    total_iterations = numel(THETA1);
    current_iteration = 0;
    
    for i = 1:size(THETA1, 1)
        for j = 1:size(THETA1, 2)
            current_iteration = current_iteration + 1;
            
            % Display progress
            completion_percent = (current_iteration / total_iterations) * 100;
            if mod(current_iteration, 100) == 1 || current_iteration == total_iterations
                fprintf('Progress: %.1f%% (%d/%d) - theta1=%.3f, theta2=%.3f\n', ...
                        completion_percent, current_iteration, total_iterations, ...
                        THETA1(i, j), THETA2(i, j));
            end
            
            % Compute values for both objectives
            Z_p2(i, j) = obj_L2(N, T, dt, J, c, gamma, nu, a, M, THETA1(i, j), THETA2(i, j), ky);
            Z_inf(i, j) = obj_Linf(N, T, dt, J, c, gamma, nu, a, M, THETA1(i, j), THETA2(i, j), ky);
        end
    end
end