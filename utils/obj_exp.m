function val = obj_exp(x, u0, v0, N, a, M, T, c, dh, dt, gamma, nu, k_optim, u_init, u_ref)
    % obj_exp: Experimental objective function for SWR optimization
    %
    % Inputs:
    %   x        : [p, q] optimization variables
    %   [...params] : Physical and numerical parameters
    %   k_optim  : Number of iterations to run for the optimization step
    %              (keep this low, e.g., 5-10, for speed)
    %   u_init   : Initial random guess for the Schwarz iteration
    %   u_ref    : The reference FDTD solution (to compute error)
    
    p = x(1);
    q = x(2);

    fprintf(1, '   [Eval] p = %10.6f,  q = %10.6f ... \n', p, q);

    try
        % Run the experimental simulation
        % We use evalc to silence the print outputs from run_swr_1D during optimization
        % If your run_swr_1D is already silent, you can remove evalc.
        [~, ~, res_hist] = run_swr_1D(u0, v0, N, a, M, T, c, dh, dt, ...
                                      gamma, nu, p, q, k_optim, u_init, u_ref);
        
        % The objective value is the error at the last iteration
        val = res_hist(end);
        
        % Handling divergence: if error is NaN or Inf, return a large penalty
        if isnan(val) || isinf(val)
            val = 1e20;
        end
        
    catch
        % If the simulation crashes (e.g., instability), return large penalty
        val = 1e20;
    end
end