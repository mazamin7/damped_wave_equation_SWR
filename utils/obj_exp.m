function val = obj_exp(x, u0, v0, N, a, M, T, c, dh, dt, gamma, nu, k_optim, u_init, u_ref)
    % obj_exp: Experimental objective function with Amplification Factor scaling
    
    p = x(1);
    q = x(2);

    fprintf(1, '   [Eval] p = %10.6f,  q = %10.6f ... \n', p, q);

    % Basic constraint check
    if p < 0 
       val = 1e10; 
       fprintf(1, 'Skipped (p<0)\n'); 
       return;
    end

    try
        % --- STEP 1: Compute Amplification Factor (F) ---
        
        % Calculate initial error of the raw random guess
        err0 = max(abs(u_init(:) - u_ref(:))) / max(abs(u_ref(:)));
        
        % Run 1 iteration to get E1
        [~, ~, res_test] = run_swr_1D(u0, v0, N, a, M, T, c, dh, dt, ...
                                      gamma, nu, p, q, 1, u_init, u_ref);
        E1 = res_test(1);
        
        % Calculate F
        F = E1 / err0;
        
        % Safety check to avoid division by zero
        if F == 0, F = 1; end
        
        % --- STEP 2: Run Real Simulation with Scaled Input ---
        
        u_init_scaled = u_init / F;
        
        [~, ~, res_hist] = run_swr_1D(u0, v0, N, a, M, T, c, dh, dt, ...
                                      gamma, nu, p, q, k_optim, u_init_scaled, u_ref);
        
        % The objective value is the error at the last iteration
        val = res_hist(end);
        
        % Check for divergence
        if isnan(val) || isinf(val)
            val = 1e20;
            fprintf(1, 'Failed (NaN/Inf)\n');
        else
            fprintf(1, 'Err = %.4e\n', val);
        end
        
    catch
        % If simulation crashes
        val = 1e20;
        fprintf(1, 'Crashed\n');
    end
end