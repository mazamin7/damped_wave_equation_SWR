function val = obj_Linf(N, T, dt, J, c, gamma, nu, a, M, theta1, theta2)
    % L∞ objective function with explicit parameters
    b = a + M;
    
    % Frequency grid
    omega_min = pi / T;
    omega_max = pi / dt;
    % omega_min = 0;
    % omega_max = 20 * pi / dt;
    % omegas = linspace(omega_min, omega_max, J);
    omegas = logspace(log(omega_min), log(omega_max), J);
    s_vals = 1i * omegas;
    
    % Compute contraction factors
    r_vals = zeros(size(s_vals));
    for idx = 1:length(s_vals)
        r_vals(idx) = rho(N, s_vals(idx), theta1, theta2, c, gamma, nu, a, b);
    end
    
    % L∞ norm (maximum over frequency band)
    val = max(abs(r_vals));
end