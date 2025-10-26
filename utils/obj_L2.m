% Use trapezoidal rule to approximate integral and normalize by bandwidth
function val = obj_L2(N, T, dt, J, c, gamma, nu, a, M, theta1, theta2, ky)
    % L2 objective function with explicit parameters
    b = a + M;
    
    % Frequency grid
    omega_min = 2*pi / T;
    omega_max = pi / dt;
    omegas = linspace(omega_min, omega_max, J);
    s_vals = 1i * omegas;
    
    % Compute contraction factors
    r_vals = zeros(size(s_vals));
    for idx = 1:length(s_vals)
        r_vals(idx) = rho(N, s_vals(idx), theta1, theta2, c, gamma, nu, a, b, ky);
    end
    
    % L2 norm (RMS over frequency band)
    bandwidth = (omega_max - omega_min);
    val = sqrt(trapz(omegas, abs(r_vals).^2) / bandwidth);
end