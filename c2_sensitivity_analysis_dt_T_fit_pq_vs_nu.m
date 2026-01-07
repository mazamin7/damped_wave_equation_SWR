%% CLEAN UP
clear; close all; clc;
addpath("utils\")

%% 1. GLOBAL CONFIGURATION
% Base parameters from get_sim_params
P_base = get_sim_params();
nu_range = logspace(-7, -4, 20); % Range for fitting (small nu)

% Arrays of parameters to sweep
% NOTE: Changing dt changes omega_max (approx pi/dt)
%       Changing T  changes omega_min (approx pi/T)
dt_sweep = [0.001, 0.002, 0.004, 0.005, 0.008, 0.01]; 
T_sweep  = [2, 5, 10, 20, 50];

%% 2. SWEEP 1: VARYING DT (Omega_max)
fprintf('=== SWEEP 1: Varying dt (Frequency Upper Bound) ===\n');
Kq_dt = zeros(size(dt_sweep));
Kp_dt = zeros(size(dt_sweep));

% Fix T for this sweep
T_fixed = 5.0; 

for i = 1:length(dt_sweep)
    dt_val = dt_sweep(i);
    fprintf('  Running for dt = %.4f...\n', dt_val);
    
    % Run optimization wrapper
    [kq, kp] = get_fitted_coefficients(P_base, T_fixed, dt_val, nu_range);
    
    Kq_dt(i) = kq;
    Kp_dt(i) = kp;
end

%% 3. SWEEP 2: VARYING T (Omega_min)
fprintf('=== SWEEP 2: Varying T (Frequency Lower Bound) ===\n');
Kq_T = zeros(size(T_sweep));
Kp_T = zeros(size(T_sweep));

% Fix dt for this sweep
dt_fixed = 0.002;

for i = 1:length(T_sweep)
    T_val = T_sweep(i);
    fprintf('  Running for T = %.1f...\n', T_val);
    
    [kq, kp] = get_fitted_coefficients(P_base, T_val, dt_fixed, nu_range);
    
    Kq_T(i) = kq;
    Kp_T(i) = kp;
end

%% 4. VISUALIZATION OF DEPENDENCIES
figure('Name','Parameter Dependencies','Position',[100 100 1000 700], 'Color','w');

% --- ROW 1: DEPENDENCE ON DT ---
subplot(2, 2, 1);
loglog(dt_sweep, abs(Kq_dt), 'b-o', 'LineWidth', 2, 'MarkerFaceColor','b');
xlabel('\Delta t'); ylabel('Slope K_q');
title('Sensitivity of q to \Delta t (\omega_{max})');
grid on;

subplot(2, 2, 2);
loglog(dt_sweep, abs(Kp_dt), 'r-s', 'LineWidth', 2, 'MarkerFaceColor','r');
xlabel('\Delta t'); ylabel('Curvature K_p');
title('Sensitivity of p to \Delta t (\omega_{max})');
grid on;

% --- ROW 2: DEPENDENCE ON T ---
subplot(2, 2, 3);
loglog(T_sweep, abs(Kq_T), 'b-o', 'LineWidth', 2, 'MarkerFaceColor','b');
xlabel('T (Simulation Time)'); ylabel('Slope K_q');
title('Sensitivity of q to T (\omega_{min})');
grid on;

subplot(2, 2, 4);
loglog(T_sweep, abs(Kp_T), 'r-s', 'LineWidth', 2, 'MarkerFaceColor','r');
xlabel('T (Simulation Time)'); ylabel('Curvature K_p');
title('Sensitivity of p to T (\omega_{min})');
grid on;

%% 5. ANALYSIS OUTPUT
fprintf('\n=== SCALING ANALYSIS ===\n');
% Check power law for dt: K ~ dt^alpha
p_poly_q = polyfit(log(dt_sweep), log(abs(Kq_dt)), 1);
p_poly_p = polyfit(log(dt_sweep), log(abs(Kp_dt)), 1);

fprintf('Scaling with dt (Expect -2 for omega^2 dependence):\n');
fprintf('  K_q ~ dt^(%.2f)\n', p_poly_q(1));
fprintf('  K_p ~ dt^(%.2f)\n', p_poly_p(1));

% Check power law for T: K ~ T^alpha
p_poly_q = polyfit(log(T_sweep), log(abs(Kq_T)), 1);
p_poly_p = polyfit(log(T_sweep), log(abs(Kp_T)), 1);

fprintf('Scaling with T (Expect -2 for omega^2 dependence):\n');
fprintf('  K_q ~ T^(%.2f)\n', p_poly_q(1));
fprintf('  K_p ~ T^(%.2f)\n', p_poly_p(1));


%% ====================================================
%  HELPER FUNCTION: RUN OPTIMIZATION & FIT CURVES
%  ====================================================
function [K_q, K_p] = get_fitted_coefficients(P, T_val, dt_val, nu_vals)
    % Unpack fixed parameters
    N = P.N; c = P.c; a = P.a; M = P.M; J = P.J;
    gamma = 0;
    
    num_nu = length(nu_vals);
    q_real = zeros(1, num_nu);
    p_real = zeros(1, num_nu);
    
    % Optimization Setup
    options = optimset('Display','off', 'TolX',1e-8, 'TolFun',1e-8);
    x0 = [1/c, 0]; 
    
    for k = 1:num_nu
        nu = nu_vals(k);
        % Construct objective with current T and dt
        objfun = @(x) obj_Linf(N, T_val, dt_val, J, c, gamma, nu, a, M, x(1), x(2));
        
        [x_opt, ~] = fminsearch(objfun, x0, options);
        p_real(k) = x_opt(1);
        q_real(k) = x_opt(2);
        x0 = x_opt; % Warm start
    end
    
    % --- Fit q(nu) = K_q * nu + C ---
    % Using robust linear fit (polyfit)
    coeffs_q = polyfit(nu_vals, q_real, 1);
    K_q = coeffs_q(1); % Slope
    
    % --- Fit p(nu) = K_p * nu + C ---
    % Fit p against nu
    coeffs_p = polyfit(nu_vals, p_real, 1);
    K_p = coeffs_p(1); % Curvature (Note: usually negative)
end