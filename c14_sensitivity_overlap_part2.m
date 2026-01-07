%% CLEAN UP
clear; close all; clc;
addpath("utils\")

%% 1. CONFIGURATION
% We sweep M finely to capture the curve shape
M_values = [0, 0.001, 0.002, 0.003, 0.004, 0.005]*10; 
num_M = length(M_values);

% Physics (Fixed in valid regime to isolate M effects)
nu    = 1e-4; 
w_max = 500;   
w_min = 1.0;
c     = 1.0; 
gamma = 0;
N = 2; a = 0.5; J = 1000;

dt = pi / w_max;
T  = pi / w_min;

% Storage
results_Kq   = zeros(num_M, 1);
results_Kp   = zeros(num_M, 1);
results_Krho = zeros(num_M, 1);
results_rho  = zeros(num_M, 1);

options = optimset('Display','off', 'TolX',1e-10, 'TolFun',1e-10);

fprintf('Analyzing sensitivity to Overlap M...\n');

%% 2. SWEEP M
for k = 1:num_M
    current_M = M_values(k);
    
    % Optimize for this specific overlap
    objfun = @(x) obj_Linf(N, T, dt, J, c, gamma, nu, a, current_M, x(1), x(2));
    
    % Initial guesses (Base Theory)
    q0 = 0.5 * nu * w_min * w_max / c^3;
    p0 = 1/c - 0.25 * nu^2 * w_max^2 / c^5;
    
    [x_opt, rho_val] = fminsearch(objfun, [p0, q0], options);
    p_opt = x_opt(1);
    q_opt = x_opt(2);
    
    % --- Extract "Constants" based on Base Theory ---
    % If constants drift, we will see it here.
    
    % K_q = q_opt / (nu * w_min * w_max / c^3)
    results_Kq(k) = q_opt / ((nu * w_min * w_max) / c^3);
    
    % K_p = (1/c - p_opt) / (nu^2 * w_max^2 / c^5)
    drop_p = (1/c) - p_opt;
    results_Kp(k) = drop_p / ((nu^2 * w_max^2) / c^5);
    
    % K_rho = rho_opt / (nu * (w_max - w_min) / c^2)
    % This is the most important one!
    results_Krho(k) = rho_val / ((nu * (w_max - w_min)) / c^2);
    
    results_rho(k) = rho_val;
    
    fprintf('  M=%2d | K_q=%.3f | K_p=%.3f | K_rho=%.3f\n', ...
            current_M, results_Kq(k), results_Kp(k), results_Krho(k));
end

%% 3. CURVE FITTING (Discovery Mode)

% --- Fit K_rho vs M ---
% Hypothesis: Exponential improvement due to spatial decay
% Model: K(M) = K_base * exp(-alpha * M)
fit_rho = fit(M_values', results_Krho, 'exp1'); 
fprintf('\n=== OVERLAP SCALING DISCOVERY ===\n');
fprintf('Convergence Improvement Model: K_rho(M) = A * exp(b * M)\n');
fprintf('  Base K (at M=0): %.4f (Expect ~0.25)\n', fit_rho.a);
fprintf('  Decay Rate (b):  %.4f \n', fit_rho.b);


%% 4. VISUALIZATION
figure('Name','M-Dependence Analysis','Position',[100 100 1200 500], 'Color','w');

% --- Panel 1: Convergence Improvement ---
subplot(1, 3, 1); hold on; grid on; box on;
plot(M_values, results_Krho, 'bo', 'MarkerFaceColor', 'b', 'MarkerSize', 8);
plot(fit_rho, 'b-'); % Plot the exponential fit
yline(0.25, 'k--', 'Zero Overlap Theory (0.25)');
xlabel('Overlap $M$', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('Scaling Factor $K_{\rho}$', 'Interpreter', 'latex', 'FontSize', 14);
title('Convergence vs Overlap', 'FontSize', 16);
legend('Sim Data', 'Exp Fit', 'Location', 'best');

% --- Panel 2: Damping Parameter q ---
subplot(1, 3, 2); hold on; grid on; box on;
plot(M_values, results_Kq, 'rs-', 'LineWidth', 2, 'MarkerFaceColor', 'r');
yline(0.50, 'k--', 'Base Theory (0.50)');
xlabel('Overlap $M$', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('Scaling Factor $K_q$', 'Interpreter', 'latex', 'FontSize', 14);
title('Stability of $q$', 'FontSize', 16);
ylim([0, 0.6]);

% --- Panel 3: Transport Parameter p ---
subplot(1, 3, 3); hold on; grid on; box on;
plot(M_values, results_Kp, 'gs-', 'LineWidth', 2, 'MarkerFaceColor', 'g');
yline(0.25, 'k--', 'Base Theory (0.25)');
xlabel('Overlap $M$', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('Scaling Factor $K_p$', 'Interpreter', 'latex', 'FontSize', 14);
title('Stability of $p$', 'FontSize', 16);
ylim([0, 0.4]);