%% CLEAN UP
clear; close all; clc;
addpath("utils\")

%% 1. CONFIGURATION
% Physics: Telegrapher Regime
gamma = 0.5;   % Fixed Gamma
c     = 1.0;
nu    = 0;
N     = 2; a = 0.5; M = 0; J = 1000;

% Sweep w_min (Simulation Time T)
% We go from very long simulations (small w_min) to short ones
NumPoints = 40;
w_min_vals = logspace(-1, 0.5, NumPoints); % 0.1 to ~3.14 rad/s

% Fixed High Frequency
w_max = 100; 
dt    = pi / w_max;

% Storage
data_p = zeros(NumPoints, 1);
data_q = zeros(NumPoints, 1);

options = optimset('Display','off', 'TolX',1e-12, 'TolFun',1e-12);

fprintf('Sweeping omega_min (T) to check scaling laws...\n');

%% 2. GENERATE DATA
for i = 1:NumPoints
    w_min = w_min_vals(i);
    T_val = pi / w_min;
    
    objfun = @(x) obj_Linf(N, T_val, dt, J, c, gamma, nu, a, M, x(1), x(2));
    
    % Initial Guess
    p0 = 1/c + 0.0625 * gamma^2 / (c * w_min^2); 
    q0 = gamma / (2*c);
    
    [x_opt, ~] = fminsearch(objfun, [p0, q0], options);
    data_p(i) = x_opt(1);
    data_q(i) = x_opt(2);
    
    if mod(i, 10) == 0, fprintf('.'); end
end
fprintf('\n');

%% 3. FITTING POWER LAWS

% --- Fit p vs w_min ---
% Model: (p - 1/c) = A * w_min^B
% Log-Log Linear Regression: log(y) = log(A) + B * log(x)
Y_p = data_p - (1/c);
% Filter out numerical noise (ensure Y_p > 0)
valid_idx = Y_p > 1e-10;
log_x = log(w_min_vals(valid_idx)');
log_y = log(Y_p(valid_idx));

coeffs_p = polyfit(log_x, log_y, 1);
power_p = coeffs_p(1); % Slope B
const_p = exp(coeffs_p(2)); % Intercept A

% Extract the coefficient K from A = K * gamma^2 / c
% K = A * c / gamma^2
K_p_extracted = const_p * c / gamma^2;


% --- Fit q vs w_min ---
% Check standard deviation to see if it's constant
q_mean = mean(data_q);
q_variation = std(data_q) / q_mean * 100;


%% 4. RESULTS
fprintf('\n=== SCALING LAW VERIFICATION ===\n');

fprintf('\n1. Transport Parameter (p)\n');
fprintf('   Model: (p - 1/c) ~ 1 / w_min^B\n');
fprintf('   Fitted Power B: %.4f (Expect -2.00)\n', power_p);
fprintf('   Extracted K:    %.4f (Expect ~0.0625)\n', K_p_extracted);

fprintf('\n2. Damping Parameter (q)\n');
fprintf('   Mean Value:     %.4f (Theory: %.4f)\n', q_mean, gamma/(2*c));
fprintf('   Variation:      %.2f%%\n', q_variation);
if q_variation < 1.0
    fprintf('   => Conclusion: q is CONSTANT wrt w_min.\n');
else
    fprintf('   => Conclusion: q has a dependency.\n');
end

%% 5. VISUALIZATION
figure('Name','Frequency Scaling Check','Position',[100 100 1000 500], 'Color','w');

% --- P Scaling Plot ---
subplot(1, 2, 1); grid on; box on;
loglog(w_min_vals, Y_p, 'bo', 'MarkerFaceColor','b', 'DisplayName','Simulation Data');
hold on;
% Plot the theoretical line with slope -2
y_theory = (0.0625 * gamma^2 / c) .* (w_min_vals.^(-2));
loglog(w_min_vals, y_theory, 'r-', 'LineWidth', 2, 'DisplayName','Theory (Slope -2, K=1/16)');

xlabel('$\omega_{\min} = \pi/T$', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('$(p - 1/c)$', 'Interpreter', 'latex', 'FontSize', 14);
title('Scaling of $p$ correction', 'FontSize', 16);
legend('Location','best');

% --- Q Stability Plot ---
subplot(1, 2, 2); grid on; box on;
semilogx(w_min_vals, data_q, 'rs-', 'MarkerFaceColor','r', 'LineWidth', 1.5);
hold on;
yline(gamma/(2*c), 'k--', 'Theory \gamma/2c', 'LineWidth', 2);

xlabel('$\omega_{\min} = \pi/T$', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('Optimal $q$', 'Interpreter', 'latex', 'FontSize', 14);
title('Stability of $q$', 'FontSize', 16);
ylim([0, max(data_q)*1.5]);