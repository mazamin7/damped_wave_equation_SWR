%% CLEAN UP
clear; close all; clc;
addpath("utils\")

%% 1. CONFIGURATION
% Physics
c = 1.0;
nu = 0; % Pure Telegrapher
N = 2; a = 0.5; M = 0; J = 1000;

% We vary Gamma to generate the curve
NumPoints = 50;
gamma_vals = logspace(-4, -2, NumPoints); % 0.01 to 1.0

% Fixed Frequency Limits (to isolate the coefficient)
w_min = 1.0; 
w_max = 100.0;
dt    = pi / w_max;
T     = pi / w_min;

% Storage
data_p = zeros(NumPoints, 1);
data_q = zeros(NumPoints, 1);

options = optimset('Display','off', 'TolX',1e-12, 'TolFun',1e-12);

fprintf('Generating Optimization Data (Gamma sweep)...\n');

%% 2. GENERATE DATA
for i = 1:NumPoints
    g_val = gamma_vals(i);
    
    objfun = @(x) obj_Linf(N, T, dt, J, c, g_val, nu, a, M, x(1), x(2));
    
    % Initial Guess
    p0 = 1/c + g_val^2 / (16 * c * w_min^2); % Your hypothesis
    q0 = g_val / (2*c);
    
    [x_opt, ~] = fminsearch(objfun, [p0, q0], options);
    data_p(i) = x_opt(1);
    data_q(i) = x_opt(2);
    
    if mod(i, 10) == 0, fprintf('.'); end
end
fprintf('\nData generation complete.\n');

%% 3. FITTING THE COEFFICIENTS

% --- Fit q ---
% Model: q = K_q * (gamma / c)
X_q = gamma_vals(:) ./ c;
Y_q = data_q(:);
K_q_fit = X_q \ Y_q; % Linear regression (slope)

% --- Fit p ---
% Model: p = 1/c + K_p * [ gamma^2 / (c * w_min^2) ]
% We subtract the base 1/c to isolate the correction
X_p = (gamma_vals(:).^2) ./ (c * w_min^2);
Y_p = data_p(:) - (1/c);
K_p_fit = X_p \ Y_p; % Linear regression (slope)

%% 4. RESULTS
fprintf('\n=== FITTING RESULTS ===\n');

fprintf('\n1. Damping Parameter (q)\n');
fprintf('   Model: q = K * (gamma / c)\n');
fprintf('   Fitted K:      %.4f\n', K_q_fit);
fprintf('   Theoretical:   0.5000 (1/2)\n');
fprintf('   Error:         %.2f%%\n', abs(K_q_fit - 0.5)*200);

fprintf('\n2. Transport Parameter (p)\n');
fprintf('   Model: p = 1/c + K * [ gamma^2 / (c * w_min^2) ]\n');
fprintf('   Fitted K:      %.4f\n', K_p_fit);
fprintf('   Strict Taylor: 0.1250 (1/8)\n');
fprintf('   Minimax Hyp:   0.0625 (1/16 = 1/2 * 1/8)\n');
fprintf('   Error vs 1/16: %.2f%%\n', abs(K_p_fit - 0.0625)/0.0625 * 100);


%% 5. VISUALIZATION
figure('Name','Telegrapher Coefficients Fit','Position',[100 100 1000 500], 'Color','w');

% --- Q Plot ---
subplot(1, 2, 1); hold on; grid on; box on;
plot(X_q, Y_q, 'bo', 'MarkerFaceColor','b', 'DisplayName','Optimal q');
plot(X_q, K_q_fit * X_q, 'r-', 'LineWidth', 2, 'DisplayName', sprintf('Fit (K=%.4f)', K_q_fit));
xlabel('Scaling Term $\gamma / c$', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('Optimal $q$', 'Interpreter', 'latex', 'FontSize', 14);
title('Damping Parameter Fit', 'FontSize', 16);
legend('Location','best');

% --- P Plot ---
subplot(1, 2, 2); hold on; grid on; box on;
plot(X_p, Y_p, 'ro', 'MarkerFaceColor','r', 'DisplayName','Optimal p correction');
plot(X_p, K_p_fit * X_p, 'b-', 'LineWidth', 2, 'DisplayName', sprintf('Fit (K=%.4f)', K_p_fit));
plot(X_p, 0.125 * X_p, 'k--', 'LineWidth', 1, 'DisplayName', 'Strict Taylor (1/8)');
plot(X_p, 0.0625 * X_p, 'g--', 'LineWidth', 2, 'DisplayName', 'Minimax (1/16)');

xlabel('Scaling Term $\gamma^2 / (c \omega_{\min}^2)$', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('Correction $(p - 1/c)$', 'Interpreter', 'latex', 'FontSize', 14);
title('Transport Parameter Fit', 'FontSize', 16);
legend('Location','best');