%% CLEAN UP
clear; close all; clc;
addpath("utils\")

%% 1. CONFIGURATION
c  = 1.0;
nu = 0; % Telegrapher
N  = 2; a = 0.5; M = 0; J = 1000;

% We generate a 2D grid of (gamma, w_min)
% specifically targeting the "Transition Regime" where terms are comparable
fprintf('Generating Data to fit A and B...\n');

num_g = 20;
num_w = 20;
gamma_vals = logspace(-4.0, -2, num_g); % 0.1 to 3.0
wmin_vals  = logspace(-1.0, 0.5, num_w); % 0.1 to 3.0

w_max = 100;
dt = pi / w_max;

% Storage
X_gamma = [];
X_wmin  = [];
Y_rho   = [];

options = optimset('Display','off', 'TolX',1e-12, 'TolFun',1e-12);

%% 2. GENERATE GROUND TRUTH
counter = 0;
for i = 1:num_g
    g_val = gamma_vals(i);
    for k = 1:num_w
        w_val = wmin_vals(k);
        T_val = pi / w_val;
        
        % Run Optimization
        objfun = @(x) obj_Linf(N, T_val, dt, J, c, g_val, nu, a, M, x(1), x(2));
        
        % Initial Guess (Telegrapher Optimized)
        p0 = 1/c + g_val^2 / (16 * c * w_val^2);
        q0 = g_val / (2*c);
        
        [~, rho_val] = fminsearch(objfun, [p0, q0], options);
        
        X_gamma = [X_gamma; g_val];
        X_wmin  = [X_wmin;  w_val];
        Y_rho   = [Y_rho;   rho_val];
        
        counter = counter + 1;
    end
    if mod(i, 5) == 0, fprintf('Generated %d points...\n', counter); end
end

%% 3. FIT 2 PARAMETERS (A, B)
% Fixed denominator constant C = 64
C_fixed = 64.0;

% Model: rho = gamma^2 * sqrt(A*g^2 + B*w^2) / (64 * w^3)
func_2param = @(p, x) (x(:,1).^2) .* ...
    sqrt(p(1)*x(:,1).^2 + p(2)*x(:,2).^2) ./ (C_fixed * x(:,2).^3);

X_data = [X_gamma, X_wmin];

% Initial guess [A, B] = [9, 4]
p_init = [9, 4];
lb = [0, 0];
ub = [100, 100];

[p_fit, resnorm, ~, ~, ~] = lsqcurvefit(func_2param, p_init, X_data, Y_rho, lb, ub, optimset('Display','off'));

%% 4. ANALYZE RESULTS
A_fit = p_fit(1);
B_fit = p_fit(2);

Y_pred = func_2param(p_fit, X_data);
mean_err = mean(abs(Y_pred - Y_rho) ./ Y_rho * 100);

fprintf('\n=== 2-PARAMETER FIT RESULTS ===\n');
fprintf('Model: rho = gamma^2 * sqrt(A*gamma^2 + B*wmin^2) / (64 * wmin^3)\n');
fprintf('--------------------------------------------------\n');
fprintf('Parameter | Expected | Fitted   | Deviation\n');
fprintf('--------------------------------------------------\n');
fprintf('    A     | 9.00     | %.4f   | %.2f%%\n', A_fit, abs(A_fit-9)/9*100);
fprintf('    B     | 4.00     | %.4f   | %.2f%%\n', B_fit, abs(B_fit-4)/4*100);
fprintf('--------------------------------------------------\n');
fprintf('Mean Relative Error of Model: %.4f%%\n', mean_err);

%% 5. VISUALIZATION
figure('Name','2-Param Coefficient Fit','Position',[100 100 1000 500], 'Color','w');

% Plot 1: Correlation
subplot(1, 2, 1); hold on; grid on; box on;
plot(Y_rho, Y_pred, 'bo', 'MarkerFaceColor', 'b', 'MarkerSize', 4);
plot([min(Y_rho) max(Y_rho)], [min(Y_rho) max(Y_rho)], 'k--', 'LineWidth', 2);
xlabel('Actual Simulation $\rho$', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('Predicted $\rho$ (Formula)', 'Interpreter', 'latex', 'FontSize', 14);
title('Model Accuracy', 'FontSize', 16);
axis square;

% Plot 2: Residuals vs Regime (Gamma/w_min)
subplot(1, 2, 2); hold on; grid on; box on;
ratio = X_gamma ./ X_wmin;
residuals = (Y_pred - Y_rho) ./ Y_rho * 100;
semilogx(ratio, residuals, 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 4);
yline(0, 'k-');
xlabel('Regime Ratio $\gamma / \omega_{\min}$', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('Relative Error (%)', 'Interpreter', 'latex', 'FontSize', 14);
title('Residuals across Regimes', 'FontSize', 16);