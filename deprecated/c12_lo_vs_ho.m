%% CLEAN UP
clear; close all; clc;
addpath("utils\")

%% 1. CONFIGURATION
gamma = 0;
c     = 1.0;   
N = 2; a = 0.5; M = 0; J = 1000;

%% 2. GENERATE VALIDATION DATASET
NumPoints = 200; 
fprintf('Generating %d validation points...\n', NumPoints);

% Ranges designed to stress-test the approximation
range_nu    = [1e-5, 1e-3];   
range_w_max = [15, 40]; 
range_w_min = [0.1, 1.0];

% Storage
data_p_opt = zeros(NumPoints, 1);
data_nu    = zeros(NumPoints, 1);
data_wmax  = zeros(NumPoints, 1);

options = optimset('Display','off', 'TolX',1e-12, 'TolFun',1e-12);

for i = 1:NumPoints
    % 1. Random Parameters
    nu_val = range_nu(1) + rand()*diff(range_nu);
    w_max_val = range_w_max(1) + rand()*diff(range_w_max);
    w_min_val = range_w_min(1) + rand()*diff(range_w_min);
    
    dt_val = pi / w_max_val;
    T_val  = pi / w_min_val;
    
    data_nu(i)   = nu_val;
    data_wmax(i) = w_max_val;
    
    % 2. Ground Truth (Numerical Optimization)
    objfun = @(x) obj_Linf(N, T_val, dt_val, J, c, gamma, nu_val, a, M, x(1), x(2));
    
    % Good guess speeds up convergence
    guess_p = 1/c - 0.25 * nu_val^2 * w_max_val^2/c^5;
    guess_q = 0.5 * nu_val * w_min_val * w_max_val/c^3;
    
    [x_opt, ~] = fminsearch(objfun, [guess_p, guess_q], options);
    data_p_opt(i) = x_opt(1);
    
    if mod(i, 10) == 0, fprintf('.'); end
end
fprintf('\n');

%% 3. DEFINE MODELS & FIT COEFFICIENTS

% The variable we are modeling is the "Drop" from 1/c
% Y = (1/c - p_opt)
Y_obs = (1/c) - data_p_opt;

% Feature 1: Leading Order Term (nu^2 * w^2) / c^5
X_LO = (data_nu.^2 .* data_wmax.^2) / c^5;

% Feature 2: Higher Order Term (nu^4 * w^4) / c^9
X_HO = (data_nu.^4 .* data_wmax.^4) / c^9;

% --- MODEL 1: Standard (Low Order) ---
% Hypothesis: Y = K1 * X_LO
% We fit K1 to give the LO model the best possible chance
K_LO_best = X_LO \ Y_obs; 
Y_pred_LO = K_LO_best * X_LO;

% --- MODEL 2: Refined (High Order) ---
% Hypothesis: Y = K1 * X_LO - K2 * X_HO
% Note: We subtract K2 because the correction adds to p, reducing the drop.
% We perform a multivariate regression to find optimal K1 and K2 together.
DesignMatrix = [X_LO, -X_HO];
Coeffs_HO    = DesignMatrix \ Y_obs;

K_HO_1 = Coeffs_HO(1); % Coefficient for nu^2 term
K_HO_2 = Coeffs_HO(2); % Coefficient for nu^4 term

Y_pred_HO = DesignMatrix * Coeffs_HO;

%% 4. COMPUTE ACCURACY METRICS

% Convert back to p values for physical error checking
p_pred_LO = (1/c) - Y_pred_LO;
p_pred_HO = (1/c) - Y_pred_HO;

% Relative Error (%)
err_LO = abs(p_pred_LO - data_p_opt) ./ abs(data_p_opt) * 100;
err_HO = abs(p_pred_HO - data_p_opt) ./ abs(data_p_opt) * 100;

fprintf('\n=== MODEL COMPARISON ===\n');
fprintf('Model 1: Standard (Low Order)\n');
fprintf('  Formula: p ~ 1/c - %.4f * (nu*w)^2\n', K_LO_best);
fprintf('  Max Relative Error: %.4f%%\n', max(err_LO));
fprintf('  Mean Relative Error: %.4f%%\n', mean(err_LO));

fprintf('\nModel 2: Refined (High Order)\n');
fprintf('  Formula: p ~ 1/c - %.4f * (nu*w)^2 + %.4f * (nu*w)^4\n', K_HO_1, K_HO_2);
fprintf('  Max Relative Error: %.4f%%\n', max(err_HO));
fprintf('  Mean Relative Error: %.4f%%\n', mean(err_HO));

improvement = mean(err_LO) / mean(err_HO);
fprintf('\n=> Accuracy Improvement Factor: %.2fx\n', improvement);

%% 5. VISUALIZATION
figure('Name','LO vs HO Model Comparison','Position',[100 100 1000 500], 'Color','w');

% --- Panel 1: Residuals vs Perturbation Strength ---
subplot(1, 2, 1); hold on; grid on; box on;
x_axis = X_LO; % Use the LO term as the x-axis "strength" metric

% Plot Errors
plot(x_axis, err_LO, 'bo', 'MarkerFaceColor','b', 'MarkerSize', 6, 'DisplayName', 'Low Order Error');
plot(x_axis, err_HO, 'rs', 'MarkerFaceColor','r', 'MarkerSize', 6, 'DisplayName', 'High Order Error');

xlabel('Perturbation Strength $(\nu \omega_{\max} / c^2)^2$', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('Relative Error in $p$ (%)', 'Interpreter', 'latex', 'FontSize', 14);
title('Error vs Damping Strength', 'FontSize', 16);
legend('Location','northwest');
axis square;

% --- Panel 2: Visualizing the Curve Fit ---
subplot(1, 2, 2); hold on; grid on; box on;

% Sort data for clean lines
[sorted_x, idx] = sort(x_axis);
sorted_y_obs = Y_obs(idx);
sorted_y_LO  = Y_pred_LO(idx);
sorted_y_HO  = Y_pred_HO(idx);

% Plot Observed Drop
plot(sorted_x, sorted_y_obs, 'ko', 'MarkerSize', 6, 'LineWidth', 1.5, 'DisplayName', 'Simulation Data');

% Plot Models
plot(sorted_x, sorted_y_LO, 'b--', 'LineWidth', 2, 'DisplayName', ['LO Model (Linear in \nu^2)']);
plot(sorted_x, sorted_y_HO, 'r-', 'LineWidth', 2, 'DisplayName', ['HO Model (Quadratic in \nu^2)']);

xlabel('Leading Term $(\nu \omega_{\max} / c^2)^2$', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('Parameter Drop $(1/c - p)$', 'Interpreter', 'latex', 'FontSize', 14);
title('Curve Fitting Performance', 'FontSize', 16);
legend('Location','best');
axis square;