%% CLEAN UP
clear; close all; clc;
addpath("utils\")

%% 1. CONFIGURATION
gamma = 0;
c     = 1.0;   
N = 2; a = 0.5; M = 0; J = 1000;

% Training settings
NumTrainPoints = 200;
NumTestPoints  = 200;

% Ranges (Logarithmic sampling)
range_nu    = [1e-5, 1e-1];   
range_w_max = [10, 40]; 
range_w_min = [0.1, 1.0];

% Optimization options
options = optimset('Display','off', 'TolX',1e-12, 'TolFun',1e-12);

%% ========================================================================
%% PHASE 1: TRAINING (Compute Coefficients Dynamically)
%% ========================================================================
fprintf('Phase 1: Generating Training Data (%d points)...\n', NumTrainPoints);

train_nu   = zeros(NumTrainPoints, 1);
train_wmax = zeros(NumTrainPoints, 1);
train_wmin = zeros(NumTrainPoints, 1);
train_p    = zeros(NumTrainPoints, 1);
train_q    = zeros(NumTrainPoints, 1);

for i = 1:NumTrainPoints
    % Logarithmic Random Sampling
    nu_val = 10^(log10(range_nu(1)) + rand()*diff(log10(range_nu)));
    w_max_val = 10^(log10(range_w_max(1)) + rand()*diff(log10(range_w_max)));
    w_min_val = range_w_min(1) + rand()*diff(range_w_min);
    
    dt_val = pi / w_max_val;
    T_val  = pi / w_min_val;
    
    % Store Inputs
    train_nu(i)   = nu_val;
    train_wmax(i) = w_max_val;
    train_wmin(i) = w_min_val;
    
    % Find Ground Truth (Numerical Optimization)
    objfun = @(x) obj_Linf(N, T_val, dt_val, J, c, gamma, nu_val, a, M, x(1), x(2));
    
    % Robust Initial Guess
    p_guess = 1/c - 0.25 * nu_val^2 * w_max_val^2 / c^5;
    q_guess = 0.5 * nu_val * w_min_val * w_max_val / c^3;
    
    [x_opt, ~] = fminsearch(objfun, [p_guess, q_guess], options);
    train_p(i) = x_opt(1);
    train_q(i) = x_opt(2);
    
    if mod(i, 20) == 0, fprintf('.'); end
end
fprintf('\nTraining Data Complete. Fitting Models...\n');

% --- FIT LO MODEL (Single Coefficient) ---
term_p_2 = (train_nu.^2 .* train_wmax.^2) / c^5;
term_q_1 = (train_nu .* train_wmin .* train_wmax) / c^3;
Y_p_drop = (1/c) - train_p; 
Y_q      = train_q;

K_p_LO = term_p_2 \ Y_p_drop;
K_q_LO = term_q_1 \ Y_q;

% --- FIT HO MODEL (Multi-Coefficient) ---
term_p_4 = (train_nu.^4 .* train_wmax.^4) / c^9;
term_q_3 = (train_nu.^3 .* train_wmax.^4) / c^7;

A_p = [term_p_2, -term_p_4];
coeffs_p = A_p \ Y_p_drop;
K_p_HO1 = coeffs_p(1);
K_p_HO2 = coeffs_p(2);

A_q = [term_q_1, -term_q_3];
coeffs_q = A_q \ Y_q;
K_q_HO1 = coeffs_q(1);
K_q_HO2 = coeffs_q(2);

fprintf('\n=== TRAINED COEFFICIENTS ===\n');
fprintf('LO Model: Kp=%.4f, Kq=%.4f\n', K_p_LO, K_q_LO);
fprintf('HO Model: Kp1=%.4f, Kp2=%.4f | Kq1=%.4f, Kq2=%.4f\n', ...
        K_p_HO1, K_p_HO2, K_q_HO1, K_q_HO2);


%% ========================================================================
%% PHASE 2: TOURNAMENT (Validation on New Data)
%% ========================================================================
fprintf('\nPhase 2: Running Tournament on %d NEW points...\n', NumTestPoints);

rho_LO  = zeros(NumTestPoints, 1);
rho_HO  = zeros(NumTestPoints, 1);
rho_NO  = zeros(NumTestPoints, 1); % Non-Optimized
rho_OPT = zeros(NumTestPoints, 1);
x_axis  = zeros(NumTestPoints, 1);

for i = 1:NumTestPoints
    nu_val = 10^(log10(range_nu(1)) + rand()*diff(log10(range_nu)));
    w_max_val = 10^(log10(range_w_max(1)) + rand()*diff(log10(range_w_max)));
    w_min_val = range_w_min(1) + rand()*diff(range_w_min);
    
    dt_val = pi / w_max_val;
    T_val  = pi / w_min_val;
    
    % --- EVALUATE LO MODEL ---
    p_LO = (1/c) - K_p_LO * (nu_val^2 * w_max_val^2) / c^5;
    q_LO = K_q_LO * (nu_val * w_min_val * w_max_val) / c^3;
    rho_LO(i) = obj_Linf(N, T_val, dt_val, J, c, gamma, nu_val, a, M, p_LO, q_LO);
    
    % --- EVALUATE HO MODEL ---
    p_HO = (1/c) - K_p_HO1 * (nu_val^2 * w_max_val^2) / c^5 + ...
                   K_p_HO2 * (nu_val^4 * w_max_val^4) / c^9;
    p_HO = max(0.1/c, p_HO); 
    q_HO = K_q_HO1 * (nu_val * w_min_val * w_max_val) / c^3 - ...
           K_q_HO2 * (nu_val^3 * w_max_val^4) / c^7;
       
    rho_HO(i) = obj_Linf(N, T_val, dt_val, J, c, gamma, nu_val, a, M, p_HO, q_HO);
    
    % --- EVALUATE NON-OPTIMIZED (NO) ---
    rho_NO(i) = obj_Linf(N, T_val, dt_val, J, c, gamma, nu_val, a, M, 1/c, 0);
    
    % --- EVALUATE OPTIMAL ---
    objfun = @(x) obj_Linf(N, T_val, dt_val, J, c, gamma, nu_val, a, M, x(1), x(2));
    [~, fval] = fminsearch(objfun, [p_LO, q_LO], options);
    rho_OPT(i) = fval;
    
    x_axis(i) = (nu_val * w_max_val);
    
    if mod(i, 20) == 0, fprintf('.'); end
end
fprintf('\n');

%% 3. ANALYSIS & VISUALIZATION
wins_HO = sum(rho_HO < rho_LO);
loss_HO = sum(rho_HO > rho_LO);

fprintf('\n=== TOURNAMENT RESULTS ===\n');
fprintf('High Order Wins: %d (%.1f%%)\n', wins_HO, wins_HO/NumTestPoints*100);
fprintf('Mean rho_NO:  %.4f\n', mean(rho_NO));
fprintf('Mean rho_LO:  %.4f\n', mean(rho_LO));
fprintf('Mean rho_HO:  %.4f\n', mean(rho_HO));
fprintf('Mean rho_OPT: %.4f\n', mean(rho_OPT));

figure('Name','Convergence Performance','Position',[100 100 1200 500], 'Color','w');

% --- SUBPLOT 1: Scatter Comparison ---
subplot(1, 2, 1); hold on; grid on; box on;
plot([0, 1], [0, 1], 'k--', 'LineWidth', 1.5);
scatter(rho_LO, rho_HO, 60, x_axis, 'filled');
colormap(jet); c_bar = colorbar; c_bar.Label.String = '\nu \omega_{max}';
xlabel('Low Order \rho', 'FontSize', 14);
ylabel('High Order \rho', 'FontSize', 14);
title('Head-to-Head: LO vs HO', 'FontSize', 16);
max_rho = max([rho_LO; rho_HO]);
xlim([0, max_rho*1.05]); ylim([0, max_rho*1.05]);
axis square;

% --- SUBPLOT 2: Absolute Performance vs Perturbation ---
subplot(1, 2, 2); hold on; grid on; box on;
[sorted_x, idx] = sort(x_axis);

% 1. Non-Optimized (Usually High)
plot(sorted_x, rho_NO(idx), 'k:', 'LineWidth', 3, 'DisplayName', 'Non-Optimized');

% 2. Low Order
plot(sorted_x, rho_LO(idx), 'b-', 'LineWidth', 2, 'DisplayName', 'Low Order (Standard)');

% 3. High Order
plot(sorted_x, rho_HO(idx), 'r-', 'LineWidth', 2, 'DisplayName', 'High Order (New)');

% 4. Optimal (Ground Truth)
plot(sorted_x, rho_OPT(idx), 'g--', 'LineWidth', 2, 'DisplayName', 'Numerical Optimal');

xlabel('Perturbation Strength \nu \omega_{max}', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('Convergence Factor $\rho$', 'Interpreter', 'latex', 'FontSize', 14);
title('Absolute Performance Comparison', 'FontSize', 16);
legend('Location','best');

% Optional: Use Log Scale if NO is vastly worse
% set(gca, 'YScale', 'log');