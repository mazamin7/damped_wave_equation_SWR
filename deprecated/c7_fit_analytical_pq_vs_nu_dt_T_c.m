%% CLEAN UP
clear; close all; clc;
addpath("utils\")

%% 1. GLOBAL CONSTANTS (From your sensitivity analysis)
K_q = 4.8413;   % Empirical constant for q
K_p = 2.3830;   % Empirical constant for p correction

%% 2. GENERATE RANDOM TEST SCENARIOS
% We will create 30 random simulation setups
NumTests = 30;

% Random ranges
range_nu = [1e-6, 1e-4];
range_dt = [0.001, 0.01];
range_T  = [2, 10];
range_c  = [0.5, 3.0]; % Now varying wave speed too!

% Storage
pred_q = zeros(NumTests, 1);
real_q = zeros(NumTests, 1);
pred_p = zeros(NumTests, 1);
real_p = zeros(NumTests, 1);

fprintf('=== UNIVERSAL FORMULA VALIDATION ===\n');
fprintf('Testing %d random configurations...\n', NumTests);

options = optimset('Display','off', 'TolX',1e-9, 'TolFun',1e-9);

for i = 1:NumTests
    % Randomized parameters
    nu = 10^( log10(range_nu(1)) + rand()*(log10(range_nu(2))-log10(range_nu(1))) );
    dt = 10^( log10(range_dt(1)) + rand()*(log10(range_dt(2))-log10(range_dt(1))) );
    T  = range_T(1) + rand()*(range_T(2)-range_T(1));
    c  = range_c(1) + rand()*(range_c(2)-range_c(1));
    
    % Fixed parameters
    N = 2; a = 0.3; M = 0; J = 1000; gamma = 0;
    
    % --- 1. CALCULATE ANALYTICAL PREDICTION ---
    % q ~ Kq * nu / (c^3 * dt * T)
    pred_q(i) = K_q * nu / (c^3 * dt * T);
    
    % p ~ 1/c - Kp * (nu / dt)^2 / c^5   <-- Theoretical c^5 correction scaling
    pred_p(i) = (1/c) - (K_p/c^5) * (nu/dt)^2;
    
    % --- 2. COMPUTE NUMERICAL OPTIMUM ---
    objfun = @(x) obj_Linf(N, T, dt, J, c, gamma, nu, a, M, x(1), x(2));
    x0 = [1/c, 0];
    [x_opt, ~] = fminsearch(objfun, x0, options);
    
    real_p(i) = x_opt(1);
    real_q(i) = x_opt(2);
    
    if mod(i,5)==0, fprintf('  Test %d/%d done.\n', i, NumTests); end
end

%% 3. VISUALIZATION
figure('Name','Universal Validation','Position',[100 100 1000 500], 'Color','w');

% --- Q Validation ---
subplot(1, 2, 1); hold on; grid on; box on;
plot(real_q, pred_q, 'bo', 'MarkerFaceColor','b', 'MarkerSize', 6);
% Draw perfect fit line
max_val = max([real_q; pred_q]);
min_val = min([real_q; pred_q]);
plot([min_val max_val], [min_val max_val], 'k--', 'LineWidth', 2);
xlabel('Numerical Optimization'); ylabel('Analytical Formula');
title('Validation of q Formula');
subtitle('q \approx 4.84 \cdot \nu / (c^3 \Delta t T)');
axis square;

% --- P Validation ---
subplot(1, 2, 2); hold on; grid on; box on;
plot(real_p, pred_p, 'ro', 'MarkerFaceColor','r', 'MarkerSize', 6);
% Draw perfect fit line
max_val = max([real_p; pred_p]);
min_val = min([real_p; pred_p]);
plot([min_val max_val], [min_val max_val], 'k--', 'LineWidth', 2);
xlabel('Numerical Optimization'); ylabel('Analytical Formula');
title('Validation of p Formula');
subtitle('p \approx 1/c - 2.38/c^5 \cdot (\nu/\Delta t)^2');
axis square;

%% 4. ERROR METRICS
mape_q = mean(abs((real_q - pred_q) ./ real_q)) * 100;
mape_p = mean(abs((real_p - pred_p) ./ real_p)) * 100;

fprintf('\n=== ACCURACY REPORT ===\n');
fprintf('Mean Absolute Percentage Error (MAPE):\n');
fprintf('  q Formula: %.2f%%\n', mape_q);
fprintf('  p Formula: %.2f%%\n', mape_p);