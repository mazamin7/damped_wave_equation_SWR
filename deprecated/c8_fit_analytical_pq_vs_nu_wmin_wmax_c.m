%% CLEAN UP
clear; close all; clc;
addpath("utils\")

%% 1. CONFIGURATION
% Physical constants
gamma = 0; % Fixed telegrapher damping

% Frequency definitions (Standard SWR factors)
% omega_max = pi / dt
% omega_min = pi / T

% Simulation size
N = 2; a = 0.3; M = 0; J = 1000;

%% 2. GENERATE RANDOM DATASET
NumPoints = 30; 
fprintf('Generating %d random test points...\n', NumPoints);

% Ranges for random sampling
range_nu = [1e-6, 1e-4];
range_w_max = [300, 3000];   % Corresponds roughly to dt ~ 0.01 to 0.001
range_w_min = [0.1, 1.5];    % Corresponds roughly to T ~ 30 to 2
range_c     = [0.5, 3.0];

% Storage
data_nu = zeros(NumPoints, 1);
data_wmin = zeros(NumPoints, 1);
data_wmax = zeros(NumPoints, 1);
data_c  = zeros(NumPoints, 1);
data_p  = zeros(NumPoints, 1);
data_q  = zeros(NumPoints, 1);

options = optimset('Display','off', 'TolX',1e-9, 'TolFun',1e-9);

for i = 1:NumPoints
    % 1. Sample Random Parameters
    nu_val = 10^(log10(range_nu(1)) + rand()*(log10(range_nu(2))-log10(range_nu(1))));
    c_val  = range_c(1) + rand()*(range_c(2)-range_c(1));
    
    % Sample frequencies directly, then convert to simulation params
    w_max_val = 10^(log10(range_w_max(1)) + rand()*(log10(range_w_max(2))-log10(range_w_max(1))));
    w_min_val = range_w_min(1) + rand()*(range_w_min(2)-range_w_min(1));
    
    % Convert to dt and T for the solver
    dt_val = pi / w_max_val;
    T_val  = pi / w_min_val;
    
    % Store inputs
    data_nu(i)   = nu_val;
    data_c(i)    = c_val;
    data_wmax(i) = w_max_val;
    data_wmin(i) = w_min_val;
    
    % 2. Run Optimization
    objfun = @(x) obj_Linf(N, T_val, dt_val, J, c_val, gamma, nu_val, a, M, x(1), x(2));
    
    x0 = [1/c_val, 0];
    [x_opt, ~] = fminsearch(objfun, x0, options);
    
    data_p(i) = x_opt(1);
    data_q(i) = x_opt(2);
    
    if mod(i, 10) == 0, fprintf('  Computed %d/%d...\n', i, NumPoints); end
end

%% 3. FIT FREQUENCY-BASED FORMULAS

% --- Fit q ---
% Model: q = Kq * (nu * w_max * w_min) / c^3
X_q = (data_nu .* data_wmax .* data_wmin) ./ (data_c.^3);
K_q = (X_q' * X_q) \ (X_q' * data_q);

% --- Fit p ---
% Model: p = 1/c - Kp * (nu * w_max)^2 / c^5
% Note: The previous dt term was (nu/dt)^2. Since 1/dt = w_max/pi,
% the fitted Kp here will absorb the factor (1/pi)^2.
target_p_diff = (1 ./ data_c) - data_p;
X_p = (data_nu .* data_wmax).^2 ./ (data_c.^5);
K_p = (X_p' * X_p) \ (X_p' * target_p_diff);

fprintf('\n=== FITTED FREQUENCY CONSTANTS ===\n');
fprintf('q Formula: q = %.4f * [ nu * w_max * w_min / c^3 ]\n', K_q);
fprintf('p Formula: p = 1/c - %.4f * [ (nu * w_max)^2 / c^5 ]\n', K_p);

%% 4. VISUALIZATION AND VALIDATION
figure('Name','Frequency-Based Fit','Position',[100 100 1000 500], 'Color','w');

% --- Q Validation ---
subplot(1, 2, 1); hold on; grid on; box on;
q_pred = K_q * X_q;
plot(data_q, q_pred, 'bo', 'MarkerFaceColor','b', 'MarkerSize', 6);
% 1:1 Line
max_q = max([data_q; q_pred]); min_q = min([data_q; q_pred]);
plot([min_q max_q], [min_q max_q], 'k--', 'LineWidth', 2);
xlabel('Numerical Optimization'); ylabel('Formula Prediction');
title('q vs Frequency Parameters');
subtitle(['K_q \approx ' num2str(K_q, '%.4f')]);
axis square;

% --- P Validation ---
subplot(1, 2, 2); hold on; grid on; box on;
p_pred = (1 ./ data_c) - K_p * X_p;
plot(data_p, p_pred, 'ro', 'MarkerFaceColor','r', 'MarkerSize', 6);
% 1:1 Line
max_p = max([data_p; p_pred]); min_p = min([data_p; p_pred]);
plot([min_p max_p], [min_p max_p], 'k--', 'LineWidth', 2);
xlabel('Numerical Optimization'); ylabel('Formula Prediction');
title('p vs Frequency Parameters');
subtitle(['K_p \approx ' num2str(K_p, '%.4f')]);
axis square;

%% 5. ANALYTICAL CONTEXT
% Check consistency with previous dt/T results
% We know from previous step: K_q_old approx 4.84
% Formula conversion:
% Old: q = 4.84 * nu / (dt * T * c^3)
% New: q = K_q  * nu * (pi/dt) * (pi/T) / c^3
% Therefore: K_q * pi^2 should approx 4.84
theoretical_Kq = 4.8413 / (pi^2);
fprintf('\n=== CONSISTENCY CHECK ===\n');
fprintf('Theoretical K_q (derived from previous fit): %.4f\n', theoretical_Kq);
fprintf('Actual Fitted K_q: %.4f\n', K_q);