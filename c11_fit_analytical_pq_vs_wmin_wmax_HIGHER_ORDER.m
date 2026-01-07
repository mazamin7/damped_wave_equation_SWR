%% CLEAN UP
clear; close all; clc;
addpath("utils\")

%% 1. CONFIGURATION
gamma = 0;
c     = 1.0;   
N = 2; a = 0.5; M = 0; J = 1000;

%% 2. GENERATE DATA (VALID ASYMPTOTIC REGIME)
NumPoints = 200; 
fprintf('Generating %d test points (Valid Regime: nu*w < c^2)...\n', NumPoints);

% Ranges (Wide enough to see both linear and non-linear behaviors)
range_nu    = [1e-5, 1e-3];  
range_w_max = [10, 50];      
range_w_min = [0.1, 1.0];

% Storage
data_nu = zeros(NumPoints, 1);
data_wmax = zeros(NumPoints, 1);
data_wmin = zeros(NumPoints, 1);
data_p = zeros(NumPoints, 1);
data_q = zeros(NumPoints, 1);

options = optimset('Display','off', 'TolX',1e-12, 'TolFun',1e-12);

for i = 1:NumPoints
    % --- RANDOM SAMPLING (LOGARITHMIC) ---
    % Sample nu uniformly in log-space: 10^( min + rand * (max-min) )
    nu_val = 10^(log10(range_nu(1)) + rand()*diff(log10(range_nu)));
    
    % Sample w_max uniformly in log-space (standard for frequency analysis)
    w_max_val = 10^(log10(range_w_max(1)) + rand()*diff(log10(range_w_max)));
    
    % w_min can remain linear (small dynamic range)
    w_min_val = range_w_min(1) + rand()*diff(range_w_min);
    
    dt_val = pi / w_max_val;
    T_val  = pi / w_min_val;
    
    data_nu(i) = nu_val;
    data_wmax(i) = w_max_val;
    data_wmin(i) = w_min_val;
    
    % Optimize
    objfun = @(x) obj_Linf(N, T_val, dt_val, J, c, gamma, nu_val, a, M, x(1), x(2));
    
    % Initial Guess
    q_guess = 0.5 * nu_val * w_min_val * w_max_val / c^3;
    p_guess = 1/c - 0.25 * nu_val^2 * w_max_val^2 / c^5;
    
    [x_opt, ~] = fminsearch(objfun, [p_guess, q_guess], options);
    data_p(i) = x_opt(1);
    data_q(i) = x_opt(2);
    
    if mod(i, 20) == 0, fprintf('.'); end
end
fprintf('\n');

%% 3. SIMULTANEOUS MULTIVARIATE FITTING

% --- Fit q (Damping) ---
% Model: q = K_q1 * Term_LO + K_q3 * Term_HO
% Term_LO = nu * w_min * w_max / c^3
% Term_HO = nu^3 * w_max^4 / c^7 (Hypothesis)
X_q_LO = (data_nu ./ c^3) .* (data_wmin .* data_wmax);
X_q_HO = (data_nu.^3 ./ c^7) .* data_wmax.^4;

% Design Matrix [Col1, Col2]
A_q = [X_q_LO, X_q_HO];
b_q = data_q;

% Solve A*k = b
k_q = A_q \ b_q;
K_q1_fit = k_q(1);
K_q3_fit = k_q(2);


% --- Fit p (Transport) ---
% Model: p = 1/c - K_p2 * Term_LO + K_p4 * Term_HO
% Rearrange: (1/c - p) = K_p2 * Term_LO - K_p4 * Term_HO
% Let y = (1/c - p)
% We fit y = beta1 * Term_LO + beta2 * Term_HO
% Then K_p2 = beta1, and K_p4 = -beta2

X_p_LO = (data_nu.^2 ./ c^5) .* data_wmax.^2;
X_p_HO = (data_nu.^4 ./ c^9) .* data_wmax.^4;

% Design Matrix
A_p = [X_p_LO, X_p_HO];
b_p = (1/c) - data_p;

% Solve A*k = b
k_p = A_p \ b_p;
K_p2_fit = k_p(1);      % The Leading Order coeff
K_p4_fit = -k_p(2);     % The Higher Order coeff (negative slope in y means positive addition to p)


%% 4. RESULTS DISPLAY
fprintf('\n=== SIMULTANEOUS FIT RESULTS ===\n');

fprintf('\n--- Parameter q (Damping) ---\n');
fprintf('Model: q = K1 * (nu*w*w) + K3 * (nu^3*w^4)\n');
fprintf('  K1 (Leading Order): %.4f  (Expect 0.5000)\n', K_q1_fit);
fprintf('  K3 (Higher Order):  %.4f  (Expect ~0.0000)\n', K_q3_fit);

fprintf('\n--- Parameter p (Transport) ---\n');
fprintf('Model: p = 1/c - K2 * (nu^2*w^2) + K4 * (nu^4*w^4)\n');
fprintf('  K2 (Leading Order): %.4f  (Expect 0.2500)\n', K_p2_fit);
fprintf('  K4 (Higher Order):  %.4f  (Expect ~0.10 - 0.15)\n', K_p4_fit);


%% 5. VISUALIZATION
figure('Name','Simultaneous Fit Analysis','Position',[100 100 1200 500], 'Color','w');

% --- Q Plot ---
subplot(1, 2, 1); hold on; grid on; box on;
% Calculate predicted values based on fit
q_pred = A_q * k_q;
% Plot Correlation
plot(data_q, q_pred, 'bo', 'MarkerFaceColor','b', 'MarkerSize', 5, 'DisplayName','Data');
plot([min(data_q) max(data_q)], [min(data_q) max(data_q)], 'k--', 'LineWidth', 2, 'DisplayName','Perfect Fit');
xlabel('Observed Optimal q', 'FontSize', 14);
ylabel('Predicted q (LO + HO)', 'FontSize', 14);
title(['q Fit: K_{LO}=' num2str(K_q1_fit,'%.3f') ', K_{HO}=' num2str(K_q3_fit,'%.5f')], 'FontSize', 14);
legend('Location','best');
axis square;

% --- P Plot ---
subplot(1, 2, 2); hold on; grid on; box on;
% Calculate predicted "drop" (1/c - p)
y_pred = A_p * k_p;
y_obs  = b_p;

% Plot vs the main driver (LO Term) to see deviations
plot(X_p_LO, y_obs, 'ro', 'MarkerFaceColor','r', 'MarkerSize', 5, 'DisplayName','Observed Drop');
[sorted_x, idx] = sort(X_p_LO);
plot(sorted_x, y_pred(idx), 'b-', 'LineWidth', 2, 'DisplayName','Combined Model Fit');
plot(sorted_x, K_p2_fit * sorted_x, 'k:', 'LineWidth', 2, 'DisplayName','Linear Only (LO)');

xlabel('Leading Driver $\nu^2 \omega_{\max}^2 / c^5$', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('Value of $(1/c - p)$', 'Interpreter', 'latex', 'FontSize', 14);
title(['p Fit: K_{LO}=' num2str(K_p2_fit,'%.3f') ', K_{HO}=' num2str(K_p4_fit,'%.3f')], 'FontSize', 14);
legend('Location','best');
axis square;