%% CLEAN UP
clear; close all; clc;
addpath("utils\")

%% 1. CONFIGURATION
gamma = 0;
nu    = 1e-6;  
c     = 1.0;   

% Simulation Geometry
N = 2; a = 0.5; M = 0; J = 1000;

%% 2. GENERATE FREQUENCY DATASET
NumPoints = 100; 
fprintf('Generating %d test points...\n', NumPoints);

% Ranges
range_w_max = [300, 3000];
range_w_min = [0.1, 2.0];

data_wmax = zeros(NumPoints, 1);
data_wmin = zeros(NumPoints, 1);
data_p    = zeros(NumPoints, 1);
data_q    = zeros(NumPoints, 1);

options = optimset('Display','off', 'TolX',1e-9, 'TolFun',1e-9);

for i = 1:NumPoints
    % 1. Randomize Frequencies
    w_max_val = 10^(log10(range_w_max(1)) + rand()*diff(log10(range_w_max)));
    w_min_val = range_w_min(1) + rand()*diff(range_w_min);
    
    dt_val = pi / w_max_val;
    T_val  = pi / w_min_val;
    
    data_wmax(i) = w_max_val;
    data_wmin(i) = w_min_val;
    
    % 2. Run Optimization
    objfun = @(x) obj_Linf(N, T_val, dt_val, J, c, gamma, nu, a, M, x(1), x(2));
    [x_opt, ~] = fminsearch(objfun, [1/c, 0], options);
    
    data_p(i) = x_opt(1);
    data_q(i) = x_opt(2);
    
    if mod(i, 10) == 0, fprintf('  Computed %d/%d...\n', i, NumPoints); end
end

%% 3. ROBUST FITTING (IGNORING OUTLIERS)

% --- Fit q (Usually stable, but we filter just in case) ---
term_q = (nu / c^3) * (data_wmin .* data_wmax);
K_q_dist = data_q ./ term_q;

% --- Fit p (High outlier risk due to resonance) ---
term_p = (nu^2 / c^5) * (data_wmax.^2);
diff_p = (1/c) - data_p;
K_p_dist = diff_p ./ term_p;

% Detect outliers in the calculated constants
% 'quartiles' method is robust to extreme values
idx_out_p = isoutlier(K_p_dist, 'quartiles');
idx_in_p  = ~idx_out_p;

% Calculate mean only on INLIERS
K_q_mean = mean(K_q_dist); % q is usually stable, use all or filter if needed
K_p_mean = mean(K_p_dist(idx_in_p));

fprintf('\n=== FINAL ROBUST CONSTANTS ===\n');
fprintf('Target Model for q: q = K_q * [ nu * w_min * w_max / c^3 ]\n');
fprintf('  -> Fitted K_q: %.4f  (Expect ~0.50)\n', K_q_mean);

fprintf('\nTarget Model for p: p = 1/c - K_p * [ (nu * w_max)^2 / c^5 ]\n');
fprintf('  -> Fitted K_p: %.4f  (Calculated from %d inliers, ignored %d outliers)\n', ...
        K_p_mean, sum(idx_in_p), sum(idx_out_p));

%% 4. VISUALIZATION
figure('Name','Robust Frequency Fit','Position',[100 100 1000 500], 'Color','w');

% --- Q Plot ---
subplot(1, 2, 1); hold on; grid on; box on;
plot(term_q, data_q, 'bo', 'MarkerFaceColor','b', 'MarkerSize', 6);
plot([min(term_q) max(term_q)], K_q_mean * [min(term_q) max(term_q)], 'r--', 'LineWidth', 2);
xlabel('Theoretical Term: $\frac{\nu \omega_{\min} \omega_{\max}}{c^3}$', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('Optimal q', 'FontSize', 14);
title(['q Fit (K_q = ' num2str(K_q_mean, '%.3f') ')'], 'FontSize', 14);
axis square;

% --- P Plot (With Outliers Highlighted) ---
subplot(1, 2, 2); hold on; grid on; box on;
x_axis_p = data_wmax.^2;
y_axis_p = (1/c) - data_p;

% Plot Outliers (Gray X)
if any(idx_out_p)
    plot(x_axis_p(idx_out_p), y_axis_p(idx_out_p), 'x', 'Color', [0.7 0.7 0.7], ...
         'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'Outliers (Ignored)');
end

% Plot Inliers (Red Circles)
plot(x_axis_p(idx_in_p), y_axis_p(idx_in_p), 'ro', 'MarkerFaceColor','r', ...
     'MarkerSize', 6, 'DisplayName', 'Valid Data');

% Plot Robust Fit Line
x_fit = linspace(min(x_axis_p), max(x_axis_p), 100);
y_fit = K_p_mean * (nu^2/c^5) * x_fit;
plot(x_fit, y_fit, 'b--', 'LineWidth', 2, 'DisplayName', 'Robust Fit');

xlabel('$\omega_{\max}^2$', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('Correction: $(1/c - p)$', 'Interpreter', 'latex', 'FontSize', 14);
title(['p Fit (K_p = ' num2str(K_p_mean, '%.3f') ')'], 'FontSize', 14);
legend('Location', 'best');
axis square;