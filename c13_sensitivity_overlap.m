%% CLEAN UP
clear; close all; clc;
addpath("utils\")

%% 1. CONFIGURATION
% We will sweep M from 0 (no overlap) to 32 (large overlap)
M_values = [0, 0.001, 0.002, 0.005, 0.01]; 
num_M = length(M_values);

% Physics constraints (Strictly valid regime)
range_nu    = [1e-6, 1e-4]; 
range_w_max = [100, 1000];   
range_w_min = [0.1, 2.0];
c = 1.0; 
gamma = 0;
N = 2; a = 0.5; J = 1000;

% Simulation size per M
NumPointsPerM = 50; 
fprintf('Analyzing sensitivity to Overlap M (%d variations, %d points each)...\n', num_M, NumPointsPerM);

% Storage for fitted constants
K_q_results   = zeros(num_M, 1);
K_p_results   = zeros(num_M, 1);
K_rho_results = zeros(num_M, 1);
rho_abs_vals  = zeros(num_M, 1); % To track if rho actually gets better

options = optimset('Display','off', 'TolX',1e-9, 'TolFun',1e-9);

%% 2. SENSITIVITY LOOP
for k = 1:num_M
    current_M = M_values(k);
    fprintf('  Processing M = %d... ', current_M);
    
    % Temp storage for this batch
    batch_Kq = zeros(NumPointsPerM, 1);
    batch_Kp = zeros(NumPointsPerM, 1);
    batch_Krho = zeros(NumPointsPerM, 1);
    batch_rho = zeros(NumPointsPerM, 1);
    
    for i = 1:NumPointsPerM
        % Randomize Physics
        nu = 10^(log10(range_nu(1)) + rand()*diff(log10(range_nu)));
        w_max = 10^(log10(range_w_max(1)) + rand()*diff(log10(range_w_max)));
        w_min = range_w_min(1) + rand()*diff(range_w_min);
        
        dt = pi / w_max;
        T  = pi / w_min;
        
        % Optimize
        objfun = @(x) obj_Linf(N, T, dt, J, c, gamma, nu, a, current_M, x(1), x(2));
        
        % Guess (using standard constants)
        q0 = 0.5 * nu * w_min * w_max / c^3;
        p0 = 1/c - 0.25 * nu^2 * w_max^2 / c^5;
        
        [x_opt, rho_val] = fminsearch(objfun, [p0, q0], options);
        p_opt = x_opt(1);
        q_opt = x_opt(2);
        
        % --- Calculate Local Constants ---
        
        % 1. K_q (Expected 0.50)
        % Formula: q = K * nu * w_min * w_max / c^3
        term_q = (nu * w_min * w_max) / c^3;
        batch_Kq(i) = q_opt / term_q;
        
        % 2. K_p (Expected 0.25)
        % Formula: (1/c - p) = K * nu^2 * w_max^2 / c^5
        term_p = (nu^2 * w_max^2) / c^5;
        drop_p = (1/c) - p_opt;
        batch_Kp(i) = drop_p / term_p;
        
        % 3. K_rho (Expected 0.25)
        % Formula: rho = K * nu * (w_max - w_min) / c^2
        term_rho = (nu * (w_max - w_min)) / c^2;
        batch_Krho(i) = rho_val / term_rho;
        
        batch_rho(i) = rho_val;
    end
    
    % Store Robust Means
    K_q_results(k)   = mean(batch_Kq);
    K_p_results(k)   = mean(batch_Kp);
    K_rho_results(k) = mean(batch_Krho);
    rho_abs_vals(k)  = mean(batch_rho);
    
    fprintf('Done. K_q=%.2f, K_p=%.2f, K_rho=%.2f\n', ...
            K_q_results(k), K_p_results(k), K_rho_results(k));
end

%% 3. VISUALIZATION
figure('Name','Sensitivity to Overlap M','Position',[100 100 1200 400], 'Color','w');

% --- K_q Plot ---
subplot(1, 3, 1); hold on; grid on; box on;
plot(M_values, K_q_results, 'bo-', 'LineWidth', 2, 'MarkerFaceColor', 'b');
yline(0.50, 'k--', 'Theoretical (0.50)', 'LineWidth', 2);
xlabel('Overlap M (grid points)', 'FontSize', 12);
ylabel('Fitted Constant K_q', 'FontSize', 12);
title('Stability of Damping Param q', 'FontSize', 14);
ylim([0.4, 0.6]);

% --- K_p Plot ---
subplot(1, 3, 2); hold on; grid on; box on;
plot(M_values, K_p_results, 'rs-', 'LineWidth', 2, 'MarkerFaceColor', 'r');
yline(0.25, 'k--', 'Theoretical (0.25)', 'LineWidth', 2);
xlabel('Overlap M (grid points)', 'FontSize', 12);
ylabel('Fitted Constant K_p', 'FontSize', 12);
title('Stability of Transport Param p', 'FontSize', 14);
ylim([0.2, 0.3]);

% --- K_rho Plot ---
subplot(1, 3, 3); hold on; grid on; box on;
[ax, h1, h2] = plotyy(M_values, K_rho_results, M_values, rho_abs_vals);

% Style Left Axis (Constant K)
set(h1, 'LineWidth', 2, 'Marker', 'd', 'MarkerFaceColor', 'g', 'Color', [0 0.6 0]);
ylabel(ax(1), 'Fitted Constant K_{\rho}', 'FontSize', 12, 'Color', [0 0.6 0]);
set(ax(1), 'YColor', [0 0.6 0], 'YLim', [0.2, 0.3]);
yline(ax(1), 0.25, 'k--', 'LineWidth', 2);

% Style Right Axis (Absolute Rho)
set(h2, 'LineWidth', 2, 'LineStyle', '--', 'Color', 'k');
ylabel(ax(2), 'Absolute \rho Value (Mean)', 'FontSize', 12, 'Color', 'k');
set(ax(2), 'YColor', 'k');

xlabel('Overlap M (grid points)', 'FontSize', 12);
title('Convergence vs Overlap', 'FontSize', 14);

legend([h1, h2], {'Scaling Constant K', 'Actual \rho Value'}, 'Location', 'north');