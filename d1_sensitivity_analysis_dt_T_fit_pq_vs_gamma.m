%% CLEAN UP
clear; close all; clc;
addpath("utils\")

%% 1. CONFIGURATION
% Physics: Telegrapher Regime (nu=0, gamma>0)
gamma = 1e-1;   
c     = 1.0;
nu    = 0;     
N     = 2; a = 0.5; M = 0; J = 1000;

% Optimization Options
options = optimset('Display','off', 'TolX',1e-12, 'TolFun',1e-12);

%% 2. EXPERIMENT A: SENSITIVITY TO W_MIN (Changing T)
% We sweep w_min from 0.05 to 5.0 rad/s
% Fixed w_max at a "safe" high value
fprintf('Experiment A: Varying w_min (Sweep T) while holding w_max fixed...\n');

num_points = 40;
w_min_vals = logspace(-1.3, 0.7, num_points); % approx 0.05 to 5.0
fixed_w_max = 200;

p_opt_A = zeros(num_points, 1);
q_opt_A = zeros(num_points, 1);

for i = 1:num_points
    w_min = w_min_vals(i);
    T_val = pi / w_min;
    dt_val = pi / fixed_w_max;
    
    objfun = @(x) obj_Linf(N, T_val, dt_val, J, c, gamma, nu, a, M, x(1), x(2));
    
    % Initial Guess (Theoretical)
    p_guess = 1/c + gamma^2 / (8*c*w_min^2);
    q_guess = gamma / (2*c);
    
    [x_opt, ~] = fminsearch(objfun, [p_guess, q_guess], options);
    p_opt_A(i) = x_opt(1);
    q_opt_A(i) = x_opt(2);
    
    if mod(i, 10)==0, fprintf('.'); end
end
fprintf('\n');

%% 3. EXPERIMENT B: SENSITIVITY TO W_MAX (Changing dt)
% We sweep w_max from 20 to 1000 rad/s
% Fixed w_min at a moderate value
fprintf('Experiment B: Varying w_max (Sweep dt) while holding w_min fixed...\n');

w_max_vals = logspace(1.3, 3, num_points); % approx 20 to 1000
fixed_w_min = 0.5;

p_opt_B = zeros(num_points, 1);
q_opt_B = zeros(num_points, 1);

for i = 1:num_points
    w_max = w_max_vals(i);
    T_val = pi / fixed_w_min;
    dt_val = pi / w_max;
    
    objfun = @(x) obj_Linf(N, T_val, dt_val, J, c, gamma, nu, a, M, x(1), x(2));
    
    % Initial Guess
    p_guess = 1/c + gamma^2 / (8*c*fixed_w_min^2);
    q_guess = gamma / (2*c);
    
    [x_opt, ~] = fminsearch(objfun, [p_guess, q_guess], options);
    p_opt_B(i) = x_opt(1);
    q_opt_B(i) = x_opt(2);
    
    if mod(i, 10)==0, fprintf('.'); end
end
fprintf('\n');

%% 4. VISUALIZATION
figure('Name','Telegrapher Sensitivity Analysis','Position',[100 100 1200 600], 'Color','w');

% --- Panel 1: Sensitivity to Simulation Time (w_min) ---
subplot(1, 2, 1); hold on; grid on; box on;
yyaxis left
plot(w_min_vals, p_opt_A, 'b-o', 'LineWidth', 2, 'MarkerFaceColor', 'b');
ylabel('Optimal $p$ (Transport)', 'Interpreter', 'latex', 'FontSize', 14, 'Color', 'b');
set(gca, 'YColor', 'b');

yyaxis right
plot(w_min_vals, q_opt_A, 'r-s', 'LineWidth', 2, 'MarkerFaceColor', 'r');
yline(gamma/(2*c), 'k--', 'Theory \gamma/2c', 'LineWidth', 2);
ylabel('Optimal $q$ (Damping)', 'Interpreter', 'latex', 'FontSize', 14, 'Color', 'r');
set(gca, 'YColor', 'r');
ylim([0, max(q_opt_A)*1.2]);

xlabel('Lower Frequency $\omega_{\min}$ ($\pi/T$)', 'Interpreter', 'latex', 'FontSize', 14);
title('Sensitivity to Simulation Time ($T$)', 'FontSize', 16);
set(gca, 'XScale', 'log');


% --- Panel 2: Sensitivity to Time Step (w_max) ---
subplot(1, 2, 2); hold on; grid on; box on;
yyaxis left
plot(w_max_vals, p_opt_B, 'b-o', 'LineWidth', 2, 'MarkerFaceColor', 'b');
ylabel('Optimal $p$ (Transport)', 'Interpreter', 'latex', 'FontSize', 14, 'Color', 'b');
set(gca, 'YColor', 'b');
% Force Y-limits to match Panel 1 range (roughly) to show scale invariance
% if p doesn't change, the line will look flat.
ylim([0.9*mean(p_opt_B), 1.1*mean(p_opt_B)]); 

yyaxis right
plot(w_max_vals, q_opt_B, 'r-s', 'LineWidth', 2, 'MarkerFaceColor', 'r');
yline(gamma/(2*c), 'k--', 'Theory \gamma/2c', 'LineWidth', 2);
ylabel('Optimal $q$ (Damping)', 'Interpreter', 'latex', 'FontSize', 14, 'Color', 'r');
set(gca, 'YColor', 'r');
ylim([0, max(q_opt_A)*1.2]); % Use same scale as left plot

xlabel('Upper Frequency $\omega_{\max}$ ($\pi/dt$)', 'Interpreter', 'latex', 'FontSize', 14);
title('Sensitivity to Time Step ($dt$)', 'FontSize', 16);
set(gca, 'XScale', 'log');

%% 5. QUANTITATIVE REPORT
fprintf('\n=== SENSITIVITY REPORT (Telegrapher Regime) ===\n');

% Calculate % change
change_p_wmin = (max(p_opt_A) - min(p_opt_A)) / min(p_opt_A) * 100;
change_q_wmin = (max(q_opt_A) - min(q_opt_A)) / min(q_opt_A) * 100;

change_p_wmax = (max(p_opt_B) - min(p_opt_B)) / mean(p_opt_B) * 100;
change_q_wmax = (max(q_opt_B) - min(q_opt_B)) / mean(q_opt_B) * 100;

fprintf('Experiment A (Varying T / w_min):\n');
fprintf('  p change: %.2f%% (Expect HUGE change)\n', change_p_wmin);
fprintf('  q change: %.2f%% (Expect small/moderate change)\n', change_q_wmin);

fprintf('\nExperiment B (Varying dt / w_max):\n');
fprintf('  p change: %.2f%% (Expect ~0%%)\n', change_p_wmax);
fprintf('  q change: %.2f%% (Expect ~0%%)\n', change_q_wmax);

fprintf('\nConclusion:\n');
if change_p_wmin > 50 && change_p_wmax < 1
    fprintf('  => p is driven entirely by w_min (Low Frequency / Time Horizon).\n');
end
if change_q_wmax < 1
    fprintf('  => q is independent of w_max (High Frequency).\n');
end