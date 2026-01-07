%% CLEAN UP
clear; close all; clc;
addpath("utils\")

%% 1. CONFIGURATION
P_base = get_sim_params();
c = P_base.c; 
gamma = 0; % Fixed telegrapher damping
M = P_base.M; % Keep overlap fixed (standard SWR requirement)

% Global Fitted Constants to Test against
Cq_global = 4.8413;
Cp_global = 2.3830;

% Parameter Ranges
a_vals = 0.2 : 0.1 : 2;  % Sweep interface position from 20% to 70% of domain
nu_range = logspace(-8, -4, 50); % Small nu range for fitting

% Preallocate storage for sensitivity results
Cq_sensitivity = zeros(size(a_vals));
Cp_sensitivity = zeros(size(a_vals));

fprintf('=== SENSITIVITY ANALYSIS: VARYING PARAMETER a ===\n');
fprintf('Testing stability of Global Constants (Cq=%.2f, Cp=%.2f)\n', Cq_global, Cp_global);

%% 2. MAIN LOOP OVER POSITIONS 'a'
options = optimset('Display','off', 'TolX',1e-8, 'TolFun',1e-8);
x0 = [1/c, 0]; 

for i = 1:length(a_vals)
    a_curr = a_vals(i);
    
    % Temporary storage for this specific 'a'
    p_opt_local = zeros(size(nu_range));
    q_opt_local = zeros(size(nu_range));
    
    % Compute optima for a range of nu at this position
    for k = 1:length(nu_range)
        nu = nu_range(k);
        % Construct objective with current 'a'
        objfun = @(x) obj_Linf(P_base.N, P_base.T, P_base.dt, P_base.J, ...
                               c, gamma, nu, a_curr, M, x(1), x(2));
        
        [x_opt, ~] = fminsearch(objfun, x0, options);
        p_opt_local(k) = x_opt(1);
        q_opt_local(k) = x_opt(2);
        x0 = x_opt; % Warm start
    end
    
    % --- RE-FIT CONSTANTS FOR THIS POSITION ---
    % Formula: q = Cq * [nu / (dt * T)]
    % Therefore: Cq = q / [nu / (dt * T)]
    denominator_q = nu_range ./ (P_base.dt * P_base.T);
    % Use mean of point-wise division for robustness (or linear regression slope)
    Cq_sensitivity(i) = mean(q_opt_local ./ denominator_q);
    
    % Formula: p = 1/c - Cp * [nu / dt]^2
    % Therefore: Cp = (1/c - p) / [nu / dt]^2
    denominator_p = (nu_range ./ P_base.dt).^2;
    delta_p = (1/c) - p_opt_local;
    Cp_sensitivity(i) = mean(delta_p ./ denominator_p);
    
    fprintf('  a=%.2f -> Cq_local=%.4f, Cp_local=%.4f\n', a_curr, Cq_sensitivity(i), Cp_sensitivity(i));
end

%% 3. VISUALIZATION
figure('Name','Sensitivity to Position a','Position',[100 100 1000 500], 'Color','w');

% --- Cq Sensitivity Plot ---
subplot(1, 2, 1); hold on; grid on; box on;
plot(a_vals, Cq_sensitivity, 'b-o', 'LineWidth', 2, 'MarkerFaceColor', 'b');
yline(Cq_global, 'k--', 'Global Fit (4.84)', 'LineWidth', 2);
xlabel('Interface Position a');
ylabel('Fitted Constant C_q');
title('Sensitivity of C_q to Geometry');
ylim([min(Cq_sensitivity)*0.9, max(Cq_sensitivity)*1.1]);

% --- Cp Sensitivity Plot ---
subplot(1, 2, 2); hold on; grid on; box on;
plot(a_vals, Cp_sensitivity, 'r-s', 'LineWidth', 2, 'MarkerFaceColor', 'r');
yline(Cp_global, 'k--', 'Global Fit (2.38)', 'LineWidth', 2);
xlabel('Interface Position a');
ylabel('Fitted Constant C_p');
title('Sensitivity of C_p to Geometry');
ylim([min(Cp_sensitivity)*0.9, max(Cp_sensitivity)*1.1]);

%% 4. STATISTICAL SUMMARY
mean_Cq = mean(Cq_sensitivity);
std_Cq  = std(Cq_sensitivity);
mean_Cp = mean(Cp_sensitivity);
std_Cp  = std(Cp_sensitivity);

fprintf('\n=== SUMMARY ===\n');
fprintf('Cq: Mean=%.4f, Std=%.4f (Variation: %.2f%%)\n', ...
        mean_Cq, std_Cq, (std_Cq/mean_Cq)*100);
fprintf('Cp: Mean=%.4f, Std=%.4f (Variation: %.2f%%)\n', ...
        mean_Cp, std_Cp, (std_Cp/mean_Cp)*100);

if (std_Cq/mean_Cq) < 0.05
    fprintf('\nCONCLUSION: Parameters are largely INDEPENDENT of position a.\n');
else
    fprintf('\nCONCLUSION: Parameters show significant dependence on position a.\n');
end