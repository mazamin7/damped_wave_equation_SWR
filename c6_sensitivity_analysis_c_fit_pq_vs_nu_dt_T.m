%% CLEAN UP
clear; close all; clc;
addpath("utils\")

%% 1. CONFIGURATION
P_base = get_sim_params();
gamma = 0;
nu_fit = 1e-5; % Fix a small nu for consistent fitting

% Sweep Wave Speed c
c_vals = [0.5, 1.0, 1.5, 2.0, 3.0, 5.0];

% Storage
opt_p = zeros(size(c_vals));
opt_q = zeros(size(c_vals));

fprintf('=== SWEEP: VARYING WAVE SPEED c ===\n');

%% 2. RUN OPTIMIZATION
options = optimset('Display','off', 'TolX',1e-9, 'TolFun',1e-9);

for i = 1:length(c_vals)
    c_test = c_vals(i);
    
    % Update guess based on theory (p ~ 1/c)
    x0 = [1/c_test, 0];
    
    % Objective function with varying c
    objfun = @(x) obj_Linf(P_base.N, P_base.T, P_base.dt, P_base.J, ...
                           c_test, gamma, nu_fit, P_base.a, P_base.M, x(1), x(2));
    
    [x_opt, ~] = fminsearch(objfun, x0, options);
    
    opt_p(i) = x_opt(1);
    opt_q(i) = x_opt(2);
    
    fprintf('  c=%.1f -> p=%.5f, q=%.5f\n', c_test, opt_p(i), opt_q(i));
end

%% 3. SCALING ANALYSIS
% We expect:
% q ~ c^alpha
% p ~ c^beta

% Fit Power Laws
poly_q = polyfit(log(c_vals), log(opt_q), 1);
alpha_q = poly_q(1);

poly_p = polyfit(log(c_vals), log(opt_p), 1);
alpha_p = poly_p(1);

fprintf('\n=== SCALING RESULTS ===\n');
fprintf('q scales as c^(%.2f)  (Theory expects -3)\n', alpha_q);
fprintf('p scales as c^(%.2f)  (Theory expects -1)\n', alpha_p);

%% 4. VISUALIZATION
figure('Name','Wave Speed Scaling','Position',[100 100 1000 500], 'Color','w');

% --- q vs c ---
subplot(1, 2, 1);
loglog(c_vals, opt_q, 'b-o', 'LineWidth', 2, 'MarkerFaceColor','b');
xlabel('Wave Speed c'); ylabel('Optimal q');
title(['q Scaling: c^{' num2str(alpha_q, '%.2f') '}']);
grid on; axis square;

% --- p vs c ---
subplot(1, 2, 2);
loglog(c_vals, opt_p, 'r-s', 'LineWidth', 2, 'MarkerFaceColor','r');
xlabel('Wave Speed c'); ylabel('Optimal p');
title(['p Scaling: c^{' num2str(alpha_p, '%.2f') '}']);
grid on; axis square;

%% 5. PROPOSED UNIVERSAL FORMULA
% Based on previous results:
% q_old ~ 4.84 * nu / (dt * T)  (at c=1)
% New guess: q ~ K * nu / (dt * T * c^3) ?

% Let's test the constant K with the c-correction
% K = q * (dt * T * c^3) / nu
predicted_power_q = round(abs(alpha_q));
K_universal_q = mean(opt_q .* (P_base.dt * P_base.T * c_vals.^predicted_power_q) / nu_fit);

fprintf('\n=== FINAL PROPOSED CONSTANT ===\n');
fprintf('If q scales as c^%d, then:\n', round(alpha_q));
fprintf('q(nu, dt, T, c) = %.4f * [ nu / (dt * T * c^%d) ]\n', ...
        K_universal_q, predicted_power_q);