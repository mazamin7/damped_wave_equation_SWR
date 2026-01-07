%% SCRIPT 3: TELEGRAPHER SCALING VALIDATION (ROBUST)
clear; close all; clc; addpath("utils\");

c = 1.0; nu = 0; N = 2; a = 0.5; M = 0; J = 1000;

% FIX: Tight tolerances for precise scaling checks
options = optimset('Display','off', 'TolX',1e-12, 'TolFun',1e-12);

fprintf('--- TELEGRAPHER SCALING CHECKS ---\n');

%% TEST 1: DEPENDENCE ON GAMMA (Fix Frequencies)
fprintf('1. Checking Scaling wrt Gamma...\n');
gamma_vals = logspace(-4, -2, 20); % 0.01 to 1.0
w_min = 1.0; 
w_max = 50.0;
T = pi/w_min; 
dt = pi/w_max;

data_p_corr = zeros(size(gamma_vals));
data_q      = zeros(size(gamma_vals));
data_rho    = zeros(size(gamma_vals));

for i=1:length(gamma_vals)
    gamma = gamma_vals(i);
    
    % Robust Initial Guess based on refined theory
    p_guess = 1/c + gamma^2/(16*c*w_min^2);
    q_guess = gamma/(2*c);
    
    obj = @(x) obj_Linf(N,T,dt,J,c,gamma,nu,a,M,x(1),x(2));
    [res, fval] = fminsearch(obj, [p_guess, q_guess], options);
    
    data_p_corr(i) = abs(res(1) - 1/c); % Positive addition
    data_q(i)      = abs(res(2));
    data_rho(i)    = fval;
end

% Robust Fitting
valid_q = data_q > 1e-15;
coeffs_q = polyfit(log(gamma_vals(valid_q)), log(data_q(valid_q)), 1);

valid_p = data_p_corr > 1e-15;
coeffs_p = polyfit(log(gamma_vals(valid_p)), log(data_p_corr(valid_p)), 1);

valid_rho = data_rho > 1e-15;
coeffs_rho = polyfit(log(gamma_vals(valid_rho)), log(data_rho(valid_rho)), 1);

fprintf('   q power wrt gamma:   %.2f (Expect 1.0)\n', coeffs_q(1));
fprintf('   p corr power wrt gamma: %.2f (Expect 2.0)\n', coeffs_p(1));
fprintf('   rho power wrt gamma: %.2f (Expect 2.0)\n', coeffs_rho(1));

%% TEST 2: DEPENDENCE ON W_MIN (Fix Gamma)
fprintf('2. Checking Scaling wrt w_min...\n');
gamma_fix = 1e-4;
w_min_vals = logspace(-1, 0.5, 20); % 0.1 to 3
w_max_fix = 100;

data_p_w_corr = zeros(size(w_min_vals));
data_q_w      = zeros(size(w_min_vals));
data_rho_w    = zeros(size(w_min_vals));

for i=1:length(w_min_vals)
    w_min = w_min_vals(i);
    T_var = pi/w_min;
    
    % Theoretical Guess
    p_guess = 1/c + gamma_fix^2/(16*c*w_min^2);
    q_guess = gamma_fix/(2*c);
    
    obj = @(x) obj_Linf(N,T_var,pi/w_max_fix,J,c,gamma_fix,nu,a,M,x(1),x(2));
    [res, fval] = fminsearch(obj, [p_guess, q_guess], options);
    
    data_p_w_corr(i) = abs(res(1) - 1/c);
    data_q_w(i)      = abs(res(2));
    data_rho_w(i)    = fval;
end

% Robust Fitting
valid_p_w = data_p_w_corr > 1e-15;
coeffs_p_w = polyfit(log(w_min_vals(valid_p_w)), log(data_p_w_corr(valid_p_w)), 1);

valid_q_w = data_q_w > 1e-15;
% Note: q should be constant, so slope ~0. We check variation instead.
q_slope = polyfit(log(w_min_vals(valid_q_w)), log(data_q_w(valid_q_w)), 1);

valid_rho_w = data_rho_w > 1e-15;
coeffs_rho_w = polyfit(log(w_min_vals(valid_rho_w)), log(data_rho_w(valid_rho_w)), 1);

fprintf('   p corr power wrt w_min: %.2f (Expect -2.0)\n', coeffs_p_w(1));
fprintf('   q power wrt w_min:      %.2f (Expect ~0.0)\n', q_slope(1));
fprintf('   rho power wrt w_min:    %.2f (Expect -2.0)\n', coeffs_rho_w(1));

%% VISUALIZATION WITH THEORY CURVES
figure('Name','Tele Scaling Validation','Color','w', 'Position', [100 100 1400 500]);

% --- SUBPLOT 1: SCALING WITH GAMMA ---
subplot(1,3,1); 
% Data
loglog(gamma_vals, data_p_corr, 'bo', 'MarkerFaceColor', 'b', 'DisplayName','p_{corr} Data');
hold on;
loglog(gamma_vals, data_q, 'rs', 'MarkerFaceColor', 'r', 'DisplayName','q Data');

% Theory
% p_corr ~ gamma^2
% q ~ gamma^1
theory_p_g = data_p_corr(end) * (gamma_vals / gamma_vals(end)).^2;
theory_q_g = data_q(end) * (gamma_vals / gamma_vals(end)).^1;

loglog(gamma_vals, theory_p_g, 'b--', 'LineWidth', 1.5, 'DisplayName','Theory \gamma^2');
loglog(gamma_vals, theory_q_g, 'r--', 'LineWidth', 1.5, 'DisplayName','Theory \gamma^1');

title('Scaling wrt \gamma'); xlabel('\gamma'); 
legend('Location','northwest'); grid on; box on;

% --- SUBPLOT 2: P & Q SCALING WITH W_MIN ---
subplot(1,3,2); 
% Data
loglog(w_min_vals, data_p_w_corr, 'bo', 'MarkerFaceColor', 'b', 'DisplayName','p_{corr} Data');
hold on;
loglog(w_min_vals, data_q_w, 'rs', 'MarkerFaceColor', 'r', 'DisplayName','q Data');

% Theory
% p_corr ~ w_min^-2
% q ~ constant
theory_p_w = data_p_w_corr(1) * (w_min_vals / w_min_vals(1)).^(-2);
theory_q_w = repmat(mean(data_q_w), size(w_min_vals));

loglog(w_min_vals, theory_p_w, 'b--', 'LineWidth', 1.5, 'DisplayName','Theory \omega_{min}^{-2}');
loglog(w_min_vals, theory_q_w, 'r--', 'LineWidth', 1.5, 'DisplayName','Theory Const');

title('Parameters wrt \omega_{min}'); xlabel('\omega_{min}'); 
legend('Location','southwest'); grid on; box on;

% --- SUBPLOT 3: RHO SCALING WITH W_MIN ---
subplot(1,3,3); 
% Data
loglog(w_min_vals, data_rho_w, 'kd', 'MarkerFaceColor', 'k', 'DisplayName','\rho Data');
hold on;

% Theory: rho ~ w_min^-2 (Assuming gamma << w_min regime holds roughly, or refined formula)
% In refined formula: rho ~ sqrt(4w^2)/w^3 ~ 1/w^2.
theory_rho_w = data_rho_w(1) * (w_min_vals / w_min_vals(1)).^(-2);

loglog(w_min_vals, theory_rho_w, 'g--', 'LineWidth', 2, 'DisplayName','Theory \omega_{min}^{-2}');

title('Convergence \rho wrt \omega_{min}'); xlabel('\omega_{min}'); 
legend('Location','southwest'); grid on; box on;