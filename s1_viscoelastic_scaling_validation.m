%% SCRIPT 1: VISCOELASTIC SCALING VALIDATION (ROBUST)
clear; close all; clc; addpath("utils\");

c = 1.0; gamma = 0; N = 2; a = 0.5; M = 0; J = 1000;

% FIX: Tight tolerances are crucial for small parameter values
options = optimset('Display','off', 'TolX',1e-12, 'TolFun',1e-12);

fprintf('--- VISCOELASTIC SCALING CHECKS ---\n');

%% TEST 1: DEPENDENCE ON NU (Fix Frequencies)
fprintf('1. Checking Scaling wrt Nu...\n');
% Use a safe range for nu to avoid machine precision issues
nu_vals = logspace(-4, -2, 20); 
w_min = 1.0; 
w_max = 50.0;
T = pi/w_min; 
dt = pi/w_max;

data_p_corr = zeros(size(nu_vals));
data_q      = zeros(size(nu_vals));
data_rho    = zeros(size(nu_vals));

for i=1:length(nu_vals)
    nu = nu_vals(i);
    
    % Use Theoretical Guess to help optimizer
    q_guess = 0.5 * nu * w_min * w_max / c^3;
    p_guess = 1/c - 0.25 * nu^2 * w_max^2 / c^5;
    
    obj = @(x) obj_Linf(N,T,dt,J,c,gamma,nu,a,M,x(1),x(2));
    [res, fval] = fminsearch(obj, [p_guess, q_guess], options);
    
    % Store absolute drop (1/c - p) and q
    data_p_corr(i) = abs(1/c - res(1)); 
    data_q(i)      = abs(res(2));
    data_rho(i)    = fval;
end

% Robust Fitting (Log-Log Slope)
coeffs_q = polyfit(log(nu_vals), log(data_q), 1);
coeffs_p = polyfit(log(nu_vals), log(data_p_corr), 1);
coeffs_rho = polyfit(log(nu_vals), log(data_rho), 1);

fprintf('   q power wrt nu:   %.2f (Expect 1.0)\n', coeffs_q(1));
fprintf('   p corr power wrt nu:   %.2f (Expect 2.0)\n', coeffs_p(1));
fprintf('   rho power wrt nu: %.2f (Expect 1.0)\n', coeffs_rho(1));

%% TEST 2: DEPENDENCE ON W_MAX (Fix Nu, Fix w_min)
fprintf('2. Checking Scaling wrt w_max...\n');
nu_fix = 1e-3; 
w_max_vals = logspace(1, 2.5, 20); % 10 to ~300
w_min_fix = 1.0; 

data_p_w_corr = zeros(size(w_max_vals));
data_q_w      = zeros(size(w_max_vals));
data_rho_w    = zeros(size(w_max_vals));

for i=1:length(w_max_vals)
    w_max = w_max_vals(i);
    dt_var = pi/w_max;
    
    q_guess = 0.5 * nu_fix * w_min_fix * w_max / c^3;
    p_guess = 1/c - 0.25 * nu_fix^2 * w_max^2 / c^5;
    
    obj = @(x) obj_Linf(N,pi/w_min_fix,dt_var,J,c,gamma,nu_fix,a,M,x(1),x(2));
    [res, fval] = fminsearch(obj, [p_guess, q_guess], options);
    
    data_p_w_corr(i) = abs(1/c - res(1));
    data_q_w(i)      = abs(res(2));
    data_rho_w(i)    = fval;
end

coeffs_q_w = polyfit(log(w_max_vals), log(data_q_w), 1);
coeffs_p_w = polyfit(log(w_max_vals), log(data_p_w_corr), 1);
coeffs_rho_w = polyfit(log(w_max_vals), log(data_rho_w), 1);

fprintf('   q power wrt w_max:   %.2f (Expect 1.0)\n', coeffs_q_w(1));
fprintf('   p corr power wrt w_max:   %.2f (Expect 2.0)\n', coeffs_p_w(1));
fprintf('   rho power wrt w_max: %.2f (Expect 1.0)\n', coeffs_rho_w(1));

%% VISUALIZATION WITH THEORY CURVES
figure('Name','Visco Scaling Validation','Color','w', 'Position', [100 100 1400 500]);

% --- SUBPLOT 1: SCALING WITH NU ---
subplot(1,3,1); 
% Data Points
loglog(nu_vals, data_p_corr, 'bo', 'MarkerFaceColor', 'b', 'DisplayName','p_{corr} Data');
hold on;
loglog(nu_vals, data_q, 'rs', 'MarkerFaceColor', 'r', 'DisplayName','q Data');

% Theory Curves
% p_corr ~ K * nu^2
% q ~ K * nu
theory_p = data_p_corr(end) * (nu_vals / nu_vals(end)).^2; 
theory_q = data_q(end) * (nu_vals / nu_vals(end)).^1;

loglog(nu_vals, theory_p, 'b--', 'LineWidth', 1.5, 'DisplayName','Theory \nu^2');
loglog(nu_vals, theory_q, 'r--', 'LineWidth', 1.5, 'DisplayName','Theory \nu^1');

title('Scaling wrt \nu'); 
xlabel('\nu (Viscosity)'); 
grid on; legend('Location','northwest'); box on;

% --- SUBPLOT 2: SCALING WITH W_MAX ---
subplot(1,3,2); 
% Data Points
loglog(w_max_vals, data_p_w_corr, 'bo', 'MarkerFaceColor', 'b', 'DisplayName','p_{corr} Data');
hold on;
loglog(w_max_vals, data_q_w, 'rs', 'MarkerFaceColor', 'r', 'DisplayName','q Data');

% Theory Curves
% p_corr ~ K * w_max^2
% q ~ K * w_max^1
theory_p_w = data_p_w_corr(end) * (w_max_vals / w_max_vals(end)).^2;
theory_q_w = data_q_w(end) * (w_max_vals / w_max_vals(end)).^1;

loglog(w_max_vals, theory_p_w, 'b--', 'LineWidth', 1.5, 'DisplayName','Theory \omega_{max}^2');
loglog(w_max_vals, theory_q_w, 'r--', 'LineWidth', 1.5, 'DisplayName','Theory \omega_{max}^1');

title('Scaling wrt \omega_{max}'); 
xlabel('\omega_{max}'); 
grid on; legend('Location','northwest'); box on;

% --- SUBPLOT 3: RHO SCALING ---
subplot(1,3,3);
% Data Points
loglog(nu_vals, data_rho, 'kd', 'MarkerFaceColor', 'k', 'DisplayName', '\rho Data');
hold on;
% Theory Curve
% rho ~ K * nu
theory_rho = data_rho(end) * (nu_vals / nu_vals(end)).^1;
loglog(nu_vals, theory_rho, 'g--', 'LineWidth', 2, 'DisplayName', 'Theory \nu^1');

title('Convergence Rate \rho wrt \nu'); 
xlabel('\nu'); 
grid on; legend('Location','northwest'); box on;