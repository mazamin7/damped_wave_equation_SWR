%% CLEAN UP
clear; close all; clc;
addpath("utils\")

%% 1. SIMULATION SETUP
P = get_sim_params();
N = P.N; T = P.T; dt = P.dt; c = P.c; a = P.a; M = P.M; J = P.J;
gamma = 0; % Fixed telegrapher damping

% Focus on very small nu where the Taylor expansion is valid
% (e.g., 1e-5 to 1e-2)
nu_vals = logspace(-7, -4, 100); 
num_nu = length(nu_vals);

%% 2. COMPUTE "REAL" OPTIMAL PARAMETERS (Numerical Optimization)
p_real = zeros(1, num_nu);
q_real = zeros(1, num_nu);

fprintf('Computing numerical optima...\n');
options = optimset('Display','off', 'TolX',1e-9, 'TolFun',1e-9);
x0 = [1/c, 0]; 

for i = 1:num_nu
    nu = nu_vals(i);
    objfun = @(x) obj_Linf(N, T, dt, J, c, gamma, nu, a, M, x(1), x(2));
    [x_opt, ~] = fminsearch(objfun, x0, options);
    
    p_real(i) = x_opt(1);
    q_real(i) = x_opt(2);
    x0 = x_opt; % Update guess
end

%% 3. FIT ANALYTICAL APPROXIMATIONS

% --- Fit q(nu) [Linear] ---
% q = p1*nu + p2  (Allowing small intercept offset)
Coeffs_q = polyfit(nu_vals, q_real, 1); 
q_approx = polyval(Coeffs_q, nu_vals);

% --- Fit p(nu) [Quadratic vs Linear Check] ---
% We try fitting p = A*nu^2 + C (Quadratic decay)
% We do this by fitting a line between p and nu^2
Coeffs_p = polyfit(nu_vals, p_real, 1);
p_approx = polyval(Coeffs_p, nu_vals);

% For comparison, let's look at the Intercept and Slope
p_intercept = Coeffs_p(2);
p_curvature = Coeffs_p(1);

fprintf('Fitted Models:\n');
fprintf('  q(nu) = %.4f * nu + %.4e\n', Coeffs_q(1), Coeffs_q(2));
fprintf('  p(nu) = %.4f + (%.4f) * nu\n', p_intercept, p_curvature);

%% 4. VISUALIZATION
figure('Name','Analytical Fits','Position',[100 100 1200 500], 'Color','w');

% --- Subplot 1: q vs nu ---
subplot(1, 3, 1); hold on; grid on; box on;
plot(nu_vals, q_real, 'ko', 'MarkerFaceColor','k', 'DisplayName', 'Data');
plot(nu_vals, q_approx, 'b-', 'LineWidth', 2, 'DisplayName', 'Linear Fit');
xlabel('\nu'); ylabel('q');
title('q vs \nu (Linear)');
legend('Location','best');
axis square;

% --- Subplot 2: p vs nu (Standard View) ---
subplot(1, 3, 2); hold on; grid on; box on;
plot(nu_vals, p_real, 'ko', 'MarkerFaceColor','k', 'DisplayName', 'Data');
plot(nu_vals, p_approx, 'r-', 'LineWidth', 2, 'DisplayName', 'Linear Fit');
xlabel('\nu'); ylabel('p');
title('p vs \nu (Linear)');
legend('Location','best');
axis square;

% --- Subplot 3: p vs nu^2 (Linearization Check) ---
% If p depends on nu^2, this plot should look like a straight line
subplot(1, 3, 3); hold on; grid on; box on;
plot(nu_vals, p_real, 'ko', 'MarkerFaceColor','k', 'DisplayName', 'Data');
plot(nu_vals, p_approx, 'r-', 'LineWidth', 2, 'DisplayName', 'Fit');
xlabel('\nu^2'); ylabel('p');
title('Linearization Check (p vs \nu)');
legend('Location','best');
axis square;