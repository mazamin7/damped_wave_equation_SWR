%% CLEAN UP
clear; close all; clc;
addpath("utils\")

%% 1. CONFIGURATION
P_base = get_sim_params();
c = P_base.c; 
gamma = 0; % Fixed

% Define the Parameter Space to explore
% We create a "cloud" of random test points
NumPoints = 40; 

% Random ranges (Log scale for nu and dt)
range_nu = [1e-6, 1e-4];
range_dt = [0.001, 0.01];
range_T  = [2, 20];

% Preallocate Data Arrays
data_nu = zeros(NumPoints, 1);
data_dt = zeros(NumPoints, 1);
data_T  = zeros(NumPoints, 1);
data_q  = zeros(NumPoints, 1);
data_p  = zeros(NumPoints, 1);

%% 2. GENERATE DATA (Randomized Sweep)
fprintf('Generating %d sample points across parameter space...\n', NumPoints);
options = optimset('Display','off', 'TolX',1e-8, 'TolFun',1e-8);

for i = 1:NumPoints
    % 1. Pick random parameters
    nu_val = 10^( log10(range_nu(1)) + rand() * (log10(range_nu(2)) - log10(range_nu(1))) );
    dt_val = 10^( log10(range_dt(1)) + rand() * (log10(range_dt(2)) - log10(range_dt(1))) );
    T_val  = range_T(1) + rand() * (range_T(2) - range_T(1));
    
    % Store inputs
    data_nu(i) = nu_val;
    data_dt(i) = dt_val;
    data_T(i)  = T_val;
    
    % 2. Run Optimization
    % Note: N is derived from Length L=1 usually, but here we assume fixed spatial grid
    % or we scale N? Let's assume fixed spatial parameters from P_base
    N = P_base.N; a = P_base.a; M = P_base.M; J = P_base.J;
    
    objfun = @(x) obj_Linf(N, T_val, dt_val, J, c, gamma, nu_val, a, M, x(1), x(2));
    
    % Use simple warm start
    x0 = [1/c, 0];
    [x_opt, ~] = fminsearch(objfun, x0, options);
    
    % Store outputs
    data_p(i) = x_opt(1);
    data_q(i) = x_opt(2);
    
    fprintf('.');
    if mod(i,20)==0, fprintf('\n'); end
end
fprintf('Done.\n');

%% 3. FIT GLOBAL CONSTANTS
% Model q: q = Cq * (nu / (dt * T))
% Linear regression for Cq
X_q = data_nu ./ (data_dt .* data_T);
C_q = (X_q' * X_q) \ (X_q' * data_q);

% Model p: p = 1/c - Cp * (nu / dt)^2
% Linear regression for Cp: Cp * (nu/dt)^2 = (1/c - p)
X_p = (data_nu ./ data_dt).^2;
Y_p = (1/c) - data_p;
C_p = (X_p' * X_p) \ (X_p' * Y_p);

fprintf('\n=== GLOBAL FITTED FORMULAS ===\n');
fprintf('q(nu, dt, T) = %.4f * [ nu / (dt * T) ]\n', C_q);
fprintf('p(nu, dt)    = 1/c - %.4f * [ nu / dt ]^2\n', C_p);

%% 4. VALIDATION PLOTS
figure('Name','Global Formula Validation','Position',[100 100 1000 500], 'Color','w');

% --- VALIDATION FOR Q ---
subplot(1, 2, 1); hold on; grid on; box on;
q_pred = C_q * X_q;
plot(data_q, q_pred, 'bo', 'MarkerFaceColor', 'b', 'MarkerSize', 6);
plot([min(data_q) max(data_q)], [min(data_q) max(data_q)], 'k--', 'LineWidth', 2);
xlabel('Numerical Optimal q');
ylabel('Formula Predicted q');
title({'Validation of q Formula', 'q \approx C_q \cdot \nu / (\Delta t \cdot T)'});
axis square;

% --- VALIDATION FOR P ---
subplot(1, 2, 2); hold on; grid on; box on;
p_pred = 1/c - C_p * X_p;
plot(data_p, p_pred, 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 6);
plot([min(data_p) max(data_p)], [min(data_p) max(data_p)], 'k--', 'LineWidth', 2);
xlabel('Numerical Optimal p');
ylabel('Formula Predicted p');
title({'Validation of p Formula', 'p \approx 1/c - C_p \cdot (\nu / \Delta t)^2'});
axis square;