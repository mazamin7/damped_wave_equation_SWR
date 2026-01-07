%% CLEAN UP
clear; close all; clc;
addpath("utils\")

%% 1. SETUP
P = get_sim_params();
nu_fixed = 1e-5; % Fix a small nu where we expect parabolic behavior
gamma = 0;
M = P.M;
c = P.c;

% High resolution sweep of 'a' to capture oscillations
a_sweep = 0.2 : 0.01 : 1.8; 
p_results = zeros(size(a_sweep));

fprintf('Scanning interface positions...\n');
options = optimset('Display','off', 'TolX',1e-9, 'TolFun',1e-9);
x0 = [1/c, 0];

%% 2. SWEEP
for i = 1:length(a_sweep)
    a_val = a_sweep(i);
    
    % Objective function
    objfun = @(x) obj_Linf(P.N, P.T, P.dt, P.J, c, gamma, nu_fixed, a_val, M, x(1), x(2));
    
    [x_opt, ~] = fminsearch(objfun, x0, options);
    p_results(i) = x_opt(1);
    
    % Warm start next step
    x0 = x_opt; 
end

%% 3. VISUALIZATION
figure('Name','Resonance Check','Position',[100 100 800 500], 'Color','w');
hold on; grid on; box on;

% Plot raw p values
plot(a_sweep, p_results, 'r-o', 'LineWidth', 1.5, 'MarkerSize', 4);

% Draw reference line at 1/c
yline(1/c, 'k--', 'Wave Speed limit (1/c)', 'LineWidth', 2);

xlabel('Interface Position a');
ylabel('Optimal p');
title(['Oscillation of p with Geometry (Fixed \nu = ' num2str(nu_fixed) ')']);

% Add context
subtitle('Oscillations suggest grid alignment/resonance effects');