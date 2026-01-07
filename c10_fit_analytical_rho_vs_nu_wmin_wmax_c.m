%% CLEAN UP
clear; close all; clc;
addpath("utils\")

%% 1. CONFIGURATION
% Fixed geometry params
gamma = 0; N = 2; a = 0.3; M = 0; J = 1000;

%% 2. GENERATE RANDOM DATA
NumPoints = 40;
fprintf('Generating %d test points to analyze Contraction Factor rho...\n', NumPoints);

% Parameter Ranges
range_nu    = [1e-6, 1e-4];
range_c     = [0.5, 3.0];
range_wmax  = [300, 3000]; % High freq
range_wmin  = [0.1, 2.0];  % Low freq

data_term = zeros(NumPoints, 1); % The theoretical term: nu * (w_max - w_min) / c^2
data_rho  = zeros(NumPoints, 1); % The actual optimal rho

options = optimset('Display','off', 'TolX',1e-9, 'TolFun',1e-9);

for i = 1:NumPoints
    % 1. Random Parameters
    nu = 10^(log10(range_nu(1)) + rand()*diff(log10(range_nu)));
    c  = range_c(1) + rand()*diff(range_c);
    w_max = 10^(log10(range_wmax(1)) + rand()*diff(log10(range_wmax)));
    w_min = range_wmin(1) + rand()*diff(range_wmin);
    
    dt = pi / w_max;
    T  = pi / w_min;
    
    % 2. Calculate Theoretical Scaling Term
    % Hypothesis: rho ~ K * nu * (w_max - w_min) / c^2
    data_term(i) = (nu * (w_max - w_min)) / c^2;
    
    % 3. Run Optimization
    objfun = @(x) obj_Linf(N, T, dt, J, c, gamma, nu, a, M, x(1), x(2));
    
    % Use our derived analytical formulas for the guess! (Faster convergence)
    q_guess = 0.5 * (nu/c^3) * w_min * w_max;
    p_guess = 1/c - 0.25 * (nu^2/c^5) * w_max^2;
    x0 = [p_guess, q_guess];
    
    [~, rho_val] = fminsearch(objfun, x0, options);
    
    data_rho(i) = rho_val;
    
    if mod(i,10)==0, fprintf('.'); end
end
fprintf('\n');

%% 3. FIT THE CONSTANT
% Model: rho = K * term
K_rho = mean(data_rho ./ data_term);

fprintf('\n=== CONTRACTION RATE FORMULA ===\n');
fprintf('Target Model: rho = K * [ nu * (w_max - w_min) / c^2 ]\n');
fprintf('Fitted Constant K: %.4f\n', K_rho);
fprintf('Theoretical Prediction: 0.2500\n');
fprintf('Error: %.2f%%\n', abs(K_rho - 0.25)/0.25 * 100);

%% 4. VISUALIZATION
figure('Name','Rho Scaling','Position',[200 200 600 500], 'Color','w');
plot(data_term, data_rho, 'bo', 'MarkerFaceColor','b', 'MarkerSize', 6);
hold on;

% Plot K=0.25 Line
x_lin = linspace(min(data_term), max(data_term), 100);
plot(x_lin, 0.25 * x_lin, 'r--', 'LineWidth', 2);

xlabel('Theoretical Term: $\frac{\nu (\omega_{\max} - \omega_{\min})}{c^2}$', 'Interpreter','latex', 'FontSize', 14);
ylabel('Optimal Contraction $\rho_{opt}$', 'Interpreter','latex', 'FontSize', 14);
title(['Scaling of Contraction Factor (K \approx ' num2str(K_rho, '%.3f') ')']);
legend('Numerical Data', 'Theory (K=0.25)', 'Location','best');
grid on; axis square;