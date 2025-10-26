clear all; close all; clc;

addpath("utils\")

%% Common parameters
N = 2;  % Number of subdomains
% N = 8;
c = 1.0;
M = 0.1;
a = 0.3;
b = a + M;
T = 5;

ky = 0; % 1D
% ky = 5;

dh = 0.01;
dt = 0.01;
% dh = 0.001;
% dt = 0.001;
% J = 500; % number of frequency steps from 2*pi/T to pi/dt
J = 50;

tol = 1e-6;
% tol = 1e-9;

%% Define parameter sets for different cases (updated as per advisor)
% Cases: Varying gamma with different fixed nu values
gamma_values = linspace(0,6,15);  % Varying gamma from 0 to 3
nu_fixed_values = 2*[0, 0.05, 0.1, 0.5, 1, 3];  % Fixed nu values (replaced 5 with 3)
num_gamma_cases = length(nu_fixed_values);

% Cases: Varying nu with different fixed gamma values  
nu_values = linspace(0,6,15);  % Varying nu from 0 to 3
gamma_fixed_values = 2*[0, 0.05, 0.1, 0.5, 1, 3];  % Fixed gamma values (replaced 5 with 3)
num_nu_cases = length(gamma_fixed_values);

%% Initialize storage arrays (only L∞)
% For gamma-varying cases
optimal_pinf_gamma = cell(num_gamma_cases, 1);
contraction_pinf_gamma = cell(num_gamma_cases, 1);

% For nu-varying cases
optimal_pinf_nu = cell(num_nu_cases, 1);
contraction_pinf_nu = cell(num_nu_cases, 1);

%% Process gamma-varying cases using a loop (only L∞)
fprintf('=== Processing gamma-Varying Cases (L∞ only) ===\n');

for case_idx = 1:num_gamma_cases
    nu_fixed = nu_fixed_values(case_idx);
    case_info = sprintf('Varying gamma (nu = %.2f)', nu_fixed);
    fprintf('\n=== Case %d: %s ===\n', case_idx, case_info);
    
    % Initialize arrays for this case
    optimal_pinf_gamma{case_idx} = zeros(length(gamma_values), 2);
    contraction_pinf_gamma{case_idx} = zeros(length(gamma_values), 1);
    
    % Initial guesses for this case
    x0_pinf = [1/c, 0];
    
    % Compute optimal parameters (only L∞)
    for i = 1:length(gamma_values)
        gamma = gamma_values(i);
        fprintf('Computing L∞ optimal parameters for gamma = %.3f (%d/%d)\n', gamma, i, length(gamma_values));
        
        % L∞ optimization only
        objfun_pinf = @(x) obj_Linf(N, T, dt, J, c, gamma, nu_fixed, a, M, x(1), x(2), ky);
        options = optimset('Display', 'off', 'TolX', tol, 'TolFun', tol);
        [x_opt_pinf, fval_pinf] = fminsearch(objfun_pinf, x0_pinf, options);
        optimal_pinf_gamma{case_idx}(i,:) = x_opt_pinf;
        contraction_pinf_gamma{case_idx}(i) = fval_pinf;
        x0_pinf = x_opt_pinf;
    end
end

%% Process nu-varying cases using a loop (only L∞)
fprintf('\n=== Processing nu-Varying Cases (L∞ only) ===\n');

for case_idx = 1:num_nu_cases
    gamma_fixed = gamma_fixed_values(case_idx);
    case_info = sprintf('Varying nu (gamma = %.2f)', gamma_fixed);
    fprintf('\n=== Case %d: %s ===\n', case_idx + num_gamma_cases, case_info);
    
    % Initialize arrays for this case
    optimal_pinf_nu{case_idx} = zeros(length(nu_values), 2);
    contraction_pinf_nu{case_idx} = zeros(length(nu_values), 1);
    
    % Initial guesses for this case
    x0_pinf = [1/c, 0];
    
    % Compute optimal parameters (only L∞)
    for i = 1:length(nu_values)
        nu = nu_values(i);
        fprintf('Computing L∞ optimal parameters for nu = %.3f (%d/%d)\n', nu, i, length(nu_values));
        
        % L∞ optimization only
        objfun_pinf = @(x) obj_Linf(N, T, dt, J, c, gamma_fixed, nu, a, M, x(1), x(2), ky);
        options = optimset('Display', 'off', 'TolX', tol, 'TolFun', tol);
        [x_opt_pinf, fval_pinf] = fminsearch(objfun_pinf, x0_pinf, options);
        optimal_pinf_nu{case_idx}(i,:) = x_opt_pinf;
        contraction_pinf_nu{case_idx}(i) = fval_pinf;
        x0_pinf = x_opt_pinf;
    end
end

%% Figure 1: gamma cases (Viscous damping)
figure('Position', [100, 100, 1200, 500]);

% Define colors for better distinction with more cases
colors = lines(max(num_gamma_cases, num_nu_cases));

% Subplot 1: gamma-varying cases - theta1-theta2 plane
subplot(1,2,1);
hold on;
for case_idx = 1:num_gamma_cases
    theta1_vals = optimal_pinf_gamma{case_idx}(:,1);
    theta2_vals = optimal_pinf_gamma{case_idx}(:,2);
    
    % Plot the trajectory
    plot(theta1_vals, theta2_vals, '-', 'LineWidth', 2, ...
         'DisplayName', sprintf('\\nu=%.2f', nu_fixed_values(case_idx)), ...
         'Color', colors(case_idx,:));
    
    % Mark the starting point (gamma = 0) - circle
    plot(theta1_vals(1), theta2_vals(1), 'o', 'MarkerSize', 8, ...
         'MarkerFaceColor', colors(case_idx,:), 'MarkerEdgeColor', 'k', 'HandleVisibility', 'off');
    
    % Mark the ending point (gamma = 3) - square
    plot(theta1_vals(end), theta2_vals(end), 's', 'MarkerSize', 8, ...
         'MarkerFaceColor', colors(case_idx,:), 'MarkerEdgeColor', 'k', 'HandleVisibility', 'off');
end

% Get the position of the current axes in normalized figure coordinates
ax_pos = get(gca, 'Position');

% Common tail position (center top of the plot)
tail_x = ax_pos(1) + ax_pos(3)*0.5 - 0.056;  % Center of x-axis
tail_y = ax_pos(2) + ax_pos(4)*0.85 - 0.068; % Near top of y-axis

% Add arrows using annotation (normalized coordinates) with common tail
% Vertical arrow for gamma increasing (pointing up)
annotation('arrow', [tail_x, tail_x], ...
           [tail_y, tail_y + ax_pos(4)*0.08], ...
           'LineWidth', 3, 'HeadWidth', 15, 'HeadLength', 12, 'Color', 'k');

% Horizontal arrow for nu increasing (pointing left)
annotation('arrow', [tail_x, tail_x - ax_pos(3)*0.08], ...
           [tail_y, tail_y], ...
           'LineWidth', 3, 'HeadWidth', 15, 'HeadLength', 12, 'Color', 'k');

% Add labels for the arrows
text(0.46, 7.15, '\gamma increasing', 'FontSize', 16, 'HorizontalAlignment', 'center');
text(0.46, 5.6, '\nu increasing', 'FontSize', 16, 'HorizontalAlignment', 'center');

xlabel('p', 'FontSize', 16);
ylabel('q', 'FontSize', 16);
%title('L∞ Optimal Parameters: Varying \gamma', 'FontSize', 12);

% Create custom legend with marker explanations
legend_handles = gobjects(num_gamma_cases + 2, 1);
legend_labels = cell(num_gamma_cases + 2, 1);

for case_idx = 1:num_gamma_cases
    legend_handles(case_idx) = plot(NaN, NaN, '-', 'LineWidth', 2, 'Color', colors(case_idx,:));
    legend_labels{case_idx} = sprintf('\\nu=%.2f', nu_fixed_values(case_idx));
end

% Add marker explanations to legend
legend_handles(num_gamma_cases + 1) = plot(NaN, NaN, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'k');
legend_labels{num_gamma_cases + 1} = 'Start: \gamma = 0';

legend_handles(num_gamma_cases + 2) = plot(NaN, NaN, 'ks', 'MarkerSize', 8, 'MarkerFaceColor', 'k'); 
legend_labels{num_gamma_cases + 2} = 'End: \gamma = 6';

legend(legend_handles, legend_labels, 'Location', 'NorthEast', 'FontSize', 14);
grid on;
xlim([0, 1.5]);
ylim([-4, 8]);
% Make tick labels bigger
set(gca, 'FontSize', 18);

% Subplot 2: gamma-varying cases - contraction factor vs gamma
subplot(1,2,2);
hold on;
for case_idx = 1:num_gamma_cases
    semilogy(gamma_values, contraction_pinf_gamma{case_idx}, '-', 'LineWidth', 2, ...
         'DisplayName', sprintf('\\nu=%.2f', nu_fixed_values(case_idx)), ...
         'Color', colors(case_idx,:));
end
xlabel('\gamma', 'FontSize', 16);
ylabel('Global Contraction Factor', 'FontSize', 16);
%title('Contraction Factor vs \gamma (L∞)', 'FontSize', 12);
legend('Location', 'NorthEast', 'FontSize', 14);
grid on;
ylim([1e-3, 1]);
% Make tick labels bigger
set(gca, 'FontSize', 18);

%% Figure 2: nu cases (Viscoelastic damping)
figure('Position', [100, 100, 1200, 500]);

% Subplot 1: nu-varying cases - theta1-theta2 plane
subplot(1,2,1);
hold on;
for case_idx = 1:num_nu_cases
    theta1_vals = optimal_pinf_nu{case_idx}(:,1);
    theta2_vals = optimal_pinf_nu{case_idx}(:,2);
    
    % Plot the trajectory
    plot(theta1_vals, theta2_vals, '-', 'LineWidth', 2, ...
         'DisplayName', sprintf('\\gamma=%.2f', gamma_fixed_values(case_idx)), ...
         'Color', colors(case_idx,:));
    
    % Mark the starting point (nu = 0) - circle
    plot(theta1_vals(1), theta2_vals(1), 'o', 'MarkerSize', 8, ...
         'MarkerFaceColor', colors(case_idx,:), 'MarkerEdgeColor', 'k', 'HandleVisibility', 'off');
    
    % Mark the ending point (nu = 3) - square
    plot(theta1_vals(end), theta2_vals(end), 's', 'MarkerSize', 8, ...
         'MarkerFaceColor', colors(case_idx,:), 'MarkerEdgeColor', 'k', 'HandleVisibility', 'off');
end

% Get the position of the current axes in normalized figure coordinates
ax_pos = get(gca, 'Position');

% Common tail position (center top of the plot)
tail_x = ax_pos(1) + ax_pos(3)*0.5 - 0.056;  % Center of x-axis
tail_y = ax_pos(2) + ax_pos(4)*0.85 - 0.068; % Near top of y-axis

% Add arrows using annotation (normalized coordinates) with common tail
% Vertical arrow for gamma increasing (pointing up)
annotation('arrow', [tail_x, tail_x], ...
           [tail_y, tail_y + ax_pos(4)*0.08], ...
           'LineWidth', 3, 'HeadWidth', 15, 'HeadLength', 12, 'Color', 'k');

% Horizontal arrow for nu increasing (pointing left)
annotation('arrow', [tail_x, tail_x - ax_pos(3)*0.08], ...
           [tail_y, tail_y], ...
           'LineWidth', 3, 'HeadWidth', 15, 'HeadLength', 12, 'Color', 'k');

% Add labels for the arrows
text(0.46, 7.15, '\gamma increasing', 'FontSize', 16, 'HorizontalAlignment', 'center');
text(0.46, 5.6, '\nu increasing', 'FontSize', 16, 'HorizontalAlignment', 'center');

xlabel('p', 'FontSize', 14);
ylabel('q', 'FontSize', 14);
%title('L∞ Optimal Parameters: Varying \nu', 'FontSize', 12);

% Create custom legend with marker explanations
legend_handles = gobjects(num_nu_cases + 2, 1);
legend_labels = cell(num_nu_cases + 2, 1);

for case_idx = 1:num_nu_cases
    legend_handles(case_idx) = plot(NaN, NaN, '-', 'LineWidth', 2, 'Color', colors(case_idx,:));
    legend_labels{case_idx} = sprintf('\\gamma=%.2f', gamma_fixed_values(case_idx));
end

% Add marker explanations to legend
legend_handles(num_nu_cases + 1) = plot(NaN, NaN, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'k');
legend_labels{num_nu_cases + 1} = 'Start: \nu = 0';

legend_handles(num_nu_cases + 2) = plot(NaN, NaN, 'ks', 'MarkerSize', 8, 'MarkerFaceColor', 'k');
legend_labels{num_nu_cases + 2} = 'End: \nu = 6';

legend(legend_handles, legend_labels, 'Location', 'NorthEast', 'FontSize', 14);
grid on;
xlim([0, 1.5]);
ylim([-4, 8]);
% Make tick labels bigger
set(gca, 'FontSize', 18);

% Subplot 2: nu-varying cases - contraction factor vs nu
subplot(1,2,2);
hold on;
for case_idx = 1:num_nu_cases
    semilogy(nu_values, contraction_pinf_nu{case_idx}, '-', 'LineWidth', 2, ...
         'DisplayName', sprintf('\\gamma=%.2f', gamma_fixed_values(case_idx)), ...
         'Color', colors(case_idx,:));
end
xlabel('\nu', 'FontSize', 14);
ylabel('Global Contraction Factor', 'FontSize', 14);
%title('Contraction Factor vs \nu (L∞)', 'FontSize', 12);
legend('Location', 'NorthEast', 'FontSize', 14);
grid on;
ylim([1e-3, 1]);
% Make tick labels bigger
set(gca, 'FontSize', 18);

%% Helper functions (unchanged)

function final_error = swr_final_error_for_alpha_explicit(u0, v0, Lx, N, a, M, T, c, dh, gamma, nu, theta1, theta2, k, u_init)
    % SWR final error function using explicit parameters
    u_ref = run_fdtd_reference(u0, v0, Lx, T, c, dh, dh/c, gamma, nu);

    [~, final_error, ~] = run_swr(u0, v0, N, a, M, T, c, dh, gamma, nu, theta1, theta2, k, u_init, u_ref);
end
