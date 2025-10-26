%% Compare SWR error surface with contraction factor surfaces for single T
clear all; close all; clc;

addpath("utils\")

%% Parameters
N = 2;
% N = 8;
a = 0.3;
M = 0.1;
% Lx = N*a + M; % automatically determined
T = 5;
c = 1.0;
dh = 0.01;
% dh = 0.001;
dt = 0.01;

% % Viscous damping case
% gamma = 1;
% nu = 0;
% % k = 10;
% k = 5*N;
% % k = 40;

% Viscoelastic damping case
gamma = 0;
nu = 1;
% k = 10;
k = 5*N;

% Parameter ranges
% theta1_range = linspace(0, 1.2, 25);
% theta2_range = linspace(-4, 8, 25);

theta1_range = linspace(0, 1.2, 13);
theta2_range = linspace(-4, 8, 13);

% % for debug
% theta1_range = linspace(0, 1.2, 1);
% theta2_range = linspace(-4, 8, 1);

% theta1_range = linspace(1e-10, 1.5/c, 31);
% theta2_range = linspace(0, 10/c, 21);

% theta1_range = linspace(1e-10, 1.5/c, 16);
% theta2_range = linspace(0, 10/c, 11);

[THETA1, THETA2] = meshgrid(theta1_range, theta2_range);
% J = 500; % frequency axis steps
J = 50;

%% Run single experiment

% Get SWR error surface for this experiment
fprintf('Computing SWR error surface...\n');
error_surface = swr_error_surface_1D(N, a, M, T, c, dh, dt, gamma, nu, k, THETA1, THETA2);

% Get contraction factor surfaces (deterministic, only need once)
fprintf('Computing contraction factor surfaces...\n');
Ly = 0; y_mode = 0; % 1D
[Z_p2, Z_inf] = contraction_surface(N, a, M, Ly, y_mode, T, c, dh, dt, gamma, nu, J, THETA1, THETA2);

%% Extract optimal parameters for the single experiment
[min_val, idx_min] = min(error_surface, [], 'all', 'linear');
[min_i, min_j] = ind2sub(size(error_surface), idx_min);
min_theta1_swr = THETA1(min_i, min_j);
min_theta2_swr = THETA2(min_i, min_j);

%% Find optimal points for theoretical surfaces
% p=2 contraction factor
[~, idx_min_p2] = min(Z_p2(:));
[min_i_p2, min_j_p2] = ind2sub(size(Z_p2), idx_min_p2);
min_theta1_p2 = THETA1(min_i_p2, min_j_p2);
min_theta2_p2 = THETA2(min_i_p2, min_j_p2);
min_val_p2 = Z_p2(min_i_p2, min_j_p2);

% p=∞ contraction factor
[~, idx_min_inf] = min(Z_inf(:));
[min_i_inf, min_j_inf] = ind2sub(size(Z_inf), idx_min_inf);
min_theta1_inf = THETA1(min_i_inf, min_j_inf);
min_theta2_inf = THETA2(min_i_inf, min_j_inf);
min_val_inf = Z_inf(min_i_inf, min_j_inf);

% Initial guess point
theta1_initial = 1/c;
theta2_initial = 0;

%% Create results directory and save figures
results_dir = 'analysis_results';
if ~exist(results_dir, 'dir')
    mkdir(results_dir);
end

%% Create main comparison figure
theta1_range_plot = THETA1(1, :);
theta2_range_plot = THETA2(:, 1)';

figure()
imagesc(theta1_range_plot, theta2_range_plot, log10(error_surface));
axis xy;
hold on;
plot(theta1_initial, theta2_initial, 'ks', 'MarkerSize', 10, 'MarkerFaceColor', 'k', 'DisplayName', 'Initial guess');
plot(min_theta1_swr, min_theta2_swr, 'rs', 'MarkerSize', 10, 'MarkerFaceColor', 'r', 'LineWidth', 2, 'DisplayName', 'SWR error opt.');
plot(min_theta1_p2, min_theta2_p2, 'm*', 'MarkerSize', 10, 'MarkerFaceColor', 'm', 'DisplayName', 'Spectral opt. L2');
plot(min_theta1_inf, min_theta2_inf, 'g^', 'MarkerSize', 10, 'MarkerFaceColor', 'g', 'DisplayName', 'Spectral opt. L∞');
colorbar;
xlabel('p', 'FontSize', 16);
ylabel('q', 'FontSize', 16);
title('Log-scale Error Surface');
legend('show', 'Location', 'NorthEast', 'FontSize', 14);
colormap('parula');
% clim([-7,3])
% Make tick labels bigger
set(gca, 'FontSize', 18);
saveas(gcf, fullfile(results_dir, 'surface_comparison.png'));

%%
figure()
imagesc(theta1_range_plot, theta2_range_plot, Z_p2);
axis xy;
hold on;
plot(theta1_initial, theta2_initial, 'ks', 'MarkerSize', 10, 'MarkerFaceColor', 'k', 'DisplayName', 'Initial guess');
plot(min_theta1_swr, min_theta2_swr, 'rs', 'MarkerSize', 10, 'MarkerFaceColor', 'r', 'LineWidth', 2, 'DisplayName', 'SWR error opt.');
plot(min_theta1_p2, min_theta2_p2, 'm*', 'MarkerSize', 10, 'MarkerFaceColor', 'm', 'DisplayName', 'Spectral opt. L2');
plot(min_theta1_inf, min_theta2_inf, 'g^', 'MarkerSize', 10, 'MarkerFaceColor', 'g', 'DisplayName', 'Spectral opt. L∞');
colorbar;
xlabel('p', 'FontSize', 16);
ylabel('q', 'FontSize', 16);
title('L2 global contraction factor');
legend('show', 'Location', 'NorthEast', 'FontSize', 14);
colormap('parula');
% Make tick labels bigger
set(gca, 'FontSize', 18);

%%
figure()
imagesc(theta1_range_plot, theta2_range_plot, Z_inf);
axis xy;
hold on;
plot(theta1_initial, theta2_initial, 'ks', 'MarkerSize', 10, 'MarkerFaceColor', 'k', 'DisplayName', 'Initial guess');
plot(min_theta1_swr, min_theta2_swr, 'rs', 'MarkerSize', 10, 'MarkerFaceColor', 'r', 'LineWidth', 2, 'DisplayName', 'SWR error opt.');
plot(min_theta1_p2, min_theta2_p2, 'm*', 'MarkerSize', 10, 'MarkerFaceColor', 'm', 'DisplayName', 'Spectral opt. L2');
plot(min_theta1_inf, min_theta2_inf, 'g^', 'MarkerSize', 10, 'MarkerFaceColor', 'g', 'DisplayName', 'Spectral opt. L∞');
colorbar;
xlabel('p', 'FontSize', 16);
ylabel('q', 'FontSize', 16);
title('Linf global contraction factor');
legend('show', 'Location', 'NorthEast', 'FontSize', 14);
colormap('parula');
% Make tick labels bigger
set(gca, 'FontSize', 18);

%% Display summary
fprintf('\n=== RESULTS SUMMARY ===\n');

fprintf('\nComparison of Methods:\n');
fprintf('Method         theta1    theta2\n');
fprintf('----------------------------------------------------------\n');
fprintf('SWR error opt. %8.4f     %8.4f\n', min_theta1_swr, min_theta2_swr);
fprintf('L2             %8.4f     %8.4f\n', min_theta1_p2, min_theta2_p2);
fprintf('Linf           %8.4f     %8.4f\n', min_theta1_inf, min_theta2_inf);
