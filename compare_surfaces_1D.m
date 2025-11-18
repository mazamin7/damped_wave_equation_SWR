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
% dh = 0.01;
% dt = 0.01;
dh = 0.002;
dt = 0.002;
% dh = 0.001;
% dt = 0.001;

% % Viscous damping case
% % gamma = 0.1;
% % gamma = 1;
% % gamma = 5;
% gamma = 10;
% nu = 0;
% % k = 10;
% k = 5*N;
% % k = 40;

% Viscoelastic damping case
gamma = 0;
nu = 0.01;
% nu = 0.05;
% nu = 0.5;
% nu = 1;
% k = 10;
k = 5*N;

% Parameter ranges
% theta1_range = linspace(0, 1.2, 25);
% theta2_range = linspace(-4, 8, 25);

theta1_range = linspace(0, 1.2, 13);
theta2_range = linspace(-4, 8, 13);

[THETA1, THETA2] = meshgrid(theta1_range, theta2_range);

theta1_min = min(theta1_range);
theta1_max = max(theta1_range);
theta2_min = min(theta2_range);
theta2_max = max(theta2_range);

clip = @(v, vmin, vmax) max(vmin, min(vmax, v));

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

%% Extract optimal parameters for the single experiment (SWR)
[min_val, idx_min] = min(error_surface, [], 'all', 'linear');
[min_i, min_j] = ind2sub(size(error_surface), idx_min);
min_theta1_swr = THETA1(min_i, min_j);
min_theta2_swr = THETA2(min_i, min_j);

%% Find optimal points for theoretical surfaces via fminsearch (ky = 0)
ky = 0;

% Initial guess point (also used in plots)
theta1_initial = 1/c;
theta2_initial = 0;
x0 = [theta1_initial, theta2_initial];

% fminsearch options
base_options = optimset('Display','off', ...
                        'TolX',1e-8, ...
                        'TolFun',1e-8);

% p=2 contraction factor minimization
objfun_p2 = @(x) obj_L2( N, T, dt, J, c, gamma, nu, a, M, x(1), x(2), ky );
[x_opt_p2, min_val_p2] = fminsearch(objfun_p2, x0, base_options);
min_theta1_p2 = clip(x_opt_p2(1), theta1_min, theta1_max);
min_theta2_p2 = clip(x_opt_p2(2), theta2_min, theta2_max);

% p=∞ contraction factor minimization
objfun_inf = @(x) obj_Linf( N, T, dt, J, c, gamma, nu, a, M, x(1), x(2), ky );
[x_opt_inf, min_val_inf] = fminsearch(objfun_inf, x0, base_options);
min_theta1_inf = clip(x_opt_inf(1), theta1_min, theta1_max);
min_theta2_inf = clip(x_opt_inf(2), theta2_min, theta2_max);

%% Create results directory and save figures
results_dir = 'analysis_results';
if ~exist(results_dir, 'dir')
    mkdir(results_dir);
end

%% Create main comparison figure
theta1_range_plot = THETA1(1, :);
theta2_range_plot = THETA2(:, 1)';

figure()
contourf(theta1_range_plot, theta2_range_plot, log10(error_surface), 50, 'LineStyle', 'none', 'HandleVisibility', 'off');
axis xy;
hold on;
plot(theta1_initial, theta2_initial, 'ks', 'MarkerSize', 10, 'MarkerFaceColor', 'k', 'DisplayName', 'Initial guess');
plot(min_theta1_swr, min_theta2_swr, 'rs', 'MarkerSize', 10, 'MarkerFaceColor', 'r', 'DisplayName', 'SWR error opt.');
plot(min_theta1_p2,  min_theta2_p2,  'm*', 'MarkerSize', 10, 'MarkerFaceColor', 'm', 'DisplayName', 'Spectral opt. L2');
plot(min_theta1_inf, min_theta2_inf, 'g^', 'MarkerSize', 10, 'MarkerFaceColor', 'g', 'DisplayName', 'Spectral opt. L∞');
colorbar;
xlabel('p', 'FontSize', 16);
ylabel('q', 'FontSize', 16);
% title('Log-scale Error Surface');
legend('show', 'Location', 'SouthWest', 'FontSize', 14);
colormap('parula');
if gamma == 0
    clim([min(log10(error_surface),[],'all'),3])
elseif nu == 0
    clim([min(log10(error_surface),[],'all'),8])
end
% Make tick labels bigger
set(gca, 'FontSize', 18);
saveas(gcf, fullfile(results_dir, 'surface_comparison.png'));

% %%
% figure()
% imagesc(theta1_range_plot, theta2_range_plot, Z_p2);
% axis xy;
% hold on;
% plot(theta1_initial, theta2_initial, 'ks', 'MarkerSize', 10, 'MarkerFaceColor', 'k', 'DisplayName', 'Initial guess');
% plot(min_theta1_swr, min_theta2_swr, 'rs', 'MarkerSize', 10, 'MarkerFaceColor', 'r', 'DisplayName', 'SWR error opt.');
% plot(min_theta1_p2,  min_theta2_p2,  'm*', 'MarkerSize', 10, 'MarkerFaceColor', 'm', 'DisplayName', 'Spectral opt. L2');
% plot(min_theta1_inf, min_theta2_inf, 'g^', 'MarkerSize', 10, 'MarkerFaceColor', 'g', 'DisplayName', 'Spectral opt. L∞');
% colorbar;
% xlabel('p', 'FontSize', 16);
% ylabel('q', 'FontSize', 16);
% title('L2 global contraction factor');
% legend('show', 'Location', 'NorthEast', 'FontSize', 14);
% colormap('parula');
% clim([0,1])
% % Make tick labels bigger
% set(gca, 'FontSize', 18);
% 
% %%
% figure()
% imagesc(theta1_range_plot, theta2_range_plot, Z_inf);
% axis xy;
% hold on;
% plot(theta1_initial, theta2_initial, 'ks', 'MarkerSize', 10, 'MarkerFaceColor', 'k', 'DisplayName', 'Initial guess');
% plot(min_theta1_swr, min_theta2_swr, 'rs', 'MarkerSize', 10, 'MarkerFaceColor', 'r', 'DisplayName', 'SWR error opt.');
% plot(min_theta1_p2,  min_theta2_p2,  'm*', 'MarkerSize', 10, 'MarkerFaceColor', 'm', 'DisplayName', 'Spectral opt. L2');
% plot(min_theta1_inf, min_theta2_inf, 'g^', 'MarkerSize', 10, 'MarkerFaceColor', 'g', 'DisplayName', 'Spectral opt. L∞');
% colorbar;
% xlabel('p', 'FontSize', 16);
% ylabel('q', 'FontSize', 16);
% title('Linf global contraction factor');
% legend('show', 'Location', 'NorthEast', 'FontSize', 14);
% colormap('parula');
% clim([0,1])
% % Make tick labels bigger
% set(gca, 'FontSize', 18);

%% Display summary
fprintf('\n=== RESULTS SUMMARY ===\n');

fprintf('\nComparison of Methods:\n');
fprintf('Method         theta1    theta2\n');
fprintf('----------------------------------------------------------\n');
fprintf('SWR error opt. %8.4f     %8.4f\n', min_theta1_swr, min_theta2_swr);
fprintf('L2             %8.4f     %8.4f\n', min_theta1_p2,  min_theta2_p2);
fprintf('Linf           %8.4f     %8.4f\n', min_theta1_inf, min_theta2_inf);
