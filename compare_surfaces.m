%% Compare SWR error surface with contraction factor surfaces for single T
clear all; close all; clc;
addpath("utils\")

%% VISUALIZATION SETTINGS (POSTER STYLE)
FS_AXIS  = 24;  % Font size for axis ticks
FS_LABEL = 30;  % Font size for X/Y labels
FS_TITLE = 28;  % Font size for titles
MS_MARK  = 15;  % Marker size (Increased for visibility)
LW_AXIS  = 2.0; % Line width for the axis box

%% Parameters
P = get_sim_params();
N  = P.N;
a  = P.a;
M  = P.M;
b  = P.b;
Lx = P.Lx;
c  = P.c;
% T  = P.T;
gamma = P.gamma; % ignore
nu = P.nu; % ignore
dh = P.dh;
dt = P.dt;
J  = P.J;

dist_to_bound = (1/3) * Lx;
T = 0.9 * (dist_to_bound / c); % Set T slightly less than travel time to boundary

%%
% % Viscous damping case
% gamma = 4;
% % gamma = 8;
% % gamma = 10;
% % gamma = 12;
% nu = 0;
% % k = 10;
% k = 5*N;
% % k = 40;

% Viscoelastic damping case
gamma = 0;
nu = 0.001;
% nu = 0.01;
% nu = 0.05;
% nu = 0.1;
% k = 10;
k = 5*N;
% k = 20;

% Parameter ranges
theta1_range = linspace(0, 2, 13);
theta2_range = linspace(-4, 8, 13);
[THETA1, THETA2] = meshgrid(theta1_range, theta2_range);
theta1_min = min(theta1_range);
theta1_max = max(theta1_range);
theta2_min = min(theta2_range);
theta2_max = max(theta2_range);

clip = @(v, vmin, vmax) max(vmin, min(vmax, v));

%% Get SWR error surface
% 1. Setup necessary fields for the objective function
Nx = round(Lx/dh) + 1;
Nt = round(T/dt) + 1;
x_grid = linspace(0,Lx,Nx);

% Initial conditions components
gaussian = @(r,mu,sigma) 1/(2*pi*sigma^2) * exp(-(r-mu).^2/(2*sigma^2));
    
% Create normalized components as function handles
gaussian_normalized = @(x) gaussian(x, 2/3*Lx, Lx/200);
    
% Create the final initial condition function
u0 = @(x) gaussian_normalized(x);
v0 = @(x) 0.*x;

figure('Color', 'w')
plot(x_grid,u0(x_grid), 'LineWidth', 2)
title('Initial Condition', 'FontSize', FS_TITLE, 'FontWeight', 'bold');
set(gca, 'FontSize', FS_AXIS, 'LineWidth', LW_AXIS, 'FontWeight', 'bold');
grid on;

%%
% Compute Reference solution (FDTD)
u_ref = run_fdtd(u0, v0, Lx, T, c, dh, dt, gamma, nu);

% Random initial guess for Schwarz iteration (fixed seed for consistency)
rng(42); 
u_init = rand(Nx, Nt);

% Get SWR error surface for this experiment
fprintf('Computing SWR error surface...\n');
error_surface = swr_error_surface(N, a, M, T, c, dh, dt, gamma, nu, k, THETA1, THETA2, u0, v0, u_init, u_ref);

%% Extract optimal parameters for the single experiment (SWR)
[min_val, idx_min] = min(error_surface, [], 'all', 'linear');
[min_i, min_j] = ind2sub(size(error_surface), idx_min);
min_theta1_swr = THETA1(min_i, min_j);
min_theta2_swr = THETA2(min_i, min_j);
fprintf('Grid search best: p=%.4f, q=%.4f, val=%.4e\n', min_theta1_swr, min_theta2_swr, min_val);

% =========================================================================
% REFINE EXPERIMENTAL OPTIMUM using fminsearch
% =========================================================================
fprintf('Refining experimental optimum using fminsearch...\n');
% 2. Define Experimental Objective Function
objfun_exp = @(x) obj_exp(x, u0, v0, N, a, M, T, c, dh, dt, gamma, nu, k, u_init, u_ref);

% 3. Run fminsearch starting from the Grid Best
x0_refine = [min_theta1_swr, min_theta2_swr];
options_exp = optimset('Display','iter', 'TolX',1e-4, 'TolFun',1e-4);
[x_refined, val_refined] = fminsearch(objfun_exp, x0_refine, options_exp);

% 4. Update the variables used for plotting
fprintf('Refined experimental: p=%.4f, q=%.4f, val=%.4e\n', x_refined(1), x_refined(2), val_refined);
min_theta1_swr = x_refined(1);
min_theta2_swr = x_refined(2);
min_val_swr    = val_refined;

%% DATA MODIFICATION FOR PLOTTING
% "Hack": Find the nearest point on the coarse grid and force it to take
% the value of the refined minimum. This forces contourf to color it deeply.
dist_sq = (THETA1 - min_theta1_swr).^2 + (THETA2 - min_theta2_swr).^2;
[~, idx_nearest] = min(dist_sq(:));
error_surface(idx_nearest) = min_val_swr;

%% Find optimal points for theoretical surfaces via fminsearch

% Initial guess point (also used in plots)
theta1_initial = 1/c;
theta2_initial = 0;
x0 = [theta1_initial, theta2_initial];

% fminsearch options
base_options = optimset('Display','off', 'TolX',1e-8, 'TolFun',1e-8);

% p=2 contraction factor minimization
objfun_p2 = @(x) obj_L2( N, T, dt, J, c, gamma, nu, a, M, x(1), x(2) );
[x_opt_p2, min_val_p2] = fminsearch(objfun_p2, x0, base_options);
min_theta1_p2 = clip(x_opt_p2(1), theta1_min, theta1_max);
min_theta2_p2 = clip(x_opt_p2(2), theta2_min, theta2_max);

% p=∞ contraction factor minimization
objfun_inf = @(x) obj_Linf( N, T, dt, J, c, gamma, nu, a, M, x(1), x(2) );
[x_opt_inf, min_val_inf] = fminsearch(objfun_inf, x0, base_options);
min_theta1_inf = clip(x_opt_inf(1), theta1_min, theta1_max);
min_theta2_inf = clip(x_opt_inf(2), theta2_min, theta2_max);

%% Create main comparison figure
theta1_range_plot = THETA1(1, :);
theta2_range_plot = THETA2(:, 1)';

figure('Name', 'SWR Error Surface', 'Color', 'w');
% Set figure size suitable for posters
set(gcf, 'Position', [100 100 800 600]); 

contourf(theta1_range_plot, theta2_range_plot, log10(error_surface), 50, ...
         'LineStyle', 'none', 'HandleVisibility', 'off');
axis xy;
hold on;

% Plot Markers with Large Sizes and Thick Borders
plot(theta1_initial, theta2_initial, 'ks', ...
    'MarkerSize', MS_MARK, 'MarkerFaceColor', 'k', 'LineWidth', 2, ...
    'DisplayName', 'Initial guess');

plot(min_theta1_swr, min_theta2_swr, 's', ...
    'MarkerSize', MS_MARK, 'Color', "k", 'LineWidth', 2, ...
    'MarkerFaceColor', "w", ...
    'DisplayName', 'SWR error opt.');

plot(min_theta1_p2, min_theta2_p2, 'mo', ...
    'MarkerSize', MS_MARK, 'MarkerFaceColor', 'm', 'LineWidth', 2, ...
    'DisplayName', 'Spectral opt. L_2');

plot(min_theta1_inf, min_theta2_inf, 'g^', ...
    'MarkerSize', MS_MARK, 'MarkerFaceColor', 'g', 'LineWidth', 2, ...
    'DisplayName', 'Spectral opt. L_\infty');

% Axis Labels
xlabel('p', 'FontSize', FS_LABEL, 'FontWeight', 'bold');
ylabel('q', 'FontSize', FS_LABEL, 'FontWeight', 'bold');

% Legend
% lgd = legend('show', 'Location', 'NorthEast');
% set(lgd, 'FontSize', 22, 'FontWeight', 'bold'); % Slightly smaller than axis labels

% Colorbar
colormap('parula');
cb = colorbar;
cb.FontSize = FS_AXIS;
cb.FontWeight = 'bold';
% ylabel(cb, 'log_{10}(Error)', 'FontSize', FS_LABEL, 'FontWeight', 'bold');

% Color Limits
if gamma == 0
    clim([min(log10(error_surface),[],'all'),3])
elseif nu == 0
    clim([min(log10(error_surface),[],'all'),8])
end

% Set Axis properties
set(gca, 'FontSize', FS_AXIS, 'LineWidth', LW_AXIS, 'FontWeight', 'bold');

%% Create separate figure for the Horizontal Legend ONLY
fig_lgd = figure('Name', 'Horizontal Legend', 'Color', 'w');
set(gcf, 'Position', [100 100 1400 100]); % Wider window for more space

hold on;
% Use \quad or \hspace{1cm} for explicit horizontal spacing in LaTeX
plot(nan, nan, 'ks', 'MarkerSize', MS_MARK, 'MarkerFaceColor', 'k', 'LineWidth', 2, ...
    'DisplayName', 'Initial guess \hspace{0.5cm}'); 
plot(nan, nan, 's', 'MarkerSize', MS_MARK, 'Color', "k", 'LineWidth', 2, ...
    'MarkerFaceColor', "w", ...
    'DisplayName', 'SWR error opt. \hspace{0.5cm}');
plot(nan, nan, 'mo', 'MarkerSize', MS_MARK, 'MarkerFaceColor', 'm', 'LineWidth', 2, ...
    'DisplayName', 'Spectral opt. $L_2$ \hspace{0.5cm}');
plot(nan, nan, 'g^', 'MarkerSize', MS_MARK, 'MarkerFaceColor', 'g', 'LineWidth', 2, ...
    'DisplayName', 'Spectral opt. $L_\infty$');

% Configure the legend
lgd_h = legend('show', 'Orientation', 'horizontal', 'Location', 'north');
set(lgd_h, 'FontSize', 22, 'FontWeight', 'bold', 'Interpreter', 'latex');

legend('boxoff');
axis off;

%% Display summary
fprintf('\n=== RESULTS SUMMARY ===\n');
fprintf('\nComparison of Methods:\n');
fprintf('Method          theta1     theta2\n');
fprintf('----------------------------------------------------------\n');
fprintf('SWR error opt. %8.4f     %8.4f\n', min_theta1_swr, min_theta2_swr);
fprintf('L2             %8.4f     %8.4f\n', min_theta1_p2,  min_theta2_p2);
fprintf('Linf           %8.4f     %8.4f\n', min_theta1_inf, min_theta2_inf);