%% Compare SWR error surface with contraction factor surfaces for single T
clear all; close all; clc;

addpath("utils\")

%% Parameters
P = get_sim_params_1D();

N  = P.N;
a  = P.a;
M  = P.M;
b  = P.b;
Lx = P.Lx;
T  = P.T;

c  = P.c;
gamma = P.gamma; % ignore
nu = P.nu; % ignore

dh = P.dh;
dt = P.dt;
J  = P.J;

%%

% Viscous damping case
% gamma = 4;
% gamma = 8;
gamma = 10;
% gamma = 12;
nu = 0;
% k = 10;
k = 5*N;
% k = 40;

% % Viscoelastic damping case
% gamma = 0;
% % nu = 0.001;
% % nu = 0.01;
% % nu = 0.05;
% nu = 0.1;
% % k = 10;
% k = 5*N;
% % k = 20;

% Parameter ranges
% theta1_range = linspace(0, 1.2, 25);
% theta2_range = linspace(-4, 8, 25);

% theta1_range = linspace(0, 1.2, 13);
% theta2_range = linspace(-4, 8, 13);

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
gaussian_normalized = @(x) gaussian(x, Lx/4, Lx/20);
    
% Create the final initial condition function
u0 = @(x) gaussian_normalized(x);
v0 = @(x) 0.*x;

% Compute Reference solution (FDTD)
u_ref = run_fdtd_1D(u0, v0, Lx, T, c, dh, dt, gamma, nu);

% Random initial guess for Schwarz iteration (fixed seed for consistency)
rng(42); 
u_init = rand(Nx, Nt);

% Get SWR error surface for this experiment
fprintf('Computing SWR error surface...\n');
error_surface = swr_error_surface_1D(N, a, M, T, c, dh, dt, gamma, nu, k, THETA1, THETA2, u0, v0, u_init, u_ref);

% % Get contraction factor surfaces
% fprintf('Computing contraction factor surfaces...\n');
% Ly = 0; y_mode = 0; % 1D
% [Z_p2, Z_inf] = contraction_surface(N, a, M, Ly, y_mode, T, c, dh, dt, gamma, nu, J, THETA1, THETA2);

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
% Note: We use the current 'k' defined in parameters for the optimization
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

% Calculate squared distance from refined point to all grid points
dist_sq = (THETA1 - min_theta1_swr).^2 + (THETA2 - min_theta2_swr).^2;

% Find index of the closest grid point
[~, idx_nearest] = min(dist_sq(:));

% Replace that grid point's value with the refined minimum
error_surface(idx_nearest) = min_val_swr;

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

%% Create main comparison figure
theta1_range_plot = THETA1(1, :);
theta2_range_plot = THETA2(:, 1)';

figure()
contourf(theta1_range_plot, theta2_range_plot, log10(error_surface), 50, 'LineStyle', 'none', 'HandleVisibility', 'off');
axis xy;
hold on;

plot(theta1_initial, theta2_initial, 'ks', ...
    'MarkerSize', 10, 'MarkerFaceColor', 'k', ...
    'DisplayName', 'Initial guess');

plot(min_theta1_swr, min_theta2_swr, 's', ...
    'MarkerSize', 10, 'Color', "k", ...
    'MarkerFaceColor', "w", ...
    'DisplayName', 'SWR error opt.');

plot(min_theta1_p2, min_theta2_p2, 'mo', ...
    'MarkerSize', 10, 'MarkerFaceColor', 'm', ...
    'DisplayName', 'Spectral opt. L2');

plot(min_theta1_inf, min_theta2_inf, 'g^', ...
    'MarkerSize', 10, 'MarkerFaceColor', 'g', ...
    'DisplayName', 'Spectral opt. L∞');

colorbar;
xlabel('p', 'FontSize', 16);
ylabel('q', 'FontSize', 16);
% title('Log-scale Error Surface');
% legend('show', 'Location', 'SouthWest', 'FontSize', 14);
% legend('show', 'Location', 'NorthEast', 'FontSize', 14);
legend('show', 'Location', 'NorthWest', 'FontSize', 14);
% legend('show', 'Location', 'best', 'FontSize', 14);
colormap('parula');
if gamma == 0
    clim([min(log10(error_surface),[],'all'),3])
elseif nu == 0
    clim([min(log10(error_surface),[],'all'),8])
end
% Make tick labels bigger
set(gca, 'FontSize', 18);

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
fprintf('Method          theta1     theta2\n');
fprintf('----------------------------------------------------------\n');
fprintf('SWR error opt. %8.4f     %8.4f\n', min_theta1_swr, min_theta2_swr);
fprintf('L2             %8.4f     %8.4f\n', min_theta1_p2,  min_theta2_p2);
fprintf('Linf           %8.4f     %8.4f\n', min_theta1_inf, min_theta2_inf);
