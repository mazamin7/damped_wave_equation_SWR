clear all; close all; clc;

addpath("utils\")

% SWR parameters
N = 2;
% N = 8;
a = 0.3;
% M = 0.05;
M = 0.1;
% M = 0.2;
% M = 0.4;
Lx = N*a + M;

% Parameters
c = 1;              % Wave speed
dh = 0.01;          % Spatial step in x
dt = 0.01;
% dh = 0.001;
% dt = 0.001;


% Initial conditions components
gaussian = @(r,mu,sigma) 1/(2*pi*sigma^2) * exp(-(r-mu).^2/(2*sigma^2));

% Define n_max
n_max = 100;

% Compute max values for normalization
x_test = linspace(0,Lx,1000);
gaussian_max = max(gaussian(x_test, Lx/4, Lx/20));

% Compute max of sine sum
sine_sum_vals = zeros(size(x_test));
for n = 1:n_max
    sine_sum_vals = sine_sum_vals + sin(n*pi*x_test/Lx);
end
sine_sum_max = max(abs(sine_sum_vals));

% Create normalized components as function handles
gaussian_normalized = @(x) gaussian(x, Lx/4, Lx/20) / gaussian_max;

% Create normalized sine sum function handle
sine_sum_normalized = @(x) arrayfun(@(xi) sum(sin((1:n_max)' * pi * xi / Lx)) / sine_sum_max, x);

% Create the final initial condition function
u0 = @(x) gaussian_normalized(x) + sine_sum_normalized(x);
v0 = @(x) 0.*x;


% Define multiple test cases

% % viscous damping case T1
% gamma = 1;
% nu = 0;
% % T = 1;
% T = 5;
% theta_sets = [
%     1/c,     0;      % Initial guess
%     0.45,    1.5;    % Numerical optimization
%     1,       1.5;    % Spectral optimization p=2
%     1,       0.5     % Spectral optimization p=∞
% ];
% % k = 10;
% % k = 80;
% k = 40*N;

% viscous damping case T1
gamma = 5;
nu = 0;
% T = 1;
T = 5;
theta_sets = [
    1/c,     0;      % Initial guess
    0.8,     0.5;      % Numerical optimization
    1,       2.5; % Spectral optimization p=2
    1.05,    2.5  % Spectral optimization p=∞
];
% k = 10;
% k = 80;
k = 40*N;

% % viscoelastic damping case T1
% gamma = 0;
% nu = 0.1;
% % T = 1;
% T = 5;
% theta_sets = [
%     1/c,     0;      % Initial guess
%     0.55,    2;   % Numerical optimization
%     0.55,    2.5; % Spectral optimization p=2
%     0.5,     2    % Spectral optimization p=∞
% ];
% % k = 10;
% % k = 80;
% % k = 40*N;
% k = 40*N;

% % viscoelastic damping case T1
% gamma = 0;
% nu = 0.5;
% % T = 1;
% T = 5;
% theta_sets = [
%     1/c,     0;      % Initial guess
%     0.15,    3.5; % Numerical optimization
%     0.1,     5;   % Spectral optimization p=2
%     0.15,    3.5  % Spectral optimization p=∞
% ];
% % k = 10;
% % k = 80;
% % k = 40*N;
% k = 40*N;


disp(T*c/M)


numSets = size(theta_sets, 1);
res_hist = cell(numSets, 1);
final_errors = zeros(numSets, 1);

%% Compute reference solution using modular function
fprintf('Computing reference FDTD solution...\n');
u_ref = run_fdtd_1D(u0, v0, Lx, T, c, dh, dt, gamma, nu);

Nx = round(Lx / dh) + 1;
Nt = round(T / dt) + 1;

u_init = rand(Nx,Nt);

%% Loop over different theta sets
for s = 1:numSets
    theta1 = theta_sets(s, 1);
    theta2 = theta_sets(s, 2);
    
    fprintf('Running simulation for theta1=%.6f, theta2=%.6f (%d/%d)\n', ...
            theta1, theta2, s, numSets);
    
    % Run SWR using modular function
    [ud, final_res, res_history] = run_swr_1D(u0, v0, N, a, M, T, c, dh, dt, gamma, nu, theta1, theta2, k, u_init, u_ref);
    
    % Store results
    res_hist{s} = res_history;
    final_errors(s) = final_res;
    
    fprintf('Finished theta set %d/%d after %d iterations (final error=%.3e)\n', ...
        s, numSets, k, final_res);
end

%% Plot results
figure;

% Define colors and line styles for consistent ordering
colors = ['k', 'r', 'm', 'g'];
line_styles = {'-', ':', '--', '-.'};
legend_labels = {'Initial guess', 'SWR error opt.', ...
                 'Spectral opt. L2', 'Spectral opt. L\infty'};

for s = 1:numSets
    semilogy(1:length(res_hist{s}), res_hist{s}, ...
             [colors(s) line_styles{s}], 'LineWidth', 2);
    hold on;
end

xlabel('Iteration','FontSize', 16);
ylabel('Error','FontSize', 16);
xlim([0,k])
xticks(0:10:k)
ylim([1e-15,1e5])
yticks([1e-15,1e-10,1e-5,1e0,1e5])
% title('WR Convergence with Different Robin Interface Parameters', 'FontSize', 24);

% legend(legend_labels{1:numSets}, 'Location', 'NorthEast', 'FontSize', 14);
legend(legend_labels{1:numSets}, 'Location', 'SouthWest', 'FontSize', 14);
grid on;

% Make tick labels bigger
set(gca, 'FontSize', 18);

%% Display final errors
fprintf('\n=== Final errors Summary ===\n');
fprintf('Theta Set\t\tFinal error (after %d iterations)\n', k);
fprintf('--------------------------------------------\n');
for s = 1:numSets
    fprintf('%s\t\t%.3e\n', ...
            legend_labels{s}, final_errors(s));
end