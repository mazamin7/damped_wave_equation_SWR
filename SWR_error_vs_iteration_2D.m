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
Ly = 0.1;

% Parameters
c = 1;              % Wave speed
dh = 0.01;          % Spatial step in x
dt = 0.7*dh/c;
% dh = 0.001;
% dt = 0.001;


% Initial conditions
u0 = @(x,y) 0.*x.*y;
v0 = @(x,y) exp(-sqrt((x-Lx/3).^2+(y-Ly/3).^2));


% Define multiple test cases

% % viscous damping case T1
% gamma = 1;
% nu = 0;
% % T = 1;
% T = 5;
% theta_sets = [
%     1/c,     0;      % Initial guess
%     0.7,     -4;      % Numerical optimization
%     0.9,       8; % Spectral optimization p=2
%     0.9,       0  % Spectral optimization p=∞
% ];
% % k = 10;
% % k = 80;
% k = 40*N;

% viscoelastic damping case T1
gamma = 0;
nu = 1;
% T = 1;
T = 5;
theta_sets = [
    1/c,     0;      % Initial guess
    0.1,     8; % Numerical optimization
    0,       8; % Spectral optimization p=2
    0,       8  % Spectral optimization p=∞
];
% k = 10;
% k = 80;
k = 10*N;
% k = 40*N;
% k = 60*N;


disp(T*c/M)


numSets = size(theta_sets, 1);
res_hist = cell(numSets, 1);
final_errors = zeros(numSets, 1);

%% Compute reference solution using modular function
fprintf('Computing reference FDTD solution...\n');
u_ref = run_fdtd_2D(u0, v0, Lx, Ly, T, c, dh, dt, gamma, nu);

Nx = round(Lx / dh) + 1;
Ny = round(Ly / dh) + 1;
Nt = floor(T / dt);

u_init = rand(Nx,Ny,Nt);

%% Loop over different theta sets
for s = 1:numSets
    theta1 = theta_sets(s, 1);
    theta2 = theta_sets(s, 2);
    
    fprintf('Running simulation for theta1=%.6f, theta2=%.6f (%d/%d)\n', ...
            theta1, theta2, s, numSets);
    
    % Run SWR using modular function
    [ud, final_res, res_history] = run_swr_2D(u0, v0, N, a, M, Ly, T, c, dh, dt, gamma, nu, theta1, theta2, k, u_init, u_ref);
    
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
xticks(0:10:k)
% title('WR Convergence with Different Robin Interface Parameters', 'FontSize', 24);

legend(legend_labels{1:numSets}, 'Location', 'NorthEast', 'FontSize', 14);
% legend(legend_labels{1:numSets}, 'Location', 'SouthWest', 'FontSize', 14);
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