% QUICK_TEST_RUN_SWR - Simple test for run_swr
clear all; close all; clc;

addpath("utils\")

fprintf('Quick test of run_swr...\n');

% Parameters
% N = 2;
N = 4;
a = 0.3;
M = 0.1;
Lx = N*a + M;
% T = 50.0;
T = 5;
c = 1.0;
dh = 0.01;
dt = 0.01;
gamma = 1;
nu = 0;
% gamma = 0;
% nu = 0.1;
theta1 = 1;
theta2 = 0;
% theta1 = 0;
% theta2 = 1;
% k = 100;  % Fixed number of iterations
k = 50*N;


% Initial conditions components
gaussian = @(r,mu,sigma) 1/(2*pi*sigma^2) * exp(-(r-mu).^2/(2*sigma^2));

% Define n_max - you can adjust this value as needed
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
v0 = @(x) 0;


% Evaluate initial conditions
x_ic = linspace(0,Lx,1000);
u0_vals = u0(x_ic);
v0_vals = v0(x_ic);

% Plot initial conditions
figure()
subplot(2,1,1)
plot(x_ic, u0_vals, 'b-', 'LineWidth', 2)
xlabel('Position (x)')
ylabel('u_0(x)')
title('Initial Displacement u_0(x)')
grid on

subplot(2,1,2)
plot(x_ic, v0_vals, 'r-', 'LineWidth', 2)
xlabel('Position (x)')
ylabel('v_0(x)')
title('Initial Velocity v_0(x)')
grid on
sgtitle('Initial Conditions')


fprintf('Computing reference...\n');
u_ref = run_fdtd_1D(u0, v0, Lx, T, c, dh, dt, gamma, nu);


Nx = round(Lx / dh) + 1;
Nt = floor(T / dt);

u_init = rand(Nx,Nt);

fprintf('Running SWR...\n');
tic;
[ud, final_res, res_history] = run_swr_1D(u0, v0, N, a, M, T, c, dh, dt, gamma, nu, theta1, theta2, k, u_init, u_ref);
time_elapsed = toc;

fprintf('Results:\n');
fprintf('  Iterations: %d\n', k);
fprintf('  Final residual: %.3e\n', final_res);
fprintf('  Time: %.2f seconds\n', time_elapsed);

% Plot residual history
figure()
semilogy(1:k, res_history, 'o-', 'LineWidth', 2)
xlabel('Iteration')
ylabel('Relative Residual')
title('SWR Convergence History')
grid on

figure()
mesh(u_ref)
title('Reference Solution');

figure()
mesh(ud)
title('SWR Solution');

% Plot error
figure()
mesh(abs(ud - u_ref)'/max(abs(u_ref), [], 'all'))
xlabel('space')
ylabel('time')
title('Absolute Error');