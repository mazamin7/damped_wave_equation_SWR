% QUICK_TEST_RUN_SWR_1D - Compare FDTD, SWR Dirichlet, and SWR Robin in 1D
clear all; close all; clc;

addpath("utils\")

fprintf('Comparing FDTD, SWR Dirichlet, and SWR Robin in 1D...\n');

% Parameters
P = get_sim_params_1D();

N  = P.N;
a  = P.a;
M  = P.M;
b  = P.b;
Lx = P.Lx;
T  = P.T;

c  = P.c;
gamma = P.gamma;
nu = P.nu;

dh = P.dh;
dt = P.dt;
J  = P.J;

%% --- snap dh,dt to the grid used everywhere ---
Nx = round(Lx/dh) + 1;
Nt = round(T/dt) + 1;

theta1 = 1/c;
theta2 = 0;
k      = 50*N;

fprintf('Number of iterations needed with Dirichlet int.cond.: %d\n', ceil(c*T/M/2*N));

% Initial conditions
gaussian = @(r,mu,sigma) 1/(2*pi*sigma^2) * exp(-(r-mu).^2/(2*sigma^2));
n_max = 100;

x_test = linspace(0,Lx,1000);
gaussian_max = max(gaussian(x_test, Lx/4, Lx/20));

sine_sum_vals = zeros(size(x_test));
for n = 1:n_max
    sine_sum_vals = sine_sum_vals + sin(n*pi*x_test/Lx);
end
sine_sum_max = max(abs(sine_sum_vals));

gaussian_normalized = @(x) gaussian(x, Lx/4, Lx/20) / gaussian_max;
sine_sum_normalized = @(x) arrayfun(@(xi) sum(sin((1:n_max)' * pi * xi / Lx)) / sine_sum_max, x);

u0 = @(x) gaussian_normalized(x) + 0*sine_sum_normalized(x);
v0 = @(x) 0.*x;

% Plot initial conditions
x_ic = linspace(0,Lx,1000);
u0_vals = u0(x_ic);
v0_vals = v0(x_ic);

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

%% Compute reference solution with FDTD
fprintf('Computing reference solution with FDTD...\n');
u_ref = run_fdtd_1D(u0, v0, Lx, T, c, dh, dt, gamma, nu);

% sizes consistent with snapped steps
Nx = size(u_ref,1);
Nt = size(u_ref,2);

rng(123);
u_init = rand(Nx,Nt);    % match sizes of u_ref

%% Run SWR with Dirichlet interfaces
fprintf('Running SWR with Dirichlet interfaces...\n');
tic;
[ud_dirichlet, final_res_dirichlet, res_history_dirichlet] = ...
    run_swr_dirichlet_1D(u0, v0, N, a, M, T, c, dh, dt, gamma, nu, k, u_init, u_ref);
time_dirichlet = toc;

fprintf('Dirichlet Results:\n');
fprintf('  Iterations: %d\n', k);
fprintf('  Final residual: %.3e\n', final_res_dirichlet);
fprintf('  Time: %.2f seconds\n', time_dirichlet);

%% Run SWR with Robin interfaces
fprintf('Running SWR with Robin interfaces...\n');
tic;
[ud_robin, final_res_robin, res_history_robin] = ...
    run_swr_1D(u0, v0, N, a, M, T, c, dh, dt, gamma, nu, theta1, theta2, k, u_init, u_ref);
time_robin = toc;

fprintf('Robin Results:\n');
fprintf('  Iterations: %d\n', k);
fprintf('  Final residual: %.3e\n', final_res_robin);
fprintf('  Time: %.2f seconds\n', time_robin);

%% Plot comparison results
x_axis = linspace(0, Lx, Nx);
t_axis = linspace(0, T, Nt);

figure()
subplot(2,3,1)
mesh(x_axis, t_axis, u_ref')
xlabel('x'); ylabel('t'); zlabel('u')
title('FDTD Reference'); colorbar

subplot(2,3,2)
mesh(x_axis, t_axis, ud_dirichlet')
xlabel('x'); ylabel('t'); zlabel('u')
title('SWR Dirichlet'); colorbar

subplot(2,3,3)
mesh(x_axis, t_axis, ud_robin')
xlabel('x'); ylabel('t'); zlabel('u')
title('SWR Robin'); colorbar

subplot(2,3,4)
error_dirichlet = abs(ud_dirichlet - u_ref);
mesh(x_axis, t_axis, error_dirichlet')
xlabel('x'); ylabel('t'); zlabel('Error')
title('Dirichlet Error'); colorbar

subplot(2,3,5)
error_robin = abs(ud_robin - u_ref);
mesh(x_axis, t_axis, error_robin')
xlabel('x'); ylabel('t'); zlabel('Error')
title('Robin Error'); colorbar

subplot(2,3,6)
semilogy(1:length(res_history_dirichlet), res_history_dirichlet, 'o-', 'LineWidth', 2, 'DisplayName', 'Dirichlet')
hold on
semilogy(1:length(res_history_robin), res_history_robin, 's-', 'LineWidth', 2, 'DisplayName', 'Robin')
xlabel('Iteration'); ylabel('Relative Residual')
title('SWR Convergence History'); legend('show','Location','best'); grid on
sgtitle('Method Comparison: FDTD vs SWR with Different Interface Conditions')

%% Plot time evolution at specific points
figure()
x_points = [Lx/4, Lx/2, 3*Lx/4];
point_indices = round(x_points / dh) + 1;   % dh is snapped

subplot(2,3,1)
for i = 1:length(point_indices)
    plot(t_axis, u_ref(point_indices(i),:), 'LineWidth', 2, 'DisplayName', sprintf('x=%.2f', x_points(i))); hold on
end
xlabel('t'); ylabel('u'); title('FDTD Reference'); legend('show'); grid on

subplot(2,3,2)
for i = 1:length(point_indices)
    plot(t_axis, ud_dirichlet(point_indices(i),:), 'LineWidth', 2, 'DisplayName', sprintf('x=%.2f', x_points(i))); hold on
end
xlabel('t'); ylabel('u'); title('SWR Dirichlet'); legend('show'); grid on

subplot(2,3,3)
for i = 1:length(point_indices)
    plot(t_axis, ud_robin(point_indices(i),:), 'LineWidth', 2, 'DisplayName', sprintf('x=%.2f', x_points(i))); hold on
end
xlabel('t'); ylabel('u'); title('SWR Robin'); legend('show'); grid on

subplot(2,3,4)
for i = 1:length(point_indices)
    plot(t_axis, abs(u_ref(point_indices(i),:) - ud_dirichlet(point_indices(i),:)), 'LineWidth', 2, ...
        'DisplayName', sprintf('x=%.2f', x_points(i))); hold on
end
xlabel('t'); ylabel('Error'); title('Dirichlet Error'); legend('show'); grid on

subplot(2,3,5)
for i = 1:length(point_indices)
    plot(t_axis, abs(u_ref(point_indices(i),:) - ud_robin(point_indices(i),:)), 'LineWidth', 2, ...
        'DisplayName', sprintf('x=%.2f', x_points(i))); hold on
end
xlabel('t'); ylabel('Error'); title('Robin Error'); legend('show'); grid on

subplot(2,3,6)
plot(x_axis, u_ref(:,end), 'k-', 'LineWidth', 3, 'DisplayName', 'FDTD Reference'); hold on
plot(x_axis, ud_dirichlet(:,end), 'r--', 'LineWidth', 2, 'DisplayName', 'SWR Dirichlet');
plot(x_axis, ud_robin(:,end), 'b:',  'LineWidth', 2, 'DisplayName', 'SWR Robin');
xlabel('x'); ylabel('u'); title('Final Time Slice Comparison'); legend('show','Location','best'); grid on
sgtitle('Time Evolution at Selected Points: Method Comparison')

%% Detailed convergence comparison
figure()
semilogy(1:length(res_history_dirichlet), res_history_dirichlet, 'o-', 'LineWidth', 2, 'DisplayName', 'Dirichlet'); hold on
semilogy(1:length(res_history_robin), res_history_robin, 's-', 'LineWidth', 2, 'DisplayName', 'Robin');
xlabel('Iteration'); ylabel('Relative Residual'); title('Convergence History'); legend('show'); grid on

%% Performance summary
fprintf('\n=== PERFORMANCE SUMMARY ===\n');
fprintf('Method            Final Residual    Time (s)    Iterations\n');
fprintf('--------------------------------------------------------\n');
fprintf('SWR Dirichlet     %.3e         %.2f        %d\n', final_res_dirichlet, time_dirichlet, k);
fprintf('SWR Robin         %.3e         %.2f        %d\n', final_res_robin, time_robin, k);

if final_res_dirichlet > 0 && final_res_robin > 0
    if final_res_dirichlet < final_res_robin
        fprintf('Dirichlet is %.2fx more accurate than Robin\n', final_res_robin / final_res_dirichlet);
    else
        fprintf('Robin is %.2fx more accurate than Dirichlet\n', final_res_dirichlet / final_res_robin);
    end
end

%% Subdomain boundaries
figure()
x_ticks_arr = zeros(2,N);
aj = @(j) a*(j-1);
bj = @(j) aj(j+1) + M;
for iii = 1:N
    x_ticks_arr(:,iii) = [aj(iii), bj(iii)];
end
x_ticks_arr = reshape(x_ticks_arr,2*N,1);

labels = cell(1, 2*N);
for iii = 1:N
    labels{2*iii-1} = sprintf('a%d', iii);
    labels{2*iii}   = sprintf('b%d', iii);
end
[x_ticks_arr, sortIdx] = sort(x_ticks_arr);
labels = labels(sortIdx);

plot(x_ticks_arr, zeros(size(x_ticks_arr)), 'rx', 'MarkerSize', 10, 'LineWidth', 2)
xlabel('x'); title('Subdomain Interface Locations'); grid on
for i = 1:length(x_ticks_arr)
    text(x_ticks_arr(i), 0.1, labels{i}, 'HorizontalAlignment', 'center')
end
ylim([-0.5, 0.5])
