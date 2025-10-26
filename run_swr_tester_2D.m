% QUICK_TEST_RUN_SWR_2D - Compare FDTD, SWR Dirichlet, and SWR Robin
clear all; close all; clc;

addpath("utils\")

fprintf('Comparing FDTD, SWR Dirichlet, and SWR Robin in 2D...\n');

% Parameters
% N = 2;
% N = 4;
% N = 8;
N = 16;

a = 0.3;
M = 0.1;
Lx = N*a + M;
Ly = 0.1;
T = 5;
c = 1.0;
% dh = 0.01;
dh = 0.05;
dt = 0.7*dh/c;

% % viscous
% gamma = 1;
% nu = 0;
% theta1 = 1/c;
% theta2 = 0;
% k = 50*N;

% viscoelastic
gamma = 0;
nu = 1;
theta1 = 0;
theta2 = 1/a;
% theta1 = 0;
% theta2 = 8;
k = 50*N;


fprintf('Number of iterations needed with Dirichlet int.cond.: %d\n', ceil(c*T/M/2*N));


% % Initial conditions
% u0 = @(x,y) sin(3*pi.*x/Lx) .* sin(5*pi.*y/Ly);
% v0 = @(x,y) 0.*x.*y;

u0 = @(x,y) 0.*x.*y;
v0 = @(x,y) exp(-sqrt((x-Lx/3).^2+(y-Ly/3).^2));

% Plot initial conditions
x_ic = linspace(0, Lx, 100);
y_ic = linspace(0, Ly, 100);
[X, Y] = meshgrid(x_ic, y_ic);
U0 = u0(X, Y);

figure()
subplot(1,2,1)
surf(X, Y, U0)
shading interp
xlabel('x')
ylabel('y')
zlabel('u_0(x,y)')
title('Initial Displacement')

subplot(1,2,2)
imagesc(x_ic, y_ic, U0)
xlabel('x')
ylabel('y')
title('Initial Displacement (top view)')
colorbar
sgtitle('Initial Conditions')

%% Compute reference solution with FDTD
fprintf('Computing reference solution with FDTD...\n');
u_ref = run_fdtd_2D(u0, v0, Lx, Ly, T, c, dh, dt, gamma, nu);
u_init = rand(size(u_ref));

%% Run SWR with Dirichlet interfaces
fprintf('Running SWR with Dirichlet interfaces...\n');
tic;
[ud_dirichlet, final_res_dirichlet, res_history_dirichlet] = run_swr_dirichlet_2D(u0, v0, N, a, M, Ly, T, c, dh, dt, gamma, nu, theta1, theta2, k, u_init, u_ref);
time_dirichlet = toc;

fprintf('Dirichlet Results:\n');
fprintf('  Iterations: %d\n', k);
fprintf('  Final residual: %.3e\n', final_res_dirichlet);
fprintf('  Time: %.2f seconds\n', time_dirichlet);

%% Run SWR with Robin interfaces
fprintf('Running SWR with Robin interfaces...\n');
tic;
[ud_robin, final_res_robin, res_history_robin] = run_swr_2D(u0, v0, N, a, M, Ly, T, c, dh, dt, gamma, nu, theta1, theta2, k, u_init, u_ref);
time_robin = toc;

fprintf('Robin Results:\n');
fprintf('  Iterations: %d\n', k);
fprintf('  Final residual: %.3e\n', final_res_robin);
fprintf('  Time: %.2f seconds\n', time_robin);

%% Plot comparison results
% Final time solutions comparison
figure()
subplot(2,3,1)
imagesc(squeeze(u_ref(:,:,end)))
title('FDTD Reference (final time)')
colorbar
axis equal tight

subplot(2,3,2)
imagesc(squeeze(ud_dirichlet(:,:,end)))
title('SWR Dirichlet (final time)')
colorbar
axis equal tight

subplot(2,3,3)
imagesc(squeeze(ud_robin(:,:,end)))
title('SWR Robin (final time)')
colorbar
axis equal tight

% Error maps
subplot(2,3,4)
error_dirichlet = abs(ud_dirichlet - u_ref);
imagesc(squeeze(error_dirichlet(:,:,end)))
title('Dirichlet Error (final time)')
colorbar
axis equal tight

subplot(2,3,5)
error_robin = abs(ud_robin - u_ref);
imagesc(squeeze(error_robin(:,:,end)))
title('Robin Error (final time)')
colorbar
axis equal tight

% Convergence history
subplot(2,3,6)
semilogy(1:length(res_history_dirichlet), res_history_dirichlet, 'o-', 'LineWidth', 2, 'DisplayName', 'Dirichlet')
hold on
semilogy(1:length(res_history_robin), res_history_robin, 's-', 'LineWidth', 2, 'DisplayName', 'Robin')
xlabel('Iteration')
ylabel('Relative Residual')
title('SWR Convergence History')
legend('show', 'Location', 'best')
grid on
sgtitle('Method Comparison: FDTD vs SWR with Different Interface Conditions')

%% Plot time evolution at center line for all methods
figure()
y_center = round(Ly/2/dh);
x_axis = linspace(0, Lx, size(u_ref,1));
t_axis = linspace(0, T, size(u_ref,3));

u_ref_center = squeeze(u_ref(:, y_center, :))';
ud_dirichlet_center = squeeze(ud_dirichlet(:, y_center, :))';
ud_robin_center = squeeze(ud_robin(:, y_center, :))';

% Reference solution
subplot(2,3,1)
surf(x_axis, t_axis, u_ref_center)
shading interp
xlabel('x')
ylabel('t')
zlabel('u')
title('FDTD Reference (center y)')

% Dirichlet solution
subplot(2,3,2)
surf(x_axis, t_axis, ud_dirichlet_center)
shading interp
xlabel('x')
ylabel('t')
zlabel('u')
title('SWR Dirichlet (center y)')

% Robin solution
subplot(2,3,3)
surf(x_axis, t_axis, ud_robin_center)
shading interp
xlabel('x')
ylabel('t')
zlabel('u')
title('SWR Robin (center y)')

% Error - Dirichlet
subplot(2,3,4)
error_dirichlet_center = abs(u_ref_center - ud_dirichlet_center);
surf(x_axis, t_axis, error_dirichlet_center)
shading interp
xlabel('x')
ylabel('t')
zlabel('Error')
title('Dirichlet Error (center y)')

% Error - Robin
subplot(2,3,5)
error_robin_center = abs(u_ref_center - ud_robin_center);
surf(x_axis, t_axis, error_robin_center)
shading interp
xlabel('x')
ylabel('t')
zlabel('Error')
title('Robin Error (center y)')

% Final time slice comparison
subplot(2,3,6)
plot(x_axis, u_ref_center(end,:), 'k-', 'LineWidth', 3, 'DisplayName', 'FDTD Reference')
hold on
plot(x_axis, ud_dirichlet_center(end,:), 'r--', 'LineWidth', 2, 'DisplayName', 'SWR Dirichlet')
plot(x_axis, ud_robin_center(end,:), 'b:', 'LineWidth', 2, 'DisplayName', 'SWR Robin')
xlabel('x')
ylabel('u')
title('Final Time Slice Comparison')
legend('show', 'Location', 'best')
grid on
sgtitle('Time Evolution at Center Line: Method Comparison')

%% Plot detailed convergence comparison
figure()
% subplot(1,2,1)
semilogy(1:length(res_history_dirichlet), res_history_dirichlet, 'o-', 'LineWidth', 2, 'DisplayName', 'Dirichlet')
hold on
semilogy(1:length(res_history_robin), res_history_robin, 's-', 'LineWidth', 2, 'DisplayName', 'Robin')
xlabel('Iteration')
ylabel('Relative Residual')
title('Convergence History')
legend('show')
grid on

% subplot(1,2,2)
% bar([final_res_dirichlet, final_res_robin])
% set(gca, 'XTickLabel', {'Dirichlet', 'Robin'})
% ylabel('Final Relative Residual')
% title('Final Error Comparison')
% grid on
% sgtitle('SWR Convergence Analysis')

%% Display performance summary
fprintf('\n=== PERFORMANCE SUMMARY ===\n');
fprintf('Method            Final Residual    Time (s)    Iterations\n');
fprintf('--------------------------------------------------------\n');
fprintf('SWR Dirichlet     %.3e         %.2f        %d\n', final_res_dirichlet, time_dirichlet, k);
fprintf('SWR Robin         %.3e         %.2f        %d\n', final_res_robin, time_robin, k);

% Calculate accuracy improvement
if final_res_dirichlet > 0 && final_res_robin > 0
    if final_res_dirichlet < final_res_robin
        improvement = final_res_robin / final_res_dirichlet;
        fprintf('Dirichlet is %.2fx more accurate than Robin\n', improvement);
    else
        improvement = final_res_dirichlet / final_res_robin;
        fprintf('Robin is %.2fx more accurate than Dirichlet\n', improvement);
    end
end

%% Plot subdomain boundaries for reference
figure()
x_ticks_arr = zeros(2,N);
aj = @(j) a*(j-1);
bj = @(j) aj(j+1) + M;
for iii = 1:N
    x_ticks_arr(:,iii) = [aj(iii),bj(iii)];
end
x_ticks_arr = reshape(x_ticks_arr,2*N,1);

% Generate labels as 'a1', 'b1', 'a2', 'b2', ...
labels = cell(1, 2*N);
for iii = 1:N
    labels{2*iii-1} = sprintf('a%d', iii);
    labels{2*iii} = sprintf('b%d', iii);
end

% Sort x_ticks_arr and apply the same order to labels
[x_ticks_arr, sortIdx] = sort(x_ticks_arr);
labels = labels(sortIdx);

plot(x_ticks_arr, zeros(size(x_ticks_arr)), 'rx', 'MarkerSize', 10, 'LineWidth', 2)
xlabel('x')
title('Subdomain Interface Locations')
grid on
for i = 1:length(x_ticks_arr)
    text(x_ticks_arr(i), 0.1, labels{i}, 'HorizontalAlignment', 'center')
end
ylim([-0.5, 0.5])