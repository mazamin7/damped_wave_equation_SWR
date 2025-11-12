clear all; close all; clc;
addpath("utils\")

% --- Geometry / discretization
N  = 2;
a  = 0.3;
M  = 0.1;
Lx = N*a + M;

c  = 1;  dh = 0.01;  dt = 0.01;
T  = 5;  k  = 40*N;

% --- Initial condition (same as before)
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
u0 = @(x) gaussian_normalized(x) + sine_sum_normalized(x);
v0 = @(x) 0.*x;

Nx = round(Lx/dh) + 1;
Nt = round(T/dt) + 1;

% --- Only L_inf-optimized cases
cases = [
    struct('label','\gamma=1,\nu=0','gamma',1,'nu',0,'theta1',1.0 ,'theta2',0.5 )
    struct('label','\gamma=5,\nu=0','gamma',5,'nu',0,'theta1',1.05,'theta2',2.5 )
    struct('label','\gamma=0,\nu=0.1','gamma',0,'nu',0.1,'theta1',0.5 ,'theta2',2.0 )
    struct('label','\gamma=0,\nu=0.5','gamma',0,'nu',0.5,'theta1',0.15,'theta2',3.5 )
];

fprintf('T*c/M = %.3f\n', T*c/M);

res_hist = cell(4,1);
final_errors = zeros(4,1);

% --- Run each case
for s = 1:4
    gamma = cases(s).gamma; nu = cases(s).nu;
    theta1 = cases(s).theta1; theta2 = cases(s).theta2;
    fprintf('\nCase %d: %s\n', s, cases(s).label);
    u_ref = run_fdtd_1D(u0,v0,Lx,T,c,dh,dt,gamma,nu);
    u_init = rand(Nx,Nt);
    [~,final_res,res_history] = run_swr_1D( ...
        u0,v0,N,a,M,T,c,dh,dt,gamma,nu,theta1,theta2,k,u_init,u_ref);
    res_hist{s} = res_history;
    final_errors(s) = final_res;
end

% --- Figure 1: gamma variations (nu=0)
figure(1);
colors = ['k','r']; styles = {'-','--'};
for s = 1:2
    semilogy(1:length(res_hist{s}), res_hist{s}, ...
        [colors(s) styles{s}], 'LineWidth',2); hold on;
end
xlabel('Iteration','FontSize',16); ylabel('Error','FontSize',16);
xlim([0,k]); ylim([1e-15,1e5]);
legend({cases(1:2).label},'Location','NorthEast','FontSize',14);
grid on; set(gca,'FontSize',18);

% --- Figure 2: nu variations (gamma=0)
figure(2);
colors = ['m','g']; styles = {':','-.'};
for s = 3:4
    semilogy(1:length(res_hist{s}), res_hist{s}, ...
        [colors(s-2) styles{s-2}], 'LineWidth',2); hold on;
end
xlabel('Iteration','FontSize',16); ylabel('Error','FontSize',16);
xlim([0,k]); ylim([1e-15,1e5]);
legend({cases(3:4).label},'Location','NorthEast','FontSize',14);
grid on; set(gca,'FontSize',18);

% --- Summary
fprintf('\n=== Final errors ===\n');
for s = 1:4
    fprintf('%-20s %.3e\n', cases(s).label, final_errors(s));
end
