% WEAK_SCALING_SWR_2D_OPT - Weak scalability with optimal Robin (p,q)
clear; close all; clc;

addpath('utils\')

% Simulation parameters
P = get_sim_params_1D();

N  = P.N; % ignore
a  = P.a;
M  = P.M;
b  = P.b;
Lx = P.Lx; % ignore
T  = P.T;

c  = P.c;
gamma = P.gamma;
nu = P.nu;

dh = P.dh;
dt = P.dt;
J  = P.J;

%%
% SWR iterations per subdomain
% k_per_dom = 50;        % viscous
k_per_dom = 50;      % viscoelastic

% test sizes
N_list = [2 4 8 16];

% initial conditions
u0_fun = @(x,y) 0.*x;

% storage
final_res = zeros(size(N_list));
time_robin = zeros(size(N_list));
Nx_local = zeros(size(N_list));
Nx_global = zeros(size(N_list));
p_opt = zeros(numel(N_list),1);
q_opt = zeros(numel(N_list),1);

% ---- optimizer setup for p,q (Robin Λ = p ∂t + q) ----
tol = 1e-6;
options = optimset('Display','off','TolX',tol,'TolFun',tol);
x0 = [1/c, 0];                 % warm start for (p,q)
ky = 0;                        % optimize along normal incidence
objfun = @(N,T,dt,J,c,gamma,nu,a,M,p,q,ky) ...
           obj_Linf(N,T,dt,J,c,gamma,nu,a,M,p,q,ky);   % user-provided

for ii = 1:numel(N_list)
    N = N_list(ii);

    % geometry
    aj  = @(j) a*(j-1);
    bj  = @(j) aj(j+1) + M;
    Lx  = bj(N);

    % grids
    Nx = round(Lx/dh) + 1;  % Ny = round(Ly/dh) + 1;
    x_axis = linspace(0,Lx,Nx);
    % y_axis = linspace(0,Ly,Ny);

    % ICs depend on domain size (Gaussian velocity bump near (Lx/3,Ly/3))
    v0_fun = @(x) exp(-sqrt((x-Lx/3).^2));

    % reference and initial guess
    fprintf('[N=%d] building reference...\n',N);
    tic;
    u_ref  = run_fdtd_1D(u0_fun, v0_fun, Lx, T, c, dh, dt, gamma, nu);
    tref = toc; %#ok<NASGU>
    u_init = rand(size(u_ref));

    % optimize (p,q) for this N
    Nt = size(u_ref,3);     % J in your notation
    J  = Nt;
    fprintf('[N=%d] optimizing (p,q)...\n',N);
    f = @(x) objfun(N,T,dt,J,c,gamma,nu,a,M,x(1),x(2),ky);
    [x_opt, ~] = fminsearch(f, x0, options);
    p = x_opt(1);  q = x_opt(2);
    % p = 1/c; q = 0;
    p_opt(ii) = p; q_opt(ii) = q;
    x0 = x_opt;    % warm-start next N

    % SWR iterations
    k = k_per_dom * N;

    % run SWR with optimized (p,q)
    fprintf('[N=%d] SWR Robin with p=%.3g, q=%.3g...\n',N,p,q);
    tic;
    [~, final_res(ii), ~] = run_swr_1D(u0_fun, v0_fun, N, a, M, T, c, dh, dt, gamma, nu, p, q, k, u_init, u_ref);
    time_robin(ii) = toc;

    % weak-scaling bookkeeping
    ajd = @(j) round(aj(j)/dh) + 1;
    bjd = @(j) round(bj(j)/dh) + 1;
    Nx_local(ii)  = bjd(1) - ajd(1) + 1;
    Nx_global(ii) = Nx;
end

% ----- plots -----
figure; 
loglog(N_list, final_res, 'o-','LineWidth',2);
grid on; xlabel('N'); ylabel('Final SWR residual');
title('Weak scalability: final SWR error vs N (optimized Robin)');

figure;
plot(N_list, Nx_local, 's-','LineWidth',2); grid on;
xlabel('N'); ylabel('Nx per subdomain'); title('Weak scaling check');

figure;
loglog(N_list, time_robin, 'd-','LineWidth',2); grid on;
xlabel('N'); ylabel('Time (s)'); title('SWR wall-clock vs N');

% report
fprintf('\nN   Lx      Nx(local) Nx(global)   p_opt        q_opt        final_res     time(s)\n');
for ii = 1:numel(N_list)
    fprintf('%-3d %-7.3f %-9d %-11d %-12.4e %-12.4e %.3e   %.2f\n', ...
        N_list(ii), a*N_list(ii)+M, Nx_local(ii), Nx_global(ii), ...
        p_opt(ii), q_opt(ii), final_res(ii), time_robin(ii));
end
