clear all; close all; clc;
addpath("utils\")

% Simulation parameters
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
Nx = round(Lx/dh) + 1;
Nt = round(T/dt) + 1;
x_grid = linspace(0,Lx,Nx);
t_grid = linspace(0,T,Nt);

% --- Modal initial condition
m_mode  = 1;
k0      = m_mode*pi/Lx;
A0      = 1.0;      % Initial displacement amplitude
v0amp   = 0.0;      % Initial velocity amplitude

u0  = @(x) A0 * sin(k0*x);
v0  = @(x) v0amp * sin(k0*x);

% --- 8 CASES
cases = [
    struct('gamma',4   ,'nu',0)
    struct('gamma',8   ,'nu',0)
    struct('gamma',10  ,'nu',0)
    struct('gamma',12  ,'nu',0)
    struct('gamma',0   ,'nu',0.001)
    struct('gamma',0   ,'nu',0.01)
    struct('gamma',0   ,'nu',0.05)
    struct('gamma',0   ,'nu',0.1)
];

nCases = numel(cases);

% Storage
res_hist        = cell(nCases,1);
final_errors    = zeros(nCases,1);
opt_params      = zeros(nCases,2);
opt_obj_vals    = zeros(nCases,1);
ref_gt_errors   = zeros(nCases,1);   % FDTD vs GT (final time) error

x0 = [1.0/c, 0];     % Initial guess (p,q)
ky = 0;

global PQ_history
PQ_history = cell(nCases,1);

base_options = optimset('Display','off','TolX',1e-8,'TolFun',1e-8);

% Random initial guess for SWR (shared across cases)
u_init = rand(Nx,Nt);

% ============================================================
% MAIN LOOP
% ============================================================
for s = 1:nCases

    gamma = cases(s).gamma;
    nu    = cases(s).nu;
    label_s = sprintf('\\gamma=%.3g, \\nu=%.3g', gamma, nu);

    fprintf('\nCase %d: %s\n', s, label_s);

    % Ground-truth analytical solution
    u_gt = analytic_solution_single_mode(x_grid, t_grid, c, gamma, nu, ...
                                         k0, A0, v0amp);

    % FDTD reference solution (as in original code)
    u_ref = run_fdtd_1D(u0, v0, Lx, T, c, dh, dt, gamma, nu);

    % Final-time relative error between FDTD and GT (L-infinity in space)
    ref_gt_errors(s) = max(abs(u_ref(:,end) - u_gt(:,end))) / max(abs(u_gt(:,end)));

    % Objective (assumed to use gamma, nu, etc., as in your current obj_Linf)
    objfun = @(x) obj_Linf(N,T,dt,J,c,gamma,nu,a,M,x(1),x(2),ky);

    % Save trajectory
    outfun = @(x,optimvalues,state) store_trajectory(x,optimvalues,state,s);
    optim_options = optimset(base_options,'OutputFcn',outfun);

    % Optimize p,q
    [x_opt,fval] = fminsearch(objfun, x0, optim_options);
    p_opt = x_opt(1);
    q_opt = x_opt(2);

    opt_params(s,:) = x_opt;
    opt_obj_vals(s) = fval;

    fprintf('  optimized p,q = (%.6f, %.6f), obj = %.3e\n', p_opt, q_opt, fval);

    % Number of SWR iterations for the "real" run
    if s < 5
        k = 60;
    else
        k = 30;
    end

    % ========================================================
    % TEST RUN: 1 SWR iteration to get amplification factor F
    % ========================================================
    k_test = 1;

    % Test run starting from common random u_init
    [~, ~, res_history_test] = run_swr_1D( ...
        u0, v0, N, a, M, T, c, dh, dt, gamma, nu, p_opt, q_opt, ...
        k_test, u_init, u_ref);

    % Error at iteration 0 for the test run
    err0_test = max(abs(u_init(:) - u_ref(:))) / max(abs(u_ref(:)));

    % Error after first SWR iteration in the test run
    E1_test = res_history_test(1);

    % Amplification factor F = E1 / E0
    F = E1_test / err0_test;

    % ========================================================
    % REAL RUN: rescaled initial condition u_init / F
    % ========================================================
    u_init_scaled = u_init / F;

    [~, final_res, res_history] = run_swr_1D( ...
        u0, v0, N, a, M, T, c, dh, dt, gamma, nu, p_opt, q_opt, ...
        k, u_init_scaled, u_ref);

    % Initial error for the scaled run (relative Linf in space-time)
    err0 = err0_test;

    % Prepend iteration-0 error to history
    res_hist{s}     = [err0; res_history(:)];
    final_errors(s) = final_res;
end

% ============================================================
% FIGURE 1: gamma (nu=0) only optimized curves
%          + horizontal line = FDTD vs GT final-time error
%          conservative color handling, including iteration 0
% ============================================================
figure(1); clf;
gamma_idx = find([cases.nu] == 0);

h_curves = gobjects(numel(gamma_idx),1);

for j = 1:numel(gamma_idx)
    s      = gamma_idx(j);
    resvec = res_hist{s};
    nIter  = numel(resvec) - 1;     % iterations 0..nIter
    its    = 0:nIter;

    % SWR error history: let MATLAB choose the color
    h = semilogy(its, resvec, 'LineWidth', 2);
    hold on;
    col = get(h,'Color');           % color actually used

    % Store handle for legend
    h_curves(j) = h;

    % Horizontal line: FDTD vs GT final-time error (same color, dashed)
    semilogy([0 nIter], ref_gt_errors(s)*[1 1], ...
             'LineWidth', 1.5, ...
             'LineStyle', '--', ...
             'Color', col, ...
             'HandleVisibility','off');
end

xlabel('Iteration','FontSize',16);
ylabel('Error','FontSize',16);
grid on; set(gca,'FontSize',18);

% ylim([1e-10, 1e10])

% Automatic labels: only gamma (nu=0)
gamma_labels = arrayfun(@(c) sprintf('\\gamma=%.3g', c.gamma), cases(gamma_idx), ...
                        'UniformOutput', false);
legend(h_curves, gamma_labels,'Location','SouthWest');

% ============================================================
% FIGURE 2: nu (gamma=0) only optimized curves
%          + horizontal line = FDTD vs GT final-time error
%          conservative color handling, including iteration 0
% ============================================================
figure(2); clf;
nu_idx = find([cases.gamma] == 0);

h_curves_nu = gobjects(numel(nu_idx),1);

for j = 1:numel(nu_idx)
    s      = nu_idx(j);
    resvec = res_hist{s};
    nIter  = numel(resvec) - 1;     % iterations 0..nIter
    its    = 0:nIter;

    % SWR error history: let MATLAB choose the color
    h = semilogy(its, resvec, 'LineWidth', 2);
    hold on;
    col = get(h,'Color');

    % Store handle for legend
    h_curves_nu(j) = h;

    % Horizontal line: FDTD vs GT final-time error (same color, dashed)
    semilogy([0 nIter], ref_gt_errors(s)*[1 1], ...
             'LineWidth', 1.5, ...
             'LineStyle', '--', ...
             'Color', col, ...
             'HandleVisibility','off');
end

xlabel('Iteration','FontSize',16);
ylabel('Error','FontSize',16);
grid on; set(gca,'FontSize',18);

% Automatic labels: only nu (gamma=0)
nu_labels = arrayfun(@(c) sprintf('\\nu=%.3g', c.nu), cases(nu_idx), ...
                     'UniformOutput', false);
legend(h_curves_nu, nu_labels,'Location','NorthEast');

% ============================================================
% TRAJECTORY FIGURES
% ============================================================
for s = 1:nCases
    pq = PQ_history{s};
    if isempty(pq), continue; end

    gamma = cases(s).gamma;
    nu    = cases(s).nu;

    objfun_pq = @(p,q) obj_Linf(N,T,dt,J,c,gamma,nu,a,M,p,q,ky);

    pmin = min(pq(:,1)); pmax = max(pq(:,1));
    qmin = min(pq(:,2)); qmax = max(pq(:,2));

    ppad = 0.2*max(1, pmax-pmin);
    qpad = 0.2*max(1, qmax-qmin);

    pr = linspace(pmin-ppad,pmax+ppad,40);
    qr = linspace(qmin-qpad,qmax+qpad,40);
    [P,Q] = meshgrid(pr,qr);

    Z = arrayfun(@(pp,qq) objfun_pq(pp,qq), P, Q);

    figure(2+s); clf; hold on;
    contourf(P, Q, log10(Z), 20, 'LineStyle','none'); colorbar;

    plot(pq(:,1), pq(:,2), '-o', 'LineWidth',2, 'Color','k');
    plot(pq(1,1),  pq(1,2),  'ws', 'MarkerSize',10,'LineWidth',2);
    plot(pq(end,1),pq(end,2),'rx', 'MarkerSize',12,'LineWidth',2);

    xlabel('p'); ylabel('q');
    title(sprintf('\\gamma=%.3g, \\nu=%.3g — fminsearch trajectory', gamma, nu));
    grid on; set(gca,'FontSize',18);
end

% ============================================================
% SUMMARY
% ============================================================
fprintf('\n=== Final errors vs GT (optimized only) ===\n');
for s = 1:nCases
    gamma = cases(s).gamma;
    nu    = cases(s).nu;
    label_s = sprintf('gamma=%.3g, nu=%.3g', gamma, nu);

    fprintf('%-25s  error_SWR = %.3e   (p,q) = (%.5f, %.5f)   error_FDTD_vs_GT_final = %.3e\n', ...
        label_s, final_errors(s), opt_params(s,1), opt_params(s,2), ref_gt_errors(s));
end

% ============================================================
% Store trajectory
% ============================================================
function stop = store_trajectory(x,optimvalues,state,idx)
    stop = false;
    global PQ_history
    if strcmp(state,'init')
        PQ_history{idx} = x(:).';
    elseif strcmp(state,'iter')
        PQ_history{idx}(end+1,:) = x(:).';
    end
end
