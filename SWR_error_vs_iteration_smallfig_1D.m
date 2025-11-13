clear all; close all; clc;
addpath("utils\")

% --- Geometry / discretization
N  = 2;
a  = 0.3;
M  = 0.1;
Lx = N*a + M;

c  = 1;  dh = 0.01;  dt = 0.01;
T  = 5;  k  = 40*N;
J  = 50;                     % number of grid points in one subdomain (assumption)

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

% --- Cases: old (p,q) included
cases = [
    struct('label','\gamma=1,\nu=0'  ,'gamma',1,'nu',0  ,'theta1_old',1.0 ,'theta2_old',0.5 )
    struct('label','\gamma=5,\nu=0'  ,'gamma',5,'nu',0  ,'theta1_old',1.05,'theta2_old',2.5 )
    struct('label','\gamma=0,\nu=0.1','gamma',0,'nu',0.1,'theta1_old',0.5 ,'theta2_old',2.0 )
    struct('label','\gamma=0,\nu=0.5','gamma',0,'nu',0.5,'theta1_old',0.15,'theta2_old',3.5 )
];

fprintf('T*c/M = %.3f\n', T*c/M);

res_hist          = cell(4,1);   % optimized parameters
res_hist_old      = cell(4,1);   % old parameters (still computed for diagnostics)
final_errors      = zeros(4,1);  % optimized
final_errors_old  = zeros(4,1);  % old
opt_params        = zeros(4,2);  % [theta1, theta2] for each case
opt_obj_vals      = zeros(4,1);  % optimal objective values
old_obj_vals      = zeros(4,1);  % objective at old parameters

% --- fminsearch base options and initial guess for [theta1, theta2]
x0 = [1.0, 1.0];  % initial guess for [theta1, theta2] (p,q)
base_options = optimset('Display','off', ...
                        'TolX',1e-4, ...
                        'TolFun',1e-8);

ky = 0;   % IMPORTANT: here ky is forced to zero

% Storage for trajectories (p,q) across iterations
global PQ_history
PQ_history = cell(4,1);

% --- Run each case
for s = 1:4
    gamma = cases(s).gamma;
    nu    = cases(s).nu;

    fprintf('\nCase %d: %s\n', s, cases(s).label);

    % Reference solution
    u_ref = run_fdtd_1D(u0, v0, Lx, T, c, dh, dt, gamma, nu);

    % Objective for this (gamma, nu)
    % obj_Linf(N, T, dt, J, c, gamma, nu, a, M, theta1, theta2, ky)
    objfun = @(x) obj_Linf(N, T, dt, J, c, gamma, nu, a, M, x(1), x(2), ky);

    % OLD parameters and objective
    theta1_old = cases(s).theta1_old;
    theta2_old = cases(s).theta2_old;
    old_obj = objfun([theta1_old, theta2_old]);
    old_obj_vals(s) = old_obj;

    % Output function to record trajectory in (p,q) plane
    outfun = @(x,optimvalues,state) store_trajectory(x,optimvalues,state,s);
    optim_options = optimset(base_options, 'OutputFcn', outfun);

    % Use the same random initial guess for OLD and OPTIMIZED runs
    u_init = rand(Nx,Nt);

    % SWR with OLD parameters (for diagnostics, not plotted in Fig 1–2)
    [~, final_res_old, res_history_old] = run_swr_1D( ...
        u0, v0, N, a, M, T, c, dh, dt, gamma, nu, theta1_old, theta2_old, k, u_init, u_ref);
    res_hist_old{s}     = res_history_old;
    final_errors_old(s) = final_res_old;

    % Optimize [theta1, theta2]
    [x_opt, fval] = fminsearch(objfun, x0, optim_options);
    theta1 = x_opt(1);
    theta2 = x_opt(2);
    opt_params(s,:)   = x_opt;
    opt_obj_vals(s)   = fval;

    fprintf('  -> optimized (theta1, theta2) = (%.6f, %.6f), obj = %.3e\n', ...
            theta1, theta2, fval);

    % SWR with optimized parameters
    [~, final_res, res_history] = run_swr_1D( ...
        u0, v0, N, a, M, T, c, dh, dt, gamma, nu, theta1, theta2, k, u_init, u_ref);
    res_hist{s}      = res_history;
    final_errors(s)  = final_res;
end

% --- Figure 1: gamma variations (nu=0) – ONLY NEW (optimized) curves
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

% --- Figure 2: nu variations (gamma=0) – ONLY NEW (optimized) curves
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

% --- Figures 3–6: trajectories + objective contour in (p,q) plane
for s = 1:4
    pq = PQ_history{s};
    if isempty(pq), continue; end

    gamma = cases(s).gamma;
    nu    = cases(s).nu;

    % Objective as function of p,q for this case
    objfun_pq = @(p,q) obj_Linf(N, T, dt, J, c, gamma, nu, a, M, p, q, ky);

    % Region in (p,q) plane to visualize: around trajectory + old + optimized
    p_all = [pq(:,1); cases(s).theta1_old; opt_params(s,1)];
    q_all = [pq(:,2); cases(s).theta2_old; opt_params(s,2)];
    pmin  = min(p_all); pmax = max(p_all);
    qmin  = min(q_all); qmax = max(q_all);

    % Add padding
    if pmax > pmin
        ppad = 0.2*(pmax - pmin);
    else
        ppad = 0.2*max(1,abs(pmin)); % fallback
    end
    if qmax > qmin
        qpad = 0.2*(qmax - qmin);
    else
        qpad = 0.2*max(1,abs(qmin)); % fallback
    end

    pr = linspace(pmin - ppad, pmax + ppad, 40);
    qr = linspace(qmin - qpad, qmax + qpad, 40);
    [P,Q] = meshgrid(pr, qr);
    Z = zeros(size(P));

    % Evaluate objective on grid
    for i = 1:numel(P)
        Z(i) = objfun_pq(P(i), Q(i));
    end

    figure(2 + s);  % Figures 3,4,5,6

    % Background contour of obj (log10 scale)
    contourf(P, Q, log10(Z), 20, 'LineStyle','none'); hold on;
    colorbar;
    colormap parula;

    % Trajectory and points
    plot(pq(:,1), pq(:,2), '-o', 'LineWidth', 2, 'Color','k'); hold on;
    plot(pq(1,1),  pq(1,2),  's', 'MarkerSize',10, 'LineWidth',2, 'Color','w'); % start
    plot(pq(end,1),pq(end,2),'x', 'MarkerSize',10, 'LineWidth',2, 'Color','r'); % end

    % Old predefined (reference)
    plot(cases(s).theta1_old, cases(s).theta2_old, 'd', ...
        'MarkerSize',10, 'LineWidth',2, 'Color','m');

    xlabel('p = \theta_1','FontSize',16);
    ylabel('q = \theta_2','FontSize',16);
    title(sprintf('fminsearch trajectory and obj contour: %s', cases(s).label), 'FontSize',16);
    grid on; set(gca,'FontSize',18);
    legend({'objective contour','trajectory','start','end','old predefined'}, ...
        'Location','best');
end

% --- Summary
fprintf('\n=== Final errors (SWR with optimized thetas) ===\n');
for s = 1:4
    fprintf('%-20s error_new = %.3e   error_old = %.3e   (theta1_new, theta2_new) = (%.6f, %.6f)\n', ...
        cases(s).label, final_errors(s), final_errors_old(s), ...
        opt_params(s,1), opt_params(s,2));
end

fprintf('\n=== Final comparison (old thetas vs optimized) ===\n');
for s = 1:4
    fprintf('%-20s  old_obj = %.3e   new_obj = %.3e   (old_p, old_q) = (%.3f, %.3f)   (new_p, new_q) = (%.3f, %.3f)\n', ...
        cases(s).label, ...
        old_obj_vals(s), ...
        opt_obj_vals(s), ...
        cases(s).theta1_old, cases(s).theta2_old, ...
        opt_params(s,1), opt_params(s,2));
end

% ============================================================
% Local function to store (p,q) trajectory from fminsearch
% ============================================================
function stop = store_trajectory(x,optimvalues,state,idx)
    stop = false;
    global PQ_history
    switch state
        case 'init'
            PQ_history{idx} = x(:).';
        case 'iter'
            PQ_history{idx}(end+1,:) = x(:).';
        case 'done'
            % nothing special
    end
end
