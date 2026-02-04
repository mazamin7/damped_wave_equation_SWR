clear all; close all; clc;
addpath("utils\")

%% VISUALIZATION SETTINGS (POSTER STYLE)
FS_AXIS  = 24;  % Font size for axis ticks
FS_LABEL = 30;  % Font size for X/Y labels
FS_TITLE = 28;  % Font size for titles
LW_BOLD  = 4.0; % Line width for curves
MS_MARK  = 12;  % Marker size for trajectory points
LW_AXIS  = 2.0; % Line width for the axis box
POS_AX_STD = [0.15, 0.15, 0.80, 0.78]; 

% Simulation parameters
P = get_sim_params();
N  = P.N;
a  = P.a;
M  = P.M;
b  = P.b;
Lx = P.Lx;
c  = P.c;

% Infinite Domain approx via Time constraint
dist_to_bound = (1/3) * Lx;
T = 0.8 * (dist_to_bound / c); 
dh = P.dh;
dt = P.dt;
J  = P.J;

% Frequency Band for Asymptotic Analysis
omega_min = pi / T;
omega_max = pi / dt;
fprintf('Frequency Band: w_min = %.2f, w_max = %.2e\n', omega_min, omega_max);

%% Grid Setup
Nx = round(Lx/dh) + 1;
Nt = round(T/dt) + 1;
x_grid = linspace(0,Lx,Nx);
t_grid = linspace(0,T,Nt);

% Initial Condition (Gaussian)
A0    = 1.0;
x_c   = (2/3) * Lx;       
sigma = Lx / 40;          
u0 = @(x) A0 * exp( - (x - x_c).^2 ./ (2*sigma^2) );
v0 = @(x) zeros(size(x)); 

% % --- 8 CASES
% cases = [
%     struct('gamma',1e-4   ,'nu',0)
%     struct('gamma',1e-3   ,'nu',0)
%     struct('gamma',1e-2  ,'nu',0)
%     struct('gamma',1e-1  ,'nu',0)
%     struct('gamma',0   ,'nu',1e-7)
%     struct('gamma',0   ,'nu',1e-6)
%     struct('gamma',0   ,'nu',1e-5)
%     struct('gamma',0   ,'nu',1e-4)
% ];
% nCases = numel(cases);

% --- 8 CASES
cases = [
    struct('gamma',1e4   ,'nu',0)
    struct('gamma',1e5   ,'nu',0)
    struct('gamma',1e6  ,'nu',0)
    struct('gamma',1e7  ,'nu',0)
    struct('gamma',0   ,'nu',1e2)
    struct('gamma',0   ,'nu',1e3)
    struct('gamma',0   ,'nu',1e4)
    struct('gamma',0   ,'nu',1e5)
];
nCases = numel(cases);

% Storage
res_hist_num    = cell(nCases,1); % Numerical Opt history
res_hist_asy    = cell(nCases,1); % Asymptotic Formula history
opt_params_num  = zeros(nCases,2);
opt_params_asy  = zeros(nCases,2);
ref_gt_errors   = zeros(nCases,1); % Reference Error (approx)

x0 = [1.0/c, 0];     % Initial guess for numerical optimizer
global PQ_history
PQ_history = cell(nCases,1);

base_options = optimset('Display','off','TolX',1e-12,'TolFun',1e-12);

% Random initial guess for SWR (shared across cases)
u_init = rand(Nx,Nt);

%% ============================================================
% MAIN LOOP
% ============================================================
for s = 1:nCases
    gamma = cases(s).gamma;
    nu    = cases(s).nu;
    label_s = sprintf('\\gamma=%.3g, \\nu=%.3g', gamma, nu);
    fprintf('\nCase %d: %s\n', s, label_s);
    
    % 1. Reference Solution (FDTD)
    u_ref = run_fdtd(u0, v0, Lx, T, c, dh, dt, gamma, nu);
    
    % 2. Asymptotic Formula p,q
    [p_asy, q_asy, regime_name] = get_formula_pq(gamma, nu, c, omega_min, omega_max);
    opt_params_asy(s,:) = [p_asy, q_asy];
    fprintf('  Formula [%s]: p=%.4f, q=%.4f\n', regime_name, p_asy, q_asy);
    
    % 3. Numerical Optimization p,q
    objfun = @(x) obj_Linf(N,T,dt,J,c,gamma,nu,a,M,x(1),x(2));
    outfun = @(x,optimvalues,state) store_trajectory(x,optimvalues,state,s);
    optim_options = optimset(base_options,'OutputFcn',outfun);
    
    [x_opt, fval] = fminsearch(objfun, x0, optim_options);
    p_num = x_opt(1);
    q_num = x_opt(2);
    opt_params_num(s,:) = x_opt;
    fprintf('  Numeric [Nelder]: p=%.4f, q=%.4f, obj=%.3e\n', p_num, q_num, fval);
    
    % 4. SWR Runs
    % Number of iterations
    if s < 5, k_iter = 500; else, k_iter = 500; end
    
    % --- A. Numerical Best Run ---
    % Test run for scaling
    [~, ~, hist_test_num] = run_swr(u0, v0, N, a, M, T, c, dh, dt, gamma, nu, p_num, q_num, 1, u_init, u_ref);
    err0_num = max(abs(u_init(:) - u_ref(:))) / max(abs(u_ref(:)));
    F_num = hist_test_num(1) / err0_num;
    
    % Real run
    [~, ~, hist_num] = run_swr(u0, v0, N, a, M, T, c, dh, dt, gamma, nu, p_num, q_num, k_iter, u_init/F_num, u_ref);
    res_hist_num{s} = [err0_num; hist_num(:)];
    
    % --- B. Asymptotic Formula Run ---
    % Test run for scaling
    [~, ~, hist_test_asy] = run_swr(u0, v0, N, a, M, T, c, dh, dt, gamma, nu, p_asy, q_asy, 1, u_init, u_ref);
    err0_asy = max(abs(u_init(:) - u_ref(:))) / max(abs(u_ref(:)));
    F_asy = hist_test_asy(1) / err0_asy;
    
    % Real run
    [~, ~, hist_asy] = run_swr(u0, v0, N, a, M, T, c, dh, dt, gamma, nu, p_asy, q_asy, k_iter, u_init/F_asy, u_ref);
    res_hist_asy{s} = [err0_asy; hist_asy(:)];
end

%% ============================================================
% PLOTTING ROUTINES
% ============================================================
% --- Helper to plot Num vs Asymp ---
% Solid = Num, Dashed = Asymp
plot_convergence_group(cases, res_hist_num, res_hist_asy, 'gamma', 'Gamma Convergence', POS_AX_STD, FS_AXIS, FS_LABEL, LW_BOLD, LW_AXIS);
plot_convergence_group(cases, res_hist_num, res_hist_asy, 'nu',    'Nu Convergence',    POS_AX_STD, FS_AXIS, FS_LABEL, LW_BOLD, LW_AXIS);

% ============================================================
% TRAJECTORY FIGURES (with Formula Point)
% ============================================================
for s = 1:nCases
    pq = PQ_history{s};
    if isempty(pq), continue; end
    gamma = cases(s).gamma;
    nu    = cases(s).nu;
    
    % Current Formula Point
    p_f = opt_params_asy(s,1);
    q_f = opt_params_asy(s,2);
    
    % Current Num Point
    p_n = opt_params_num(s,1);
    q_n = opt_params_num(s,2);
    
    % Objective Landscape
    objfun_pq = @(p,q) obj_Linf(N,T,dt,J,c,gamma,nu,a,M,p,q);
    
    % Define plot bounds to include both trajectory and formula point
    p_all = [pq(:,1); p_f];
    q_all = [pq(:,2); q_f];
    pmin = min(p_all); pmax = max(p_all);
    qmin = min(q_all); qmax = max(q_all);
    
    ppad = 0.3*max(1, pmax-pmin);
    qpad = 0.3*max(1, qmax-qmin);
    
    pr = linspace(pmin-ppad, pmax+ppad, 50);
    qr = linspace(qmin-qpad, qmax+qpad, 50);
    [P,Q] = meshgrid(pr,qr);
    
    % Evaluate landscape (can be slow, but robust)
    Z = zeros(size(P));
    for i=1:numel(P)
        Z(i) = objfun_pq(P(i), Q(i));
    end
    
    figure('Name', sprintf('Traj Case %d', s), 'Color', 'w'); clf; hold on;
    contourf(P, Q, log10(Z), 30, 'LineStyle','none'); 
    
    cb = colorbar;
    cb.FontSize = FS_AXIS;
    cb.Label.String = 'log_{10}(Error)';
    
    % Plot Optimization Path
    plot(pq(:,1), pq(:,2), '-o', 'LineWidth', 2.0, 'Color','k', 'MarkerSize', 4);
    
    % Start (White Square)
    plot(pq(1,1), pq(1,2), 'ws', 'MarkerSize', 10, 'LineWidth', 2, 'MarkerFaceColor', 'k');
    
    % Numerical Optimum (Red Star)
    plot(p_n, q_n, 'p', 'MarkerSize', 18, 'LineWidth', 2, 'MarkerFaceColor', 'r', 'MarkerEdgeColor','k');
    
    % Asymptotic Formula (Green Hexagram)
    plot(p_f, q_f, 'h', 'MarkerSize', 18, 'LineWidth', 2, 'MarkerFaceColor', 'g', 'MarkerEdgeColor','k');
    
    xlabel('p', 'FontSize', FS_LABEL, 'FontWeight', 'bold'); 
    ylabel('q', 'FontSize', FS_LABEL, 'FontWeight', 'bold');
    
    title_str = sprintf('\\gamma=%.3g, \\nu=%.3g', gamma, nu);
    title(title_str, 'FontSize', FS_TITLE, 'FontWeight', 'bold');
    
    legend({'Landscape','Path','Start','Num. Opt','Asymp. Formula'}, 'Location','best', 'FontSize', 14);
    grid on; 
    set(gca,'FontSize',FS_AXIS, 'LineWidth', LW_AXIS, 'FontWeight', 'bold');
end

%% ============================================================
% LOCAL FUNCTIONS
% ============================================================
function [p, q, regime] = get_formula_pq(gamma, nu, c, w_min, w_max, delta)
    % GET_FORMULA_PQ Computes asymptotically optimal transmission parameters.
    %
    % Inputs:
    %   gamma, nu : Damping parameters
    %   c         : Wave speed
    %   w_min, max: Frequency band
    %   delta     : (Optional) Overlap size. Defaults to 0.
    %
    % Outputs:
    %   p, q      : Transmission coefficients (Lambda = q + s*p)
    %   regime    : String description of the active regime
    
    if nargin < 6, delta = 0; end

    % Common Geometric Mean vars (for Small Delta / Equioscillation)
    w_mid = sqrt(w_min * w_max);
    z     = (w_min / w_max)^(0.25);
    Y     = sqrt(z / (1 + z + z^2));
    
    if nu > 1e-9
        % =========================================================
        % VISCOELASTIC REGIME
        % =========================================================
        % Check: Viscoelastic Diffusive if nu*w_min/c^2 >> 1
        is_diffusive = (nu * w_min / c^2) > 1.0;
        
        if ~is_diffusive
            % --- Visco Propagative ---
            regime = 'Visco Prop';
            
            % Check Overlap Condition: nu * w_min^2 * delta / c^3
            % (Transition to Big Overlap logic)
            is_big_overlap = (nu * w_min^2 * delta / (2*c^3)) > 1.0;
            
            if is_big_overlap
                % Large Overlap: Match w_min exactly
                q = (nu * w_min^2) / (2 * c^3);
                p = 1/c - (3 * nu^2 * w_min^2) / (8 * c^5);
                regime = [regime ' (Big \delta)'];
            else
                % Small Overlap: Equioscillation (Minimax)
                q = (nu * w_min * w_max) / (2 * c^3);
                p = 1/c - (nu^2 / (4*c^5)) * (w_max^2 + 0.5*w_min^2); 
            end
            
        else
            % --- Visco Diffusive ---
            regime = 'Visco Diff';
            
            % Check Overlap Condition: delta * sqrt(w_min / nu)
            is_big_overlap = (delta * sqrt(w_min / (2*nu))) > 1.0;
            
            if is_big_overlap
                % Large Overlap: Pointwise match at w_min
                q = sqrt(w_min / (2*nu));
                p = 1 / sqrt(2 * nu * w_min);
                regime = [regime ' (Big \delta)'];
            else
                % Small Overlap: 3-Point Equioscillation (w_mid)
                p = Y / sqrt(2 * nu * w_mid);
                q = p * w_mid;
            end
        end
        
    elseif gamma > 1e-9
        % =========================================================
        % TELEGRAPHER REGIME
        % =========================================================
        % Check: Telegrapher Diffusive if gamma >> w_max
        % (Conservative threshold: if gamma dominates even the highest freq)
        is_diffusive = (gamma / w_max) > 1.0;
        
        if ~is_diffusive
            % --- Tele Propagative ---
            regime = 'Tele Prop';
            
            % Check Overlap Condition: gamma * delta / c
            is_big_overlap = (gamma * delta / (2*c)) > 1.0;
            
            if is_big_overlap
                % Large Overlap: Match w_min
                q = gamma / (2 * c);
                p = 1/c + (gamma^2 / (8 * c * w_min^2));
                regime = [regime ' (Big \delta)'];
            else
                % Small Overlap: Equioscillation
                q = gamma / (2 * c);
                p = 1/c + (gamma^2 / (16*c)) * (1/w_min^2 + 1/w_max^2);
            end
            
        else
            % --- Tele Diffusive ---
            regime = 'Tele Diff';
            
            % Check Overlap Condition: delta/c * sqrt(gamma * w_min)
            is_big_overlap = (delta/c * sqrt(gamma * w_min / 2)) > 1.0;
            
            if is_big_overlap
                % Large Overlap: Pointwise match at w_min
                q = sqrt(gamma * w_min) / (c * sqrt(2));
                p = sqrt(gamma) / (c * sqrt(2 * w_min));
                regime = [regime ' (Big \delta)'];
            else
                % Small Overlap: 3-Point Equioscillation
                % Isomorphic to Visco Diff, scaled by sqrt(gamma)/c
                factor = sqrt(gamma) / c;
                p = (factor / sqrt(2 * w_mid)) * Y;
                q = p * w_mid;
            end
        end
    else
        % =========================================================
        % PURE WAVE EQUATION
        % =========================================================
        regime = 'Pure Wave';
        p = 1/c; 
        q = 0;
    end
end

function plot_convergence_group(cases, hist_num, hist_asy, type, fig_name, pos, fs_ax, fs_lbl, lw_bold, lw_ax)
    figure('Name', fig_name, 'Color', 'w'); clf;
    set(gcf, 'Position', [100 100 900 600]); 
    
    if strcmp(type, 'gamma')
        idx_list = find([cases.nu] == 0);
    else
        idx_list = find([cases.gamma] == 0);
    end
    
    h_plots = [];
    labels  = {};
    
    % Removed 'hold on' here to prevent locking linear scale before semilogy
    colors = lines(numel(idx_list));
    
    for j = 1:numel(idx_list)
        s = idx_list(j);
        col = colors(j,:);
        
        % Numerical (Solid)
        res = hist_num{s};
        its = 0:(numel(res)-1);
        hp = semilogy(its, res, '-', 'LineWidth', lw_bold, 'Color', col);
        
        hold on; % Call hold on AFTER the first semilogy
        
        h_plots(end+1) = hp;
        
        % Asymptotic (Dashed)
        res_a = hist_asy{s};
        its_a = 0:(numel(res_a)-1);
        semilogy(its_a, res_a, '--', 'LineWidth', lw_bold, 'Color', col, 'HandleVisibility','off');
        
        if strcmp(type, 'gamma')
            lbl = sprintf('\\gamma=%.3g', cases(s).gamma);
        else
            lbl = sprintf('\\nu=%.3g', cases(s).nu);
        end
        labels{end+1} = lbl;
    end
    
    xlabel('Iteration k', 'FontSize', fs_lbl, 'FontWeight', 'bold');
    ylabel('Error', 'FontSize', fs_lbl, 'FontWeight', 'bold');
    grid on;
    set(gca, 'FontSize', fs_ax, 'LineWidth', lw_ax, 'FontWeight', 'bold', 'Position', pos);
    
    % Add dummy lines for legend explanation
    h_solid = plot(nan,nan, 'k-', 'LineWidth', 2);
    h_dash  = plot(nan,nan, 'k--', 'LineWidth', 2);
    
    legend([h_plots, h_solid, h_dash], [labels, {'Num. Opt.', 'Formula'}], ...
           'Location', 'SouthWest', 'FontSize', 16);
end

function stop = store_trajectory(x,optimvalues,state,idx)
    stop = false;
    global PQ_history
    if strcmp(state,'init')
        PQ_history{idx} = x(:).';
    elseif strcmp(state,'iter')
        PQ_history{idx}(end+1,:) = x(:).';
    end
end