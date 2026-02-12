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
sigma = Lx / 160;          
u0 = @(x) A0 * exp( - (x - x_c).^2 ./ (2*sigma^2) );
v0 = @(x) zeros(size(x)); 

% % --- 8 CASES
% cases = [
%     struct('gamma',1e-5   ,'nu',0)
%     struct('gamma',1e-4   ,'nu',0)
%     struct('gamma',1e-3  ,'nu',0)
%     struct('gamma',1e-2  ,'nu',0)
%     struct('gamma',0   ,'nu',1e-7)
%     struct('gamma',0   ,'nu',1e-6)
%     struct('gamma',0   ,'nu',1e-5)
%     struct('gamma',0   ,'nu',1e-4)
% ];
% nCases = numel(cases);

% --- 8 CASES
cases = [
    struct('gamma',1e0   ,'nu',0)
    struct('gamma',1e1   ,'nu',0)
    struct('gamma',1e2  ,'nu',0)
    struct('gamma',1e3  ,'nu',0)
    struct('gamma',0   ,'nu',1e-2)
    struct('gamma',0   ,'nu',1e-1)
    struct('gamma',0   ,'nu',1e0)
    struct('gamma',0   ,'nu',1e1)
];
nCases = numel(cases);

% % --- 8 CASES
% cases = [
%     struct('gamma',1e4   ,'nu',0)
%     struct('gamma',1e5   ,'nu',0)
%     struct('gamma',1e6  ,'nu',0)
%     struct('gamma',1e7  ,'nu',0)
%     struct('gamma',0   ,'nu',1e2)
%     struct('gamma',0   ,'nu',1e3)
%     struct('gamma',0   ,'nu',1e4)
%     struct('gamma',0   ,'nu',1e5)
% ];
% nCases = numel(cases);

% Storage
res_hist_num    = cell(nCases,1); 
res_hist_asy    = cell(nCases,1); 
opt_params_num  = zeros(nCases,2);
opt_params_asy  = zeros(nCases,2);
max_rho_num     = zeros(nCases,1);
max_rho_asy     = zeros(nCases,1);

global PQ_history
PQ_history = cell(nCases,1);

base_options = optimset('Display','off','TolX',1e-12,'TolFun',1e-12);

% Random initial guess for SWR
u_init = rand(Nx,Nt);

%% ============================================================
% MAIN LOOP
% ============================================================
for s = 1:nCases
    gamma = cases(s).gamma;
    nu    = cases(s).nu;
    label_s = sprintf('\\gamma=%.3g, \\nu=%.3g', gamma, nu);
    fprintf('\nCase %d: %s\n', s, label_s);
    
    % 1. Reference Solution
    u_ref = run_fdtd(u0, v0, Lx, T, c, dh, dt, gamma, nu);
    
    % 2. Asymptotic Formula
    [p_asy, q_asy, regime_name] = get_formula_pq(gamma, nu, c, omega_min, omega_max);
    opt_params_asy(s,:) = [p_asy, q_asy];
    fprintf('  Formula [%s]: p=%.4f, q=%.4f\n', regime_name, p_asy, q_asy);
    
    % Calculate Max Rho for Asymptotic
    max_rho_asy(s) = obj_Linf(N,T,dt,J,c,gamma,nu,a,M,p_asy,q_asy);
    
    % 3. ROBUST NUMERICAL OPTIMIZATION (Multi-Start)
    % ---------------------------------------------------------
    % We try 4 diverse starting points to avoid local minima
    guesses = [
        p_asy,       q_asy;      % 1. Asymptotic Formula (Good for Mixed/Diff)
        1.0/c,       0.0;        % 2. Pure Wave limit (Good for low gamma)
        1e-5,        q_asy;      % 3. Heat Equation limit (p->0, Good for high nu)
        p_asy*1.5,   q_asy*1.5   % 4. Perturbed guess
    ];
    
    best_obj = Inf;
    best_idx = 1;
    
    % Phase A: Scan all guesses (No trajectory logging)
    for k = 1:size(guesses, 1)
        x_start = guesses(k, :);
        if any(x_start < 0), continue; end
        
        % Robust Objective (Penalize negative values)
        robust_obj = @(x) check_positivity(x, N,T,dt,J,c,gamma,nu,a,M);
        
        [x_k, fval_k] = fminsearch(robust_obj, x_start, base_options);
        
        if fval_k < best_obj
            best_obj = fval_k;
            best_idx = k;
        end
    end
    
    % Phase B: Re-run the WINNER to capture trajectory
    x_winner = guesses(best_idx, :);
    outfun = @(x,optimvalues,state) store_trajectory(x,optimvalues,state,s);
    optim_options = optimset(base_options,'OutputFcn',outfun);
    robust_obj = @(x) check_positivity(x, N,T,dt,J,c,gamma,nu,a,M);
    
    [x_opt, fval] = fminsearch(robust_obj, x_winner, optim_options);
    
    p_num = x_opt(1);
    q_num = x_opt(2);
    opt_params_num(s,:) = x_opt;
    max_rho_num(s) = fval;
    
    fprintf('  Numeric [MultiStart]: p=%.4f, q=%.4f, MaxRho=%.4f (Start: %d)\n', ...
            p_num, q_num, fval, best_idx);
    
    % 4. SWR Runs
    if s < 5, k_iter = 500; else, k_iter = 500; end
    
    % --- A. Numerical Best Run ---
    [~, ~, hist_test_num] = run_swr(u0, v0, N, a, M, T, c, dh, dt, gamma, nu, p_num, q_num, 1, u_init, u_ref);
    err0_num = max(abs(u_init(:) - u_ref(:))) / max(abs(u_ref(:)));
    F_num = hist_test_num(1) / err0_num;
    [~, ~, hist_num] = run_swr(u0, v0, N, a, M, T, c, dh, dt, gamma, nu, p_num, q_num, k_iter, u_init/F_num, u_ref);
    res_hist_num{s} = [err0_num; hist_num(:)];
    
    % --- B. Asymptotic Formula Run ---
    [~, ~, hist_test_asy] = run_swr(u0, v0, N, a, M, T, c, dh, dt, gamma, nu, p_asy, q_asy, 1, u_init, u_ref);
    err0_asy = max(abs(u_init(:) - u_ref(:))) / max(abs(u_ref(:)));
    F_asy = hist_test_asy(1) / err0_asy;
    [~, ~, hist_asy] = run_swr(u0, v0, N, a, M, T, c, dh, dt, gamma, nu, p_asy, q_asy, k_iter, u_init/F_asy, u_ref);
    res_hist_asy{s} = [err0_asy; hist_asy(:)];
end

%% ============================================================
% PLOTTING ROUTINES (CONVERGENCE HISTORY)
% ============================================================
plot_convergence_group(cases, res_hist_num, res_hist_asy, 'gamma', 'Gamma Convergence', POS_AX_STD, FS_AXIS, FS_LABEL, LW_BOLD, LW_AXIS);
plot_convergence_group(cases, res_hist_num, res_hist_asy, 'nu',    'Nu Convergence',    POS_AX_STD, FS_AXIS, FS_LABEL, LW_BOLD, LW_AXIS);

%% ============================================================
% RHO ANALYSIS PLOTTING
% ============================================================
plot_rho_group_analysis(cases, opt_params_num, opt_params_asy, omega_min, omega_max, c, 'gamma', 'Rho(w): Telegrapher Cases');
plot_rho_group_analysis(cases, opt_params_num, opt_params_asy, omega_min, omega_max, c, 'nu',    'Rho(w): Viscoelastic Cases');

%% ============================================================
% MAX RHO COMPARISON BAR CHART
% ============================================================
plot_max_rho_bars(cases, max_rho_num, max_rho_asy, 'gamma', 'Max Rho: Telegrapher', FS_AXIS, FS_LABEL, FS_TITLE);
plot_max_rho_bars(cases, max_rho_num, max_rho_asy, 'nu',    'Max Rho: Viscoelastic', FS_AXIS, FS_LABEL, FS_TITLE);

%% ============================================================
% TRAJECTORY FIGURES
% ============================================================
for s = 1:nCases
    pq = PQ_history{s};
    if isempty(pq), continue; end
    gamma = cases(s).gamma;
    nu    = cases(s).nu;
    
    p_f = opt_params_asy(s,1); q_f = opt_params_asy(s,2);
    p_n = opt_params_num(s,1); q_n = opt_params_num(s,2);
    
    objfun_pq = @(p,q) obj_Linf(N,T,dt,J,c,gamma,nu,a,M,p,q);
    
    % Define plot bounds
    p_all = [pq(:,1); p_f]; q_all = [pq(:,2); q_f];
    pmin = min(p_all); pmax = max(p_all);
    qmin = min(q_all); qmax = max(q_all);
    ppad = 0.3*max(1e-3, pmax-pmin); qpad = 0.3*max(1e-3, qmax-qmin);
    
    pr = linspace(max(0, pmin-ppad), pmax+ppad, 50);
    qr = linspace(max(0, qmin-qpad), qmax+qpad, 50);
    [P,Q] = meshgrid(pr,qr);
    
    Z = zeros(size(P));
    for i=1:numel(P), Z(i) = objfun_pq(P(i), Q(i)); end
    
    figure('Name', sprintf('Traj Case %d', s), 'Color', 'w'); clf; hold on;
    contourf(P, Q, log10(Z), 30, 'LineStyle','none'); 
    cb = colorbar; cb.FontSize = FS_AXIS; cb.Label.String = 'log_{10}(Error)';
    
    plot(pq(:,1), pq(:,2), '-o', 'LineWidth', 2.0, 'Color','k', 'MarkerSize', 4);
    plot(pq(1,1), pq(1,2), 'ws', 'MarkerSize', 10, 'LineWidth', 2, 'MarkerFaceColor', 'k');
    plot(p_n, q_n, 'p', 'MarkerSize', 18, 'LineWidth', 2, 'MarkerFaceColor', 'r', 'MarkerEdgeColor','k');
    plot(p_f, q_f, 'h', 'MarkerSize', 18, 'LineWidth', 2, 'MarkerFaceColor', 'g', 'MarkerEdgeColor','k');
    
    title_str = sprintf('\\gamma=%.3g, \\nu=%.3g', gamma, nu);
    title(title_str, 'FontSize', FS_TITLE, 'FontWeight', 'bold');
    xlabel('p', 'FontSize', FS_LABEL); ylabel('q', 'FontSize', FS_LABEL);
    legend({'Landscape','Path','Start','Num. Opt','Asymp.'}, 'Location','best', 'FontSize', 14);
    grid on; set(gca,'FontSize',FS_AXIS, 'LineWidth', LW_AXIS, 'FontWeight', 'bold');
end

%% ============================================================
% LOCAL FUNCTIONS
% ============================================================

function val = check_positivity(x, N,T,dt,J,c,gamma,nu,a,M)
    % Penalizes negative parameters to enforce constraints in fminsearch
    if x(1) < 0 || x(2) < 0
        val = 1e9; % Soft constraint penalty
    else
        val = obj_Linf(N,T,dt,J,c,gamma,nu,a,M,x(1),x(2));
    end
end

function val = obj_Linf(N, T, dt, J, c, gamma, nu, a, M, theta1, theta2)
    b = a + M;
    omega_min = pi / T;
    omega_max = pi / dt;
    omegas = logspace(log10(omega_min), log10(omega_max), J);
    s_vals = 1i * omegas;
    r_vals = zeros(size(s_vals));
    for idx = 1:length(s_vals)
        r_vals(idx) = rho(N, s_vals(idx), theta1, theta2, c, gamma, nu, a, b);
    end
    val = max(abs(r_vals));
end

function r = rho(N, s, theta1, theta2, c, gamma, nu, a, b)
    if(N == 2)
        ikappa = s/c .* sqrt((1+gamma./s) ./ (1+nu.*s/c^2));
        num = ikappa - (theta1*s + theta2);
        den = ikappa + (theta1*s + theta2);
        r = num ./ den; % NO OVERLAP TERM
    else
        r = 1.0; 
    end
end

% function [p, q, regime] = get_formula_pq(gamma, nu, c, w_min, w_max, delta)
%     if nargin < 6, delta = 0; end
%     w_mid = sqrt(w_min * w_max);
%     z     = (w_min / w_max)^(0.25);
%     Y     = sqrt(z / (1 + z + z^2));
%     if nu > 1e-9
%         is_diffusive = (nu * w_min / c^2) > 1.0;
%         if ~is_diffusive
%             regime = 'Visco Prop';
%             is_big_overlap = (nu * w_min^2 * delta / (2*c^3)) > 1.0;
%             if is_big_overlap, q = (nu * w_min^2)/(2*c^3); p = 1/c - (3*nu^2*w_min^2)/(8*c^5); regime=[regime ' (Big \delta)'];
%             else, q = (nu*w_min*w_max)/(2*c^3); p = 1/c - (nu^2/(4*c^5))*(w_max^2 + 0.5*w_min^2); end
%         else
%             regime = 'Visco Diff';
%             is_big_overlap = (delta * sqrt(w_min/(2*nu))) > 1.0;
%             if is_big_overlap, q=sqrt(w_min/(2*nu)); p=1/sqrt(2*nu*w_min); regime=[regime ' (Big \delta)'];
%             else, p = Y/sqrt(2*nu*w_mid); q = p*w_mid; end
%         end
%     elseif gamma > 1e-9
%         is_diffusive = (gamma / w_max) > 1.0;
%         if ~is_diffusive
%             regime = 'Tele Prop';
%             is_big_overlap = (gamma*delta/(2*c)) > 1.0;
%             if is_big_overlap, q=gamma/(2*c); p=1/c+(gamma^2/(8*c*w_min^2)); regime=[regime ' (Big \delta)'];
%             else, q=gamma/(2*c); p=1/c+(gamma^2/(16*c))*(1/w_min^2 + 1/w_max^2); end
%         else
%             regime = 'Tele Diff';
%             is_big_overlap = (delta/c*sqrt(gamma*w_min/2)) > 1.0;
%             if is_big_overlap, q=sqrt(gamma*w_min)/(c*sqrt(2)); p=sqrt(gamma)/(c*sqrt(2*w_min)); regime=[regime ' (Big \delta)'];
%             else, factor=sqrt(gamma)/c; p=(factor/sqrt(2*w_mid))*Y; q=p*w_mid; end
%         end
%     else
%         regime = 'Pure Wave'; p = 1/c; q = 0;
%     end
% end

function [p, q, regime] = get_formula_pq(gamma, nu, c, w_min, w_max)
    % This function now implements the Geometric Mean (Phase Opposite)
    % condition: Lambda(w1)*Lambda(w2) = ikappa(w1)*ikappa(w2)
    
    % 1. Define Physics (Wavenumber handle)
    % ikappa = (i*w/c) * sqrt( (1 + gamma/(i*w)) / (1 + i*w*nu/c^2) )
    get_Z = @(w) (1i*w/c) .* sqrt( (1 + gamma./(1i*w)) ./ (1 + (1i*w*nu)/c^2) );
    
    % 2. Evaluate at Endpoints
    Z1 = get_Z(w_min);
    Z2 = get_Z(w_max);
    
    % 3. Target Product P = A + iB
    P = Z1 * Z2;
    A = real(P);
    B = imag(P);
    
    % 4. Solve Biquadratic Equation for p
    % Formula derived from: p^4(w1*w2) + A*p^2 - (B/(w1+w2))^2 = 0
    Wsum  = w_min + w_max;
    Wprod = w_min * w_max;
    K_term = (B / Wsum)^2;
    
    Discriminant = A^2 + 4 * Wprod * K_term;
    
    % p^2 solution
    p_sq = (-A + sqrt(Discriminant)) / (2 * Wprod);
    
    % Final Parameters
    p = sqrt(p_sq);
    q = abs( B / (p * Wsum) ); % abs() ensures physical stability (positivity)
    
    % Identify Regime for labeling
    if nu > 1e-9 && gamma > 1e-9, regime = 'Full Visco-Tele';
    elseif nu > 1e-9, regime = 'Viscoelastic';
    elseif gamma > 1e-9, regime = 'Telegrapher';
    else, regime = 'Pure Wave'; end
end

function plot_convergence_group(cases, hist_num, hist_asy, type, fig_name, pos, fs_ax, fs_lbl, lw_bold, lw_ax)
    figure('Name', fig_name, 'Color', 'w'); clf; set(gcf, 'Position', [100 100 900 600]); 
    if strcmp(type, 'gamma'), idx_list = find([cases.nu] == 0); else, idx_list = find([cases.gamma] == 0); end
    h_plots = []; labels = {}; colors = lines(numel(idx_list));
    for j = 1:numel(idx_list)
        s = idx_list(j); col = colors(j,:);
        res = hist_num{s}; hp = semilogy(0:numel(res)-1, res, '-', 'LineWidth', lw_bold, 'Color', col); hold on;
        h_plots(end+1) = hp;
        res_a = hist_asy{s}; semilogy(0:numel(res_a)-1, res_a, '--', 'LineWidth', lw_bold, 'Color', col, 'HandleVisibility','off');
        if strcmp(type, 'gamma'), lbl=sprintf('\\gamma=%.1e', cases(s).gamma); else, lbl=sprintf('\\nu=%.1e', cases(s).nu); end
        labels{end+1} = lbl;
    end
    grid on; xlabel('Iteration k', 'FontSize', fs_lbl); ylabel('Error', 'FontSize', fs_lbl);
    set(gca, 'FontSize', fs_ax, 'LineWidth', lw_ax, 'FontWeight', 'bold', 'Position', pos);
    h_solid = plot(nan,nan, 'k-', 'LineWidth', 2); h_dash = plot(nan,nan, 'k--', 'LineWidth', 2);
    legend([h_plots, h_solid, h_dash], [labels, {'Num. Opt.', 'Formula'}], 'Location', 'SouthWest', 'FontSize', 16);
end

function plot_rho_group_analysis(cases, params_num, params_asy, w_min, w_max, c, type, fig_name)
    figure('Name', fig_name, 'Color', 'w'); set(gcf, 'Position', [150 150 1000 700]);
    t = tiledlayout(2,2, 'TileSpacing', 'compact', 'Padding', 'compact');
    if strcmp(type, 'gamma'), idx_list = find([cases.nu] == 0); else, idx_list = find([cases.gamma] == 0); end
    w_eval = logspace(log10(w_min), log10(w_max), 500);
    for j = 1:min(4, numel(idx_list))
        s = idx_list(j); gamma = cases(s).gamma; nu = cases(s).nu;
        if nu > 1e-9, denom = sqrt(1 + 1i.*w_eval.*nu./c^2); ik_x = (1i.*w_eval./c)./denom;
        else, term = sqrt(1 + gamma./(1i.*w_eval)); ik_x = (1i.*w_eval./c).*term; end
        p_n=params_num(s,1); q_n=params_num(s,2); L_n=q_n+1i.*w_eval.*p_n; rho_n=abs((ik_x-L_n)./(ik_x+L_n));
        p_a=params_asy(s,1); q_a=params_asy(s,2); L_a=q_a+1i.*w_eval.*p_a; rho_a=abs((ik_x-L_a)./(ik_x+L_a));
        nexttile; hold on; plot(w_eval, rho_n, 'r-', 'LineWidth', 2.5); plot(w_eval, rho_a, 'b--', 'LineWidth', 2.5);
        set(gca, 'XScale', 'log', 'YScale', 'log'); xlim([w_min, w_max]); ylim([1e-5, 1.1]); grid on;
        if strcmp(type, 'gamma'), title(sprintf('\\gamma = %.1e', gamma), 'FontSize', 14); else, title(sprintf('\\nu = %.1e', nu), 'FontSize', 14); end
    end
    lg = legend({'Numeric', 'Formula'}, 'Orientation', 'horizontal'); lg.Layout.Tile = 'north'; title(t, fig_name, 'FontSize', 16, 'FontWeight', 'bold');
end

function plot_max_rho_bars(cases, rho_num, rho_asy, type, fig_name, fs_ax, fs_lbl, fs_title)
    figure('Name', fig_name, 'Color', 'w'); set(gcf, 'Position', [200 200 800 500]);
    if strcmp(type, 'gamma'), idx_list = find([cases.nu] == 0); x_cats = categorical(cellstr(num2str([cases(idx_list).gamma]','%.0e'))); xlab='Damping \gamma';
    else, idx_list = find([cases.gamma] == 0); x_cats = categorical(cellstr(num2str([cases(idx_list).nu]','%.0e'))); xlab='Viscosity \nu'; end
    data = [rho_num(idx_list), rho_asy(idx_list)];
    b = bar(x_cats, data, 'grouped'); b(1).FaceColor = 'r'; b(2).FaceColor = 'b';
    ylabel('Max |\rho|', 'FontSize', fs_lbl); xlabel(xlab, 'FontSize', fs_lbl); title(fig_name, 'FontSize', fs_title);
    legend({'Numeric Opt', 'Asymptotic Formula'}, 'Location', 'best', 'FontSize', 16); grid on; set(gca, 'FontSize', fs_ax, 'FontWeight', 'bold');
end

function stop = store_trajectory(x,optimvalues,state,idx)
    stop = false; global PQ_history;
    if strcmp(state,'init'), PQ_history{idx} = x(:).'; elseif strcmp(state,'iter'), PQ_history{idx}(end+1,:) = x(:).'; end
end