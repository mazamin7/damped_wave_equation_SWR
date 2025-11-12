clear; close all; clc;

addpath("utils\")

%% Common parameters
N = 2;          % Number of subdomains
c = 1.0;
M = 0.1;
a = 0.3;
b = a + M;
T = 5;
ky = 0;         % 1D case

% Discretization parameters
dh = 0.01;
dt = 0.01;
J = 50;         % Number of frequency steps

% Optimization parameters
tol = 1e-9;
tightness = 0.5;

%% Define parameter sweep cases
% Use log-spaced values for better coverage
gamma_values = logspace(-1, 1, 61);
nu_values = logspace(-2, 0, 61);
nu_fixed_values = logspace(-2, 0, 11);
gamma_fixed_values = logspace(-1, 1, 11);

num_gamma_cases = length(nu_fixed_values);
num_nu_cases = length(gamma_fixed_values);

%% Initialize storage arrays
optimal_pinf_gamma = cell(num_gamma_cases, 1);
contraction_pinf_gamma = cell(num_gamma_cases, 1);
optimal_pinf_nu = cell(num_nu_cases, 1);
contraction_pinf_nu = cell(num_nu_cases, 1);

%% Optimization settings
optim_options = optimset('Display', 'off', 'TolX', tol, 'TolFun', tol, ...
                        'MaxFunEvals', 2e4, 'MaxIter', 2e4);

%% Helper functions
find_nearest_index = @(v, x) find(abs(v - x) == min(abs(v - x)), 1, 'first');
x_fallback = [1/c, 0];

%% Process parameter sweeps
fprintf('=== Processing gamma-Varying Cases (L∞ only) ===\n');
for case_idx = 1:num_gamma_cases
    nu_fixed = nu_fixed_values(case_idx);
    fprintf('\n=== Case %d: Varying gamma (nu = %.2f) ===\n', case_idx, nu_fixed);
    
    optimal_pinf_gamma{case_idx} = zeros(length(gamma_values), 2);
    contraction_pinf_gamma{case_idx} = zeros(length(gamma_values), 1);
    x0 = [1/c, 0];
    
    for i = 1:length(gamma_values)
        gamma = gamma_values(i);
        fprintf('Computing L∞ for gamma = %.3f (%d/%d)\n', gamma, i, length(gamma_values));
        
        objfun = @(x) obj_Linf(N, T, dt, J, c, gamma, nu_fixed, a, M, x(1), x(2), ky);
        [x_opt, fval] = fminsearch(objfun, x0, optim_options);
        
        optimal_pinf_gamma{case_idx}(i,:) = x_opt;
        contraction_pinf_gamma{case_idx}(i) = fval;
        x0 = x_opt;
    end
end

fprintf('\n=== Processing nu-Varying Cases (L∞ only) ===\n');
for case_idx = 1:num_nu_cases
    gamma_fixed = gamma_fixed_values(case_idx);
    fprintf('\n=== Case %d: Varying nu (gamma = %.2f) ===\n', case_idx, gamma_fixed);
    
    optimal_pinf_nu{case_idx} = zeros(length(nu_values), 2);
    contraction_pinf_nu{case_idx} = zeros(length(nu_values), 1);
    x0 = [1/c, 0];
    
    for i = 1:length(nu_values)
        nu = nu_values(i);
        fprintf('Computing L∞ for nu = %.3f (%d/%d)\n', nu, i, length(nu_values));
        
        objfun = @(x) obj_Linf(N, T, dt, J, c, gamma_fixed, nu, a, M, x(1), x(2), ky);
        [x_opt, fval] = fminsearch(objfun, x0, optim_options);
        
        optimal_pinf_nu{case_idx}(i,:) = x_opt;
        contraction_pinf_nu{case_idx}(i) = fval;
        x0 = x_opt;
    end
end

%% Cross-initialization refinement
fprintf('\n=== Cross-initialization Refinement ===\n');

% Refine gamma-varying cases using nu-sweep results
for case_idx = 1:num_gamma_cases
    nu_fixed = nu_fixed_values(case_idx);
    nu_index = find_nearest_index(nu_values, nu_fixed);
    
    for i = 1:length(gamma_values)
        gamma = gamma_values(i);
        nu_case_idx = find_nearest_index(gamma_fixed_values, gamma);
        x0 = optimal_pinf_nu{nu_case_idx}(nu_index, :);
        
        if any(~isfinite(x0)), x0 = x_fallback; end
        
        objfun = @(x) obj_Linf(N, T, dt, J, c, gamma, nu_fixed, a, M, x(1), x(2), ky);
        [x_opt, fval] = fminsearch(objfun, x0, optim_options);
        
        optimal_pinf_gamma{case_idx}(i,:) = x_opt;
        contraction_pinf_gamma{case_idx}(i) = fval;
    end
end

% Refine nu-varying cases using gamma-sweep results
for case_idx = 1:num_nu_cases
    gamma_fixed = gamma_fixed_values(case_idx);
    gamma_index = find_nearest_index(gamma_values, gamma_fixed);
    
    for i = 1:length(nu_values)
        nu = nu_values(i);
        gamma_case_idx = find_nearest_index(nu_fixed_values, nu);
        x0 = optimal_pinf_gamma{gamma_case_idx}(gamma_index, :);
        
        if any(~isfinite(x0)), x0 = x_fallback; end
        
        objfun = @(x) obj_Linf(N, T, dt, J, c, gamma_fixed, nu, a, M, x(1), x(2), ky);
        [x_opt, fval] = fminsearch(objfun, x0, optim_options);
        
        optimal_pinf_nu{case_idx}(i,:) = x_opt;
        contraction_pinf_nu{case_idx}(i) = fval;
    end
end

%% Reconcile overlapping points
fprintf('\n=== Reconciling Overlapping Points ===\n');

tol_match = 1e-12;
[gamma_match, gamma_indices] = ismembertol(gamma_fixed_values, gamma_values, tol_match);
[nu_match, nu_indices] = ismembertol(nu_fixed_values, nu_values, tol_match);

for gamma_fix_idx = 1:length(gamma_fixed_values)
    if ~gamma_match(gamma_fix_idx), continue; end
    
    for nu_fix_idx = 1:length(nu_fixed_values)
        if ~nu_match(nu_fix_idx), continue; end
        
        gamma_sweep_idx = gamma_indices(gamma_fix_idx);
        nu_sweep_idx = nu_indices(nu_fix_idx);
        
        x_gamma = optimal_pinf_gamma{nu_fix_idx}(gamma_sweep_idx, :);
        rho_gamma = contraction_pinf_gamma{nu_fix_idx}(gamma_sweep_idx);
        x_nu = optimal_pinf_nu{gamma_fix_idx}(nu_sweep_idx, :);
        rho_nu = contraction_pinf_nu{gamma_fix_idx}(nu_sweep_idx);
        
        gamma_valid = all(isfinite(x_gamma)) && isfinite(rho_gamma);
        nu_valid = all(isfinite(x_nu)) && isfinite(rho_nu);
        
        if ~gamma_valid && ~nu_valid
            continue;
        elseif gamma_valid && ~nu_valid
            optimal_pinf_nu{gamma_fix_idx}(nu_sweep_idx, :) = x_gamma;
            contraction_pinf_nu{gamma_fix_idx}(nu_sweep_idx) = rho_gamma;
        elseif ~gamma_valid && nu_valid
            optimal_pinf_gamma{nu_fix_idx}(gamma_sweep_idx, :) = x_nu;
            contraction_pinf_gamma{nu_fix_idx}(gamma_sweep_idx) = rho_nu;
        else
            if rho_nu < rho_gamma
                optimal_pinf_gamma{nu_fix_idx}(gamma_sweep_idx, :) = x_nu;
                contraction_pinf_gamma{nu_fix_idx}(gamma_sweep_idx) = rho_nu;
            elseif rho_gamma < rho_nu
                optimal_pinf_nu{gamma_fix_idx}(nu_sweep_idx, :) = x_gamma;
                contraction_pinf_nu{gamma_fix_idx}(nu_sweep_idx) = rho_gamma;
            else
                if norm(x_nu) < norm(x_gamma)
                    optimal_pinf_gamma{nu_fix_idx}(gamma_sweep_idx, :) = x_nu;
                    contraction_pinf_gamma{nu_fix_idx}(gamma_sweep_idx) = rho_nu;
                else
                    optimal_pinf_nu{gamma_fix_idx}(nu_sweep_idx, :) = x_gamma;
                    contraction_pinf_nu{gamma_fix_idx}(nu_sweep_idx) = rho_gamma;
                end
            end
        end
    end
end

%% Create fine grid for global contraction factor
fprintf('\n=== Creating Fine Grid for Global Contraction Factor ===\n');

% Collect all data points for interpolation
all_gamma = []; all_nu = []; all_contraction = []; all_p = []; all_q = [];

% Add gamma-sweep data
for case_idx = 1:num_gamma_cases
    nu_fixed = nu_fixed_values(case_idx);
    for i = 1:length(gamma_values)
        if isfinite(contraction_pinf_gamma{case_idx}(i))
            all_gamma = [all_gamma; gamma_values(i)];
            all_nu = [all_nu; nu_fixed];
            all_contraction = [all_contraction; contraction_pinf_gamma{case_idx}(i)];
            all_p = [all_p; optimal_pinf_gamma{case_idx}(i, 1)];
            all_q = [all_q; optimal_pinf_gamma{case_idx}(i, 2)];
        end
    end
end

% Add nu-sweep data
for case_idx = 1:num_nu_cases
    gamma_fixed = gamma_fixed_values(case_idx);
    for i = 1:length(nu_values)
        if isfinite(contraction_pinf_nu{case_idx}(i))
            all_gamma = [all_gamma; gamma_fixed];
            all_nu = [all_nu; nu_values(i)];
            all_contraction = [all_contraction; contraction_pinf_nu{case_idx}(i)];
            all_p = [all_p; optimal_pinf_nu{case_idx}(i, 1)];
            all_q = [all_q; optimal_pinf_nu{case_idx}(i, 2)];
        end
    end
end

% Create scattered interpolants
F_contraction = scatteredInterpolant(all_gamma, all_nu, all_contraction, 'natural', 'none');
F_p = scatteredInterpolant(all_gamma, all_nu, all_p, 'natural', 'none');
F_q = scatteredInterpolant(all_gamma, all_nu, all_q, 'natural', 'none');

% Create meshgrid for fine interpolation
[GAMMA, NU] = meshgrid(gamma_values, nu_values);

% Interpolate onto fine grid
fprintf('Performing interpolation...\n');
contraction_fine = F_contraction(GAMMA, NU);
p_fine = F_p(GAMMA, NU);
q_fine = F_q(GAMMA, NU);

% Find minimum contraction factor
[min_contraction, min_idx] = min(contraction_fine(:));
[min_row, min_col] = ind2sub(size(contraction_fine), min_idx);
min_gamma = gamma_values(min_col);
min_nu = nu_values(min_row);
min_p = p_fine(min_row, min_col);
min_q = q_fine(min_row, min_col);

fprintf('Minimum contraction factor: %.6f\n', min_contraction);
fprintf('At gamma=%.3f, nu=%.3f\n', min_gamma, min_nu);
fprintf('Optimal parameters: p=%.6f, q=%.6f\n', min_p, min_q);

%% Plot results
plot_basic_figures(gamma_values, nu_fixed_values, optimal_pinf_gamma, contraction_pinf_gamma, ...
                   nu_values, gamma_fixed_values, optimal_pinf_nu, contraction_pinf_nu);

plot_joint_view_figures(gamma_values, nu_fixed_values, optimal_pinf_gamma, contraction_pinf_gamma, ...
                        nu_values, gamma_fixed_values, optimal_pinf_nu, contraction_pinf_nu, ...
                        gamma_values, nu_values, contraction_fine, tightness);

%% Helper function for basic plotting (legend markers label the SWEPT parameter)
function plot_basic_figures(gamma_vals, nu_fixed_vals, opt_gamma, contr_gamma, ...
                           nu_vals, gamma_fixed_vals, opt_nu, contr_nu)
    
    num_gamma_cases = length(nu_fixed_vals);
    num_nu_cases    = length(gamma_fixed_vals);
    colors = lines(max(num_gamma_cases, num_nu_cases));
    
    % -------------------------
    % Figure 1: Varying gamma (fixed ν)
    % -------------------------
    figure('Position', [100, 100, 1200, 500], 'Name', 'Gamma-Varying Cases');
    
    % Subplot 1: (p,q) trajectories
    subplot(1, 2, 1); hold on;
    hLines = gobjects(num_gamma_cases,1);
    names  = cell(num_gamma_cases,1);
    for case_idx = 1:num_gamma_cases
        p_vals = opt_gamma{case_idx}(:, 1);
        q_vals = opt_gamma{case_idx}(:, 2);
        
        names{case_idx} = sprintf('\\nu=%.2f', nu_fixed_vals(case_idx));
        hLines(case_idx) = plot(p_vals, q_vals, '-', 'LineWidth', 2, ...
             'DisplayName', names{case_idx}, 'Color', colors(case_idx, :));
        
        % Endpoints (hidden from legend)
        plot(p_vals(1),  q_vals(1),  'o', 'MarkerSize', 8, ...
             'MarkerFaceColor', colors(case_idx, :), 'MarkerEdgeColor', 'k', 'HandleVisibility', 'off');
        plot(p_vals(end), q_vals(end), 's', 'MarkerSize', 8, ...
             'MarkerFaceColor', colors(case_idx, :), 'MarkerEdgeColor', 'k', 'HandleVisibility', 'off');
    end
    % Dummy handles ONCE in legend, labeled by the SWEPT parameter (gamma)
    hMinSweep = plot(nan, nan, 'ok', 'MarkerSize', 8, 'MarkerFaceColor', 'w'); % circle -> gamma min
    hMaxSweep = plot(nan, nan, 'sk', 'MarkerSize', 8, 'MarkerFaceColor', 'w'); % square -> gamma max
    
    xlabel('p', 'FontSize', 16);
    ylabel('q', 'FontSize', 16);
    grid on; xlim([0, 1.5]); ylim([-4, 8]); set(gca, 'FontSize', 16);
    labels = [names; ...
        {sprintf('\\gamma=%.2f', min(gamma_vals))}; ...
        {sprintf('\\gamma=%.2f', max(gamma_vals))}];
    legend([hLines; hMinSweep; hMaxSweep], labels, 'Location', 'NorthEast', 'FontSize', 12);
    
    % Subplot 2: Contraction vs gamma
    subplot(1, 2, 2); hold on;
    for case_idx = 1:num_gamma_cases
        semilogy(gamma_vals, contr_gamma{case_idx}, '-', 'LineWidth', 2, ...
                 'DisplayName', sprintf('\\nu=%.2f', nu_fixed_vals(case_idx)), ...
                 'Color', colors(case_idx, :));
    end
    xlabel('\gamma', 'FontSize', 16);
    ylabel('Global Contraction Factor', 'FontSize', 16);
    grid on; ylim([1e-3, 1]); set(gca, 'FontSize', 16);
    legend('Location', 'NorthEast', 'FontSize', 12);
    
    % -------------------------
    % Figure 2: Varying ν (fixed γ)
    % -------------------------
    figure('Position', [100, 100, 1200, 500], 'Name', 'Nu-Varying Cases');
    
    % Subplot 1: (p,q) trajectories
    subplot(1, 2, 1); hold on;
    hLines = gobjects(num_nu_cases,1);
    names  = cell(num_nu_cases,1);
    for case_idx = 1:num_nu_cases
        p_vals = opt_nu{case_idx}(:, 1);
        q_vals = opt_nu{case_idx}(:, 2);
        
        names{case_idx} = sprintf('\\gamma=%.2f', gamma_fixed_vals(case_idx));
        hLines(case_idx) = plot(p_vals, q_vals, '-', 'LineWidth', 2, ...
             'DisplayName', names{case_idx}, 'Color', colors(case_idx, :));
        
        % Endpoints (hidden from legend)
        plot(p_vals(1),  q_vals(1),  'o', 'MarkerSize', 8, ...
             'MarkerFaceColor', colors(case_idx, :), 'MarkerEdgeColor', 'k', 'HandleVisibility', 'off');
        plot(p_vals(end), q_vals(end), 's', 'MarkerSize', 8, ...
             'MarkerFaceColor', colors(case_idx, :), 'MarkerEdgeColor', 'k', 'HandleVisibility', 'off');
    end
    % Dummy handles ONCE in legend, labeled by the SWEPT parameter (nu)
    hMinSweep = plot(nan, nan, 'ok', 'MarkerSize', 8, 'MarkerFaceColor', 'w'); % circle -> nu min
    hMaxSweep = plot(nan, nan, 'sk', 'MarkerSize', 8, 'MarkerFaceColor', 'w'); % square -> nu max
    
    xlabel('p', 'FontSize', 16);
    ylabel('q', 'FontSize', 16);
    grid on; xlim([0, 1.5]); ylim([-4, 8]); set(gca, 'FontSize', 16);
    labels = [names; ...
        {sprintf('\\nu=%.2f', min(nu_vals))}; ...
        {sprintf('\\nu=%.2f', max(nu_vals))}];
    legend([hLines; hMinSweep; hMaxSweep], labels, 'Location', 'NorthEast', 'FontSize', 12);
    
    % Subplot 2: Contraction vs nu
    subplot(1, 2, 2); hold on;
    for case_idx = 1:num_nu_cases
        semilogy(nu_vals, contr_nu{case_idx}, '-', 'LineWidth', 2, ...
                 'DisplayName', sprintf('\\gamma=%.2f', gamma_fixed_vals(case_idx)), ...
                 'Color', colors(case_idx, :));
    end
    xlabel('\nu', 'FontSize', 16);
    ylabel('Global Contraction Factor', 'FontSize', 16);
    grid on; ylim([1e-3, 1]); set(gca, 'FontSize', 16);
    legend('Location', 'NorthEast', 'FontSize', 12);
end

%% Helper function for joint view plotting with boundaries (legend markers label the SWEPT parameter)
function plot_joint_view_figures(gamma_vals, nu_fixed_vals, opt_gamma, contr_gamma, ...
                                nu_vals, gamma_fixed_vals, opt_nu, contr_nu, ...
                                gamma_fine, nu_fine, contraction_fine, tightness)
    % -------------------------
    % Varying gamma (fixed ν)
    % -------------------------
    V  = nu_fixed_vals(:)';        
    Nv = length(V);
    
    p_gam = zeros(Nv, numel(gamma_vals));
    q_gam = zeros(Nv, numel(gamma_vals));
    for case_idx = 1:Nv
        p_gam(case_idx, :) = opt_gamma{case_idx}(:, 1)';
        q_gam(case_idx, :) = opt_gamma{case_idx}(:, 2)';
    end
    
    figure(); hold on;
    % use a blue-ish gradient instead of grey
    cmap_gamma = parula(Nv+5);
    hLines = gobjects(Nv,1);
    names  = cell(Nv,1);
    for case_idx = 1:Nv
        this_color = cmap_gamma(case_idx, :);
        names{case_idx} = sprintf('\\nu=%.2f', V(case_idx));
        hLines(case_idx) = plot(p_gam(case_idx, :), q_gam(case_idx, :), '-', ...
            'LineWidth', 2, 'Color', this_color, 'DisplayName', names{case_idx});
        % Endpoints (hidden)
        plot(p_gam(case_idx, 1),  q_gam(case_idx, 1),  'o', 'MarkerSize', 6, ...
            'MarkerFaceColor', this_color, 'MarkerEdgeColor', 'k', 'HandleVisibility', 'off');
        plot(p_gam(case_idx, end), q_gam(case_idx, end), 's', 'MarkerSize', 6, ...
            'MarkerFaceColor', this_color, 'MarkerEdgeColor', 'k', 'HandleVisibility', 'off');
    end
    
    % Boundary (varying gamma)
    Pg_all = p_gam(:); Qg_all = q_gam(:);
    mask = isfinite(Pg_all) & isfinite(Qg_all);
    Pg_all = Pg_all(mask); Qg_all = Qg_all(mask);
    if ~isempty(Pg_all)
        try
            ash_g = alphaShape(Pg_all, Qg_all); 
            ash_g.Alpha = 0.85 * criticalAlpha(ash_g);
            plot(ash_g, 'FaceColor', 'none', 'EdgeColor', 'k', 'LineWidth', 1.4, 'HandleVisibility', 'off');
        catch
            bi_g = boundary(Pg_all, Qg_all, tightness);
            plot(Pg_all(bi_g), Qg_all(bi_g), 'r--', 'LineWidth', 1.4, 'HandleVisibility', 'off');
        end
    end
    
    xlabel('p'); ylabel('q');
    % title('Optimal (p,q) — Varying \gamma (Fixed \nu)');
    grid on; box on; xlim([0, 1.5]); ylim([-4, 8]);
    % Legend entries for the SWEPT parameter (gamma)
    hMinSweep = plot(nan, nan, 'ok', 'MarkerSize', 8, 'MarkerFaceColor', 'w'); % circle -> gamma min
    hMaxSweep = plot(nan, nan, 'sk', 'MarkerSize', 8, 'MarkerFaceColor', 'w'); % square -> gamma max
    labels = [names; ...
        {sprintf('\\gamma=%.2f', min(gamma_vals))}; ...
        {sprintf('\\gamma=%.2f', max(gamma_vals))}];
    legend([hLines; hMinSweep; hMaxSweep], labels, 'Location', 'northeast', 'FontSize', 12);

    % Make tick labels bigger
    set(gca, 'FontSize', 18);
    
    % -------------------------
    % Varying ν (fixed γ)
    % -------------------------
    Gf  = gamma_fixed_vals(:)';    
    Ngf = length(Gf);
    p_nu = zeros(Ngf, numel(nu_vals));
    q_nu = zeros(Ngf, numel(nu_vals));
    for case_idx = 1:Ngf
        p_nu(case_idx, :) = opt_nu{case_idx}(:, 1)';
        q_nu(case_idx, :) = opt_nu{case_idx}(:, 2)';
    end
    
    figure(); hold on;
    % use a red/orange gradient instead of grey
    cmap_nu = parula(Ngf+5);
    hLines = gobjects(Ngf,1);
    names  = cell(Ngf,1);
    for case_idx = 1:Ngf
        this_color = cmap_nu(case_idx, :);
        names{case_idx} = sprintf('\\gamma=%.2f', Gf(case_idx));
        hLines(case_idx) = plot(p_nu(case_idx, :), q_nu(case_idx, :), '-', ...
            'LineWidth', 2, 'Color', this_color, 'DisplayName', names{case_idx});
        % Endpoints (hidden)
        plot(p_nu(case_idx, 1),  q_nu(case_idx, 1),  'o', 'MarkerSize', 6, ...
            'MarkerFaceColor', this_color, 'MarkerEdgeColor', 'k', 'HandleVisibility', 'off');
        plot(p_nu(case_idx, end), q_nu(case_idx, end), 's', 'MarkerSize', 6, ...
            'MarkerFaceColor', this_color, 'MarkerEdgeColor', 'k', 'HandleVisibility', 'off');
    end
    
    % Boundary (varying nu)
    Pn_all = p_nu(:); Qn_all = q_nu(:);
    mask = isfinite(Pn_all) & isfinite(Qn_all);
    Pn_all = Pn_all(mask); Qn_all = Qn_all(mask);
    if ~isempty(Pn_all)
        try
            ash_n = alphaShape(Pn_all, Qn_all); 
            ash_n.Alpha = 0.85 * criticalAlpha(ash_n);
            plot(ash_n, 'FaceColor', 'none', 'EdgeColor', 'k', 'LineWidth', 1.4, 'HandleVisibility', 'off');
        catch
            bi_n = boundary(Pn_all, Qn_all, tightness);
            plot(Pn_all(bi_n), Qn_all(bi_n), 'r--', 'LineWidth', 1.4, 'HandleVisibility', 'off');
        end
    end
    
    xlabel('p'); ylabel('q');
    % title('Optimal (p,q) — Varying \nu (Fixed \gamma)');
    grid on; box on; xlim([0, 1.5]); ylim([-4, 8]);
    % Legend entries for the SWEPT parameter (nu)
    hMinSweep = plot(nan, nan, 'ok', 'MarkerSize', 8, 'MarkerFaceColor', 'w'); % circle -> nu min
    hMaxSweep = plot(nan, nan, 'sk', 'MarkerSize', 8, 'MarkerFaceColor', 'w'); % square -> nu max
    labels = [names; ...
        {sprintf('\\nu=%.2f', min(nu_vals))}; ...
        {sprintf('\\nu=%.2f', max(nu_vals))}];
    legend([hLines; hMinSweep; hMaxSweep], labels, 'Location', 'northeast', 'FontSize', 12);

    % Make tick labels bigger
    set(gca, 'FontSize', 18);
    
    % -------------------------
    % Boundaries only (Figure 3)
    % -------------------------
    if exist('bi_g','var')
        idx_g = bi_g;
    else
        idx_g = boundary(Pg_all, Qg_all, tightness);
    end
    if exist('bi_n','var')
        idx_n = bi_n;
    else
        idx_n = boundary(Pn_all, Qn_all, tightness);
    end
    
    figure(); hold on;
    plot(Pg_all(idx_g), Qg_all(idx_g), 'k-', 'LineWidth', 2, ...
        'DisplayName', 'Boundary: varying \gamma (fixed \nu)');
    plot(Pn_all(idx_n), Qn_all(idx_n), 'r--', 'LineWidth', 2, ...
        'DisplayName', 'Boundary: varying \nu (fixed \gamma)');
    xlabel('p', 'FontSize', 16); ylabel('q', 'FontSize', 16);
    % title('Boundary Comparison in (p,q) Space', 'FontSize', 12);
    grid on; box on; xlim([0, 1.5]); ylim([-4, 8]);
    legend('Location', 'northeast', 'FontSize', 16);

    % Make tick labels bigger
    set(gca, 'FontSize', 18);
    
    % -------------------------
    % Contraction contour (log-log) (Figure 4)
    % -------------------------
    figure();
    contourf(gamma_fine, nu_fine, log10(contraction_fine), 20, 'LineColor', 'none');
    colormap(parula); colorbar;
    set(gca, 'XScale', 'log', 'YScale', 'log');
    xlabel('\gamma', 'FontSize', 16); ylabel('\nu', 'FontSize', 16);
    % title('Optimal Contraction Factor - Log-Log-Log Scale', 'FontSize', 12);
    grid on; box on; set(gca, 'FontSize', 12);
    set(gca, 'XMinorGrid', 'on', 'YMinorGrid', 'on');
    xlim([1e-1, 1e1]); ylim([1e-2, 1e0]);

    % Make tick labels bigger
    set(gca, 'FontSize', 18);
end
