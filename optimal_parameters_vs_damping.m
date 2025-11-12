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
gamma_values = logspace(-2, 1, 91);
nu_values = gamma_values;
nu_fixed_values = logspace(-2, 1, 31);
gamma_fixed_values = nu_fixed_values;

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

%% Helper function for basic plotting
function plot_basic_figures(gamma_vals, nu_fixed_vals, opt_gamma, contr_gamma, ...
                           nu_vals, gamma_fixed_vals, opt_nu, contr_nu)
    
    num_gamma_cases = length(nu_fixed_vals);
    num_nu_cases = length(gamma_fixed_vals);
    colors = lines(max(num_gamma_cases, num_nu_cases));
    
    % Figure 1: Gamma-varying cases
    figure('Position', [100, 100, 1200, 500], 'Name', 'Gamma-Varying Cases');
    
    % Subplot 1: Parameter space
    subplot(1, 2, 1);
    hold on;
    for case_idx = 1:num_gamma_cases
        p_vals = opt_gamma{case_idx}(:, 1);
        q_vals = opt_gamma{case_idx}(:, 2);
        
        plot(p_vals, q_vals, '-', 'LineWidth', 2, ...
             'DisplayName', sprintf('\\nu=%.2f', nu_fixed_vals(case_idx)), ...
             'Color', colors(case_idx, :));
        
        % Mark start and end points
        plot(p_vals(1), q_vals(1), 'o', 'MarkerSize', 8, ...
             'MarkerFaceColor', colors(case_idx, :), 'MarkerEdgeColor', 'k', 'HandleVisibility', 'off');
        plot(p_vals(end), q_vals(end), 's', 'MarkerSize', 8, ...
             'MarkerFaceColor', colors(case_idx, :), 'MarkerEdgeColor', 'k', 'HandleVisibility', 'off');
    end
    
    xlabel('p', 'FontSize', 16);
    ylabel('q', 'FontSize', 16);
    legend('Location', 'NorthEast', 'FontSize', 12);
    grid on;
    xlim([0, 1.5]);
    ylim([-4, 8]);
    set(gca, 'FontSize', 14);
    
    % Subplot 2: Contraction factors
    subplot(1, 2, 2);
    hold on;
    for case_idx = 1:num_gamma_cases
        semilogy(gamma_vals, contr_gamma{case_idx}, '-', 'LineWidth', 2, ...
                 'DisplayName', sprintf('\\nu=%.2f', nu_fixed_vals(case_idx)), ...
                 'Color', colors(case_idx, :));
    end
    xlabel('\gamma', 'FontSize', 16);
    ylabel('Global Contraction Factor', 'FontSize', 16);
    legend('Location', 'NorthEast', 'FontSize', 12);
    grid on;
    ylim([1e-3, 1]);
    set(gca, 'FontSize', 14);
    
    % Figure 2: Nu-varying cases
    figure('Position', [100, 100, 1200, 500], 'Name', 'Nu-Varying Cases');
    
    % Subplot 1: Parameter space
    subplot(1, 2, 1);
    hold on;
    for case_idx = 1:num_nu_cases
        p_vals = opt_nu{case_idx}(:, 1);
        q_vals = opt_nu{case_idx}(:, 2);
        
        plot(p_vals, q_vals, '-', 'LineWidth', 2, ...
             'DisplayName', sprintf('\\gamma=%.2f', gamma_fixed_vals(case_idx)), ...
             'Color', colors(case_idx, :));
        
        % Mark start and end points
        plot(p_vals(1), q_vals(1), 'o', 'MarkerSize', 8, ...
             'MarkerFaceColor', colors(case_idx, :), 'MarkerEdgeColor', 'k', 'HandleVisibility', 'off');
        plot(p_vals(end), q_vals(end), 's', 'MarkerSize', 8, ...
             'MarkerFaceColor', colors(case_idx, :), 'MarkerEdgeColor', 'k', 'HandleVisibility', 'off');
    end
    
    xlabel('p', 'FontSize', 16);
    ylabel('q', 'FontSize', 16);
    legend('Location', 'NorthEast', 'FontSize', 12);
    grid on;
    xlim([0, 1.5]);
    ylim([-4, 8]);
    set(gca, 'FontSize', 14);
    
    % Subplot 2: Contraction factors
    subplot(1, 2, 2);
    hold on;
    for case_idx = 1:num_nu_cases
        semilogy(nu_vals, contr_nu{case_idx}, '-', 'LineWidth', 2, ...
                 'DisplayName', sprintf('\\gamma=%.2f', gamma_fixed_vals(case_idx)), ...
                 'Color', colors(case_idx, :));
    end
    xlabel('\nu', 'FontSize', 16);
    ylabel('Global Contraction Factor', 'FontSize', 16);
    legend('Location', 'NorthEast', 'FontSize', 12);
    grid on;
    ylim([1e-3, 1]);
    set(gca, 'FontSize', 14);
end

%% Helper function for joint view plotting with boundaries
function plot_joint_view_figures(gamma_vals, nu_fixed_vals, opt_gamma, contr_gamma, ...
                                nu_vals, gamma_fixed_vals, opt_nu, contr_nu, ...
                                gamma_fine, nu_fine, contraction_fine, tightness)
    
    num_gamma_cases = length(nu_fixed_vals);
    num_nu_cases = length(gamma_fixed_vals);
    
    % Gamma-sweep joint view
    G = gamma_vals(:)';
    V = nu_fixed_vals(:)';
    Ng = length(G);
    Nv = length(V);
    
    % Extract data for gamma cases
    rho_gam = zeros(Nv, Ng);
    p_gam = zeros(Nv, Ng);
    q_gam = zeros(Nv, Ng);
    
    for case_idx = 1:Nv
        rho_gam(case_idx, :) = contr_gamma{case_idx}(:)';
        p_gam(case_idx, :) = opt_gamma{case_idx}(:, 1)';
        q_gam(case_idx, :) = opt_gamma{case_idx}(:, 2)';
    end
    
    figure('Name', 'Gamma and Nu Sweep Joint View', 'Position', [100, 100, 1200, 520]);
    
    % Left: (p,q) sweep lines with boundary - USING GRAYSCALE
    subplot(1, 2, 1);
    hold on;
    
    % Create grayscale colors from light to dark
    gray_levels = 1 - linspace(0.2, 1, Nv); % Light to dark gray
    
    for case_idx = 1:Nv
        gray_color = gray_levels(case_idx) * [1 1 1]; % RGB gray
        
        plot(p_gam(case_idx, :), q_gam(case_idx, :), '-', 'LineWidth', 2, ...
            'Color', gray_color, 'DisplayName', sprintf('\\nu=%.2f', V(case_idx)));
        
        % Mark start and end points
        plot(p_gam(case_idx, 1), q_gam(case_idx, 1), 'o', 'MarkerSize', 6, ...
            'MarkerFaceColor', gray_color, 'MarkerEdgeColor', 'k', 'HandleVisibility', 'off');
        plot(p_gam(case_idx, end), q_gam(case_idx, end), 's', 'MarkerSize', 6, ...
            'MarkerFaceColor', gray_color, 'MarkerEdgeColor', 'k', 'HandleVisibility', 'off');
    end
    
    % Calculate and plot boundary
    P = p_gam(:);
    Q = q_gam(:);
    valid_points = isfinite(P) & isfinite(Q);
    P = P(valid_points);
    Q = Q(valid_points);
    
    if ~isempty(P)
        try
            % Try alphaShape first (more sophisticated)
            alpha_shape = alphaShape(P, Q);
            alpha_shape.Alpha = 0.85 * criticalAlpha(alpha_shape);
            plot(alpha_shape, 'FaceColor', 'none', 'EdgeColor', 'k', 'LineWidth', 1.4, 'HandleVisibility', 'off');
        catch
            % Fallback to boundary function
            boundary_indices = boundary(P, Q, tightness);
            plot(P(boundary_indices), Q(boundary_indices), 'r--', 'LineWidth', 1.4, 'HandleVisibility', 'off');
        end
    end
    
    xlabel('p', 'FontSize', 14);
    ylabel('q', 'FontSize', 14);
    title('Optimal (p,q) — Varying \gamma (Fixed \nu)', 'FontSize', 12);
    grid on;
    box on;
    xlim([0, 1.5]);
    ylim([-4, 8]);
    set(gca, 'FontSize', 12);
    legend('Location', 'northeast', 'FontSize', 10);
    
    % Nu-sweep joint view - USING GRAYSCALE
    Gf = gamma_fixed_vals(:)';
    Vn = nu_vals(:)';
    Ngf = length(Gf);
    Nn = length(Vn);
    
    % Extract data for nu cases
    rho_nu = zeros(Ngf, Nn);
    p_nu = zeros(Ngf, Nn);
    q_nu = zeros(Ngf, Nn);
    
    for case_idx = 1:Ngf
        rho_nu(case_idx, :) = contr_nu{case_idx}(:)';
        p_nu(case_idx, :) = opt_nu{case_idx}(:, 1)';
        q_nu(case_idx, :) = opt_nu{case_idx}(:, 2)';
    end
    
    % Right: (p,q) sweep lines with boundary - USING GRAYSCALE
    subplot(1, 2, 2);
    hold on;
    
    % Create grayscale colors from light to dark
    gray_levels_nu = 1 - linspace(0.2, 1, Ngf); % Light to dark gray
    
    for case_idx = 1:Ngf
        gray_color = gray_levels_nu(case_idx) * [1 1 1]; % RGB gray
        
        plot(p_nu(case_idx, :), q_nu(case_idx, :), '-', 'LineWidth', 2, ...
            'Color', gray_color, 'DisplayName', sprintf('\\gamma=%.2f', Gf(case_idx)));
        
        % Mark start and end points
        plot(p_nu(case_idx, 1), q_nu(case_idx, 1), 'o', 'MarkerSize', 6, ...
            'MarkerFaceColor', gray_color, 'MarkerEdgeColor', 'k', 'HandleVisibility', 'off');
        plot(p_nu(case_idx, end), q_nu(case_idx, end), 's', 'MarkerSize', 6, ...
            'MarkerFaceColor', gray_color, 'MarkerEdgeColor', 'k', 'HandleVisibility', 'off');
    end
    
    % Calculate and plot boundary
    P = p_nu(:);
    Q = q_nu(:);
    valid_points = isfinite(P) & isfinite(Q);
    P = P(valid_points);
    Q = Q(valid_points);
    
    if ~isempty(P)
        try
            % Try alphaShape first
            alpha_shape = alphaShape(P, Q);
            alpha_shape.Alpha = 0.85 * criticalAlpha(alpha_shape);
            plot(alpha_shape, 'FaceColor', 'none', 'EdgeColor', 'k', 'LineWidth', 1.4, 'HandleVisibility', 'off');
        catch
            % Fallback to boundary function
            boundary_indices = boundary(P, Q, tightness);
            plot(P(boundary_indices), Q(boundary_indices), 'r--', 'LineWidth', 1.4, 'HandleVisibility', 'off');
        end
    end
    
    xlabel('p', 'FontSize', 14);
    ylabel('q', 'FontSize', 14);
    title('Optimal (p,q) — Varying \nu (Fixed \gamma)', 'FontSize', 12);
    grid on;
    box on;
    xlim([0, 1.5]);
    ylim([-4, 8]);
    set(gca, 'FontSize', 12);
    legend('Location', 'northeast', 'FontSize', 10);
    
    % Contour plot - USING LOG-LOG SCALE
    figure('Name', 'Fine Grid Contraction Factor - Log Scale', 'Position', [100, 100, 800, 600]);
    
    % Create contour plot with log-log scale
    contourf(gamma_fine, nu_fine, log10(contraction_fine), 20, 'LineColor', 'none');
    colormap(parula);
    colorbar;
    
    % Set both axes to log scale
    set(gca, 'XScale', 'log', 'YScale', 'log');
    
    xlabel('\gamma', 'FontSize', 14);
    ylabel('\nu', 'FontSize', 14);
    title('Optimal Contraction Factor - Log-Log-Log Scale', 'FontSize', 12);
    grid on;
    box on;
    set(gca, 'FontSize', 12);
    
    % Add grid lines for better readability in log scale
    set(gca, 'XMinorGrid', 'on', 'YMinorGrid', 'on');
    
    % Set appropriate axis limits for log scale
    xlim([1e-2, 1e1]);
    ylim([1e-2, 1e1]);
end