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
J  = 50;        % Number of frequency steps

% Optimization parameters
tol = 1e-9;

padding = 0.1;

%% Fine grid in (gamma, nu)

% Gamma in (0.1, 10), nu in (0.1, 1)
Ng = 81;                        % number of gamma samples
Nn = 81;                        % number of nu samples
gamma_vals = logspace(-1, 1, Ng);   % 10^{-1} .. 10^{1}
nu_vals    = logspace(-1, 0, Nn);   % 10^{-1} .. 10^{0}

% Storage for optimal (p,q) and contraction
p_opt   = zeros(Ng, Nn);
q_opt   = zeros(Ng, Nn);
rho_opt = zeros(Ng, Nn);

% Optimization settings
optim_options = optimset('Display', 'off', ...
                         'TolX', tol, 'TolFun', tol, ...
                         'MaxFunEvals', 2e4, 'MaxIter', 2e4);

x_fallback = [1/c, 0];

fprintf('=== Computing optimal (p,q) on fine (gamma,nu) grid ===\n');
for j = 1:Nn
    nu = nu_vals(j);
    fprintf('  Row %d/%d (nu = %.3f)\n', j, Nn, nu);
    for i = 1:Ng
        gamma = gamma_vals(i);

        % Warm-start: use neighbor in gamma or nu if available
        if i > 1
            x0 = [p_opt(i-1,j), q_opt(i-1,j)];
        elseif j > 1
            x0 = [p_opt(i,j-1), q_opt(i,j-1)];
        else
            x0 = x_fallback;
        end

        if any(~isfinite(x0)), x0 = x_fallback; end

        objfun = @(x) obj_Linf(N, T, dt, J, c, gamma, nu, a, M, x(1), x(2), ky);
        [x_opt, fval] = fminsearch(objfun, x0, optim_options);

        if ~isfinite(fval) || any(~isfinite(x_opt))
            x_opt = x_fallback;
            fval  = objfun(x_opt);
        end

        p_opt(i,j)   = x_opt(1);
        q_opt(i,j)   = x_opt(2);
        rho_opt(i,j) = fval;
    end
end

%% Build 4-sided boundary in (p,q) space (envelope)

[~, i_gmin] = min(gamma_vals);
[~, i_gmax] = max(gamma_vals);
[~, j_nmin] = min(nu_vals);
[~, j_nmax] = max(nu_vals);

% 4 boundary curves in parameter space mapped to (p,q):
% gamma = gamma_min, nu from min..max
P_gmin = p_opt(i_gmin,:);  Q_gmin = q_opt(i_gmin,:);
% gamma = gamma_max, nu from max..min (reverse to close polygon)
P_gmax = fliplr(p_opt(i_gmax,:));  Q_gmax = fliplr(q_opt(i_gmax,:));
% nu = nu_min, gamma from max..min
P_nmin = p_opt(:,j_nmin);  Q_nmin = q_opt(:,j_nmin);
P_nmin = flipud(P_nmin);   Q_nmin = flipud(Q_nmin);
% nu = nu_max, gamma from min..max
P_nmax = p_opt(:,j_nmax);  Q_nmax = q_opt(:,j_nmax);

% Concatenate in correct order to form closed polygon:
P_boundary = [P_gmin,  P_nmax',  P_gmax,  P_nmin'];
Q_boundary = [Q_gmin,  Q_nmax',  Q_gmax,  Q_nmin'];

% Use this to define plotting limits (envelope box)
p_min = min(P_boundary); p_max = max(P_boundary);
q_min = min(Q_boundary); q_max = max(Q_boundary);

%% Figure 1: envelope and all points (optional) + boundary

figure('Name','Envelope in (p,q) space','Position',[100 100 800 600]);
hold on; box on; grid on;

valid = isfinite(p_opt) & isfinite(q_opt);
scatter(p_opt(valid), q_opt(valid), 10, [0.8 0.8 0.8], 'filled', ...
        'DisplayName','Optimal (p,q) for all (\gamma,\nu)');

plot(P_boundary, Q_boundary, 'k-', 'LineWidth', 2, 'HandleVisibility', 'off');%, 'DisplayName','Envelope boundary');

xlabel('p','FontSize',16);
ylabel('q','FontSize',16);
set(gca,'FontSize',16);
legend('Location','best');

% Add padding of 5% of interval
px = padding * (p_max - p_min);
qx = padding * (q_max - q_min);

xlim([p_min - px, p_max + px]);
ylim([q_min - qx, q_max + qx]);

%% Figure 2: isolines for constant nu (γ varying) inside same boundary

% Choose nu-levels (logarithmic)
nu_levels = logspace(-1, 0, 10);
nu_idx = zeros(size(nu_levels));
for k = 1:numel(nu_levels)
    [~, nu_idx(k)] = min(abs(nu_vals - nu_levels(k)));
end
nu_idx = unique(nu_idx);

figure('Name','Isolines: constant \nu in (p,q) space','Position',[150 150 800 600]);
hold on; box on; grid on;

% Draw boundary again for context
plot(P_boundary, Q_boundary, 'k-', 'LineWidth', 2, 'HandleVisibility', 'off');%, 'DisplayName','Envelope boundary');

colors = lines(numel(nu_idx));
for k = 1:numel(nu_idx)
    j = nu_idx(k);
    p_curve = p_opt(:,j);
    q_curve = q_opt(:,j);
    mask = isfinite(p_curve) & isfinite(q_curve);
    plot(p_curve(mask), q_curve(mask), '-', 'LineWidth', 2, ...
         'Color', colors(k,:), ...
         'DisplayName', sprintf('\\nu \\approx %.3f', nu_vals(j)));
end

xlabel('p','FontSize',16);
ylabel('q','FontSize',16);
set(gca,'FontSize',16);
legend('Location','best');

% Add padding of 5% of interval
px = padding * (p_max - p_min);
qx = padding * (q_max - q_min);

xlim([p_min - px, p_max + px]);
ylim([q_min - qx, q_max + qx]);

%% Figure 3: isolines for constant gamma (ν varying) inside same boundary

gamma_levels = logspace(-1, 1, 10);
gamma_idx = zeros(size(gamma_levels));
for k = 1:numel(gamma_levels)
    [~, gamma_idx(k)] = min(abs(gamma_vals - gamma_levels(k)));
end
gamma_idx = unique(gamma_idx);

figure('Name','Isolines: constant \gamma in (p,q) space','Position',[200 200 800 600]);
hold on; box on; grid on;

% Draw boundary again for context
plot(P_boundary, Q_boundary, 'k-', 'LineWidth', 2, 'HandleVisibility', 'off');%, 'DisplayName','Envelope boundary');

colors = lines(numel(gamma_idx));
for k = 1:numel(gamma_idx)
    i = gamma_idx(k);
    p_curve = p_opt(i,:);
    q_curve = q_opt(i,:);
    mask = isfinite(p_curve) & isfinite(q_curve);
    plot(p_curve(mask), q_curve(mask), '-', 'LineWidth', 2, ...
         'Color', colors(k,:), ...
         'DisplayName', sprintf('\\gamma \\approx %.3f', gamma_vals(i)));
end

xlabel('p','FontSize',16);
ylabel('q','FontSize',16);
set(gca,'FontSize',16);
legend('Location','best');

% Add padding of 5% of interval
px = padding * (p_max - p_min);
qx = padding * (q_max - q_min);

xlim([p_min - px, p_max + px]);
ylim([q_min - qx, q_max + qx]);

%% Figure 4: contour of global contraction factor in (p,q) inside boundary

% Scattered data: (p_opt, q_opt) -> rho_opt
mask_scatt = isfinite(p_opt) & isfinite(q_opt) & isfinite(rho_opt);
p_scatt = p_opt(mask_scatt);
q_scatt = q_opt(mask_scatt);
rho_scatt = rho_opt(mask_scatt);

% Interpolant over (p,q)
F_rho = scatteredInterpolant(p_scatt, q_scatt, rho_scatt, ...
                             'natural', 'none');

% Regular grid in (p,q) over the envelope box
Np = 200;
Nq = 200;
[Pgrid, Qgrid] = meshgrid(linspace(p_min, p_max, Np), ...
                          linspace(q_min, q_max, Nq));

Rgrid = F_rho(Pgrid, Qgrid);

% Mask outside the polygon boundary
[in_poly, on_poly] = inpolygon(Pgrid, Qgrid, P_boundary, Q_boundary);
inside = in_poly | on_poly;
Rgrid(~inside) = NaN;

figure('Name','Contraction factor in (p,q) space','Position',[250 250 800 600]);
hold on; box on; grid on;

% Contour plot of the global contraction factor
contourf(Pgrid, Qgrid, Rgrid, 20, 'LineStyle','none');
colormap(parula);
colorbar;

% Draw boundary
plot(P_boundary, Q_boundary, 'k-', 'LineWidth', 2, 'HandleVisibility', 'off');

xlabel('p','FontSize',16);
ylabel('q','FontSize',16);
set(gca,'FontSize',16);
% Add padding of 5% of interval
px = padding * (p_max - p_min);
qx = padding * (q_max - q_min);

xlim([p_min - px, p_max + px]);
ylim([q_min - qx, q_max + qx]);

%% Figure 5: global contraction factor rho(gamma,nu) in parameter space

figure('Name','Contraction factor in (\gamma,\nu) space','Position',[300 300 800 600]);
hold on; box on; grid on;

% We already have rho_opt(i,j) defined on the grid:
%   gamma_vals(i), nu_vals(j)

% Create grids for plotting
[GammaGrid, NuGrid] = meshgrid(gamma_vals, nu_vals);

% Note: rho_opt(i,j) needs to be transposed because meshgrid produces
%       GammaGrid(j,i), NuGrid(j,i)
R = rho_opt';   

% Option A: contour plot (recommended)
contourf(GammaGrid, NuGrid, R, 30, 'LineStyle', 'none');
colormap(parula);
colorbar;

set(gca, 'XScale', 'log', 'YScale', 'log');   % Log scale for both axes
xlabel('\gamma', 'FontSize', 16);
ylabel('\nu',    'FontSize', 16);
% title('Global contraction factor \rho(\gamma,\nu)', 'FontSize', 16);
set(gca, 'FontSize', 16);
