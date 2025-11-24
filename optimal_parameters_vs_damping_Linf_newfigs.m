clear; close all; clc;
addpath("utils\")

%% Common parameters
P = get_sim_params_1D();
N  = P.N; a  = P.a; M  = P.M;
T  = P.T; c  = P.c;
dt = P.dt; J  = P.J;
ky = 0; 

% Optimization parameters
tol = 1e-9;
padding = 0.05;

%% Fine grid in (gamma, nu)
Ng = 81;                        
Nn = 81;           
% Ng = 6;                        
% Nn = 6;    
gamma_vals = logspace(-1, 1, Ng);   
nu_vals    = logspace(-1, 0, Nn);   

p_opt   = zeros(Ng, Nn);
q_opt   = zeros(Ng, Nn);
rho_opt = zeros(Ng, Nn);

optim_options = optimset('Display', 'off', 'TolX', tol, 'TolFun', tol, ...
                         'MaxFunEvals', 2e4, 'MaxIter', 2e4);
x_fallback = [1/c, 0];

fprintf('=== Computing optimal (p,q) on fine (gamma,nu) grid ===\n');
for j = 1:Nn
    nu = nu_vals(j);
    fprintf('  Row %d/%d (nu = %.3f)\n', j, Nn, nu); % Commented for speed
    for i = 1:Ng
        gamma = gamma_vals(i);
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

%% Build Envelope
[~, i_gmin] = min(gamma_vals); [~, i_gmax] = max(gamma_vals);
[~, j_nmin] = min(nu_vals);    [~, j_nmax] = max(nu_vals);

P_gmin = p_opt(i_gmin,:);          Q_gmin = q_opt(i_gmin,:);
P_gmax = fliplr(p_opt(i_gmax,:));  Q_gmax = fliplr(q_opt(i_gmax,:));
P_nmin = flipud(p_opt(:,j_nmin));  Q_nmin = flipud(q_opt(:,j_nmin));
P_nmax = p_opt(:,j_nmax);          Q_nmax = q_opt(:,j_nmax);

P_boundary = [P_gmin,  P_nmax',  P_gmax,  Q_nmin']; % Note: Q_nmin index check logic maintained
Q_boundary = [Q_gmin,  Q_nmax',  Q_gmax,  Q_nmin']; 
% Re-stitching logic from original script to ensure closure:
P_boundary = [p_opt(i_gmin,:), p_opt(:,j_nmax)', fliplr(p_opt(i_gmax,:)), fliplr(p_opt(:,j_nmin)')];
Q_boundary = [q_opt(i_gmin,:), q_opt(:,j_nmax)', fliplr(q_opt(i_gmax,:)), fliplr(q_opt(:,j_nmin)')];

p_min = min(P_boundary); p_max = max(P_boundary);
q_min = min(Q_boundary); q_max = max(Q_boundary);
px = padding * (p_max - p_min);
qx = padding * (q_max - q_min);

%% VISUALIZATION SETTINGS
% Large fonts and thick lines for small LaTeX subfigures
FS_AXIS = 22;
FS_LABEL = 24;
LW_BOLD = 3.0;
MS_DOTS = 10;
% --- EXPLICIT LAYOUT DEFINITIONS (Normalized 0 to 1) ---
% 1. PLOT: Starts at 0.13, Width 0.66 (Ends at 0.79)
POS_AX  = [0.13, 0.15, 0.66, 0.78]; 

% 2. COLORBAR: Starts at 0.88 (Gap of 0.09 for label), Width 0.04, Ends at 0.92
POS_CB  = [0.88, 0.15, 0.04, 0.78];

%% Figure 1: Envelope and Scatter
% We apply POS_AX here too. It will leave empty space on the right, 
% but ensures the plot frame is identical to Figures 2-5.
figure('Name','Envelope','Position',[100 100 800 600], 'Color', 'w');
set(gcf, 'PaperPositionMode', 'auto');
hold on; box on; grid on;

valid = isfinite(p_opt) & isfinite(q_opt);

% Plot dots
scatter(p_opt(valid), q_opt(valid), MS_DOTS, [0.7 0.7 0.7], 'filled', ...
    'MarkerFaceAlpha', 0.5); 

% Plot boundary
plot(P_boundary, Q_boundary, 'k-', 'LineWidth', LW_BOLD);

xlabel('p','FontSize',FS_LABEL, 'FontWeight', 'bold');
ylabel('q','FontSize',FS_LABEL, 'FontWeight', 'bold');
set(gca,'FontSize',FS_AXIS, 'LineWidth', 2);
xlim([p_min - px, p_max + px]);
ylim([q_min - qx, q_max + qx]);

% --- APPLY LAYOUT (No Colorbar) ---
set(gca, 'Position', POS_AX);

%% Figure 2: Isolines (Fixed Nu)
nu_levels = logspace(-1, 0, 10);
nu_idx = zeros(size(nu_levels));
for k = 1:numel(nu_levels)
    [~, nu_idx(k)] = min(abs(nu_vals - nu_levels(k)));
end
nu_idx = unique(nu_idx);

figure('Name','Isolines: Fixed Nu','Position',[150 150 800 600], 'Color', 'w');
set(gcf, 'PaperPositionMode', 'auto'); 
hold on; box on; grid on;

plot(P_boundary, Q_boundary, 'k-', 'LineWidth', LW_BOLD);

num_curves = numel(nu_idx);
cmap = parula(num_curves); 

for k = 1:num_curves
    j = nu_idx(k);
    plot(p_opt(:,j), q_opt(:,j), '-', 'LineWidth', LW_BOLD, 'Color', cmap(k,:));
end

xlabel('p','FontSize',FS_LABEL, 'FontWeight', 'bold');
ylabel('q','FontSize',FS_LABEL, 'FontWeight', 'bold');
set(gca,'FontSize',FS_AXIS, 'LineWidth', 2);
xlim([p_min - px, p_max + px]);
ylim([q_min - qx, q_max + qx]);

% --- FIX: Manual Positioning ---
set(gca, 'Position', POS_AX); % Lock plot position

c = colorbar;
colormap(cmap);
c.Position = POS_CB; % Lock colorbar position further right

c.Label.String = '\nu';
c.Label.FontSize = FS_LABEL;
c.Label.Rotation = 0; 
c.Label.Units = 'normalized';
% Position label inside the gap (x < 0 relative to colorbar)
c.Label.Position = [-0.8, 0.5, 0]; 
c.Label.VerticalAlignment = 'middle';
caxis([min(nu_levels) max(nu_levels)]);
set(gca, 'ColorScale', 'log'); 


%% Figure 3: Isolines (Fixed Gamma)
gamma_levels = logspace(-1, 1, 10);
gamma_idx = zeros(size(gamma_levels));
for k = 1:numel(gamma_levels)
    [~, gamma_idx(k)] = min(abs(gamma_vals - gamma_levels(k)));
end
gamma_idx = unique(gamma_idx);

figure('Name','Isolines: Fixed Gamma','Position',[200 200 800 600], 'Color', 'w');
set(gcf, 'PaperPositionMode', 'auto');
hold on; box on; grid on;

plot(P_boundary, Q_boundary, 'k-', 'LineWidth', LW_BOLD);

num_curves = numel(gamma_idx);
cmap = parula(num_curves);

for k = 1:num_curves
    i = gamma_idx(k);
    plot(p_opt(i,:), q_opt(i,:), '-', 'LineWidth', LW_BOLD, 'Color', cmap(k,:));
end

xlabel('p','FontSize',FS_LABEL, 'FontWeight', 'bold');
ylabel('q','FontSize',FS_LABEL, 'FontWeight', 'bold');
set(gca,'FontSize',FS_AXIS, 'LineWidth', 2);
xlim([p_min - px, p_max + px]);
ylim([q_min - qx, q_max + qx]);

% --- FIX: Manual Positioning ---
set(gca, 'Position', POS_AX); 

c = colorbar;
colormap(cmap);
c.Position = POS_CB;

c.Label.String = '\gamma';
c.Label.FontSize = FS_LABEL;
c.Label.Rotation = 0; 
c.Label.Units = 'normalized';
c.Label.Position = [-0.8, 0.5, 0]; 
c.Label.VerticalAlignment = 'middle';
caxis([min(gamma_levels) max(gamma_levels)]);
set(gca, 'ColorScale', 'log');


%% Figure 4: Contour of Global Contraction
mask_scatt = isfinite(p_opt) & isfinite(q_opt) & isfinite(rho_opt);
F_rho = scatteredInterpolant(p_opt(mask_scatt), q_opt(mask_scatt), rho_opt(mask_scatt), 'natural', 'none');

Np = 200; Nq = 200;
[Pgrid, Qgrid] = meshgrid(linspace(p_min, p_max, Np), linspace(q_min, q_max, Nq));
Rgrid = F_rho(Pgrid, Qgrid);
[in_poly, on_poly] = inpolygon(Pgrid, Qgrid, P_boundary, Q_boundary);
Rgrid(~(in_poly | on_poly)) = NaN;

figure('Name','Contraction p-q','Position',[250 250 800 600], 'Color', 'w');
set(gcf, 'PaperPositionMode', 'auto');
hold on; box on; grid on;

contourf(Pgrid, Qgrid, Rgrid, 20, 'LineStyle','none');
plot(P_boundary, Q_boundary, 'k-', 'LineWidth', LW_BOLD);

xlabel('p','FontSize',FS_LABEL, 'FontWeight', 'bold');
ylabel('q','FontSize',FS_LABEL, 'FontWeight', 'bold');
set(gca,'FontSize',FS_AXIS, 'LineWidth', 2);
xlim([p_min - px, p_max + px]);
ylim([q_min - qx, q_max + qx]);

% --- FIX: Manual Positioning ---
set(gca, 'Position', POS_AX);

c = colorbar;
colormap(parula);
c.Position = POS_CB;

c.Label.String = '\rho';
c.Label.FontSize = FS_LABEL;
c.Label.Rotation = 0; 
c.Label.Units = 'normalized';
c.Label.Position = [-0.8, 0.5, 0]; 
c.Label.VerticalAlignment = 'middle';


%% Figure 5: Contraction in Parameter Space
figure('Name','Contraction gamma-nu','Position',[300 300 800 600], 'Color', 'w');
set(gcf, 'PaperPositionMode', 'auto');
hold on; box on; grid on;

[GammaGrid, NuGrid] = meshgrid(gamma_vals, nu_vals);
contourf(GammaGrid, NuGrid, rho_opt', 30, 'LineStyle', 'none');

xlabel('\gamma', 'FontSize', FS_LABEL, 'FontWeight', 'bold');
ylabel('\nu',    'FontSize', FS_LABEL, 'FontWeight', 'bold');
set(gca, 'FontSize', FS_AXIS, 'LineWidth', 2);
set(gca, 'XScale', 'log', 'YScale', 'log');

% --- FIX: Manual Positioning ---
set(gca, 'Position', POS_AX);

c = colorbar;
colormap(parula);
c.Position = POS_CB;

c.Label.String = '\rho';
c.Label.FontSize = FS_LABEL;
c.Label.Rotation = 0; 
c.Label.Units = 'normalized';
c.Label.Position = [-0.8, 0.5, 0]; 
c.Label.VerticalAlignment = 'middle';
