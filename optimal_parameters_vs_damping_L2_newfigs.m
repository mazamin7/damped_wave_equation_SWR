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
        objfun = @(x) obj_L2(N, T, dt, J, c, gamma, nu, a, M, x(1), x(2), ky);
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
% --- DRASTICALLY INCREASED SIZES FOR LATEX VISIBILITY ---
FS_AXIS  = 30;  % Tick numbers (was 22)
FS_LABEL = 36;  % x/y labels (was 24)
LW_BOLD  = 4.0; % Line width (was 3.0)
MS_DOTS  = 10; % Marker size (was 60)

% --- EXPLICIT LAYOUT DEFINITIONS (Normalized 0 to 1) ---
% We need large margins to fit the huge text.
% 1. PLOT: Left=0.18 (for Y-label), Bottom=0.18, Width=0.52 (Narrower)
POS_AX  = [0.18, 0.18, 0.52, 0.75]; 

% 2. COLORBAR: Start=0.78 (Gap=0.08), Width=0.04
POS_CB  = [0.78, 0.18, 0.04, 0.75]; 

%% Figure 1: Envelope and Scatter
figure('Name','Envelope','Position',[100 100 900 650], 'Color', 'w');
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
set(gca,'FontSize',FS_AXIS, 'LineWidth', 3, 'FontWeight', 'bold');
xlim([p_min - px, p_max + px]);
ylim([q_min - qx, q_max + qx]);

% Apply layout (Creates empty space on right to match other figs)
set(gca, 'Position', POS_AX); 


%% Figure 2: Isolines (Fixed Nu)
nu_levels = logspace(-1, 0, 10);
nu_idx = zeros(size(nu_levels));
for k = 1:numel(nu_levels)
    [~, nu_idx(k)] = min(abs(nu_vals - nu_levels(k)));
end
nu_idx = unique(nu_idx);

figure('Name','Isolines: Fixed Nu','Position',[150 150 900 650], 'Color', 'w');
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
set(gca,'FontSize',FS_AXIS, 'LineWidth', 3, 'FontWeight', 'bold');
xlim([p_min - px, p_max + px]);
ylim([q_min - qx, q_max + qx]);

% --- LAYOUT & COLORBAR ---
set(gca, 'Position', POS_AX);
c = colorbar; colormap(cmap); c.Position = POS_CB;

c.Label.String = '\nu'; 
c.Label.FontSize = FS_LABEL; 
c.Label.Rotation = 0; 
c.Label.Units = 'normalized';
% Move label further left (-2.0) to clear the larger tick numbers
c.Label.Position = [-1.0, 0.5, 0]; 
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

figure('Name','Isolines: Fixed Gamma','Position',[200 200 900 650], 'Color', 'w');
set(gcf, 'PaperPositionMode', 'auto');
hold on; box on; grid on;

plot(P_boundary, Q_boundary, 'k-', 'LineWidth', LW_BOLD);
for k = 1:numel(gamma_idx)
    i = gamma_idx(k);
    plot(p_opt(i,:), q_opt(i,:), '-', 'LineWidth', LW_BOLD, 'Color', cmap(k,:));
end

xlabel('p','FontSize',FS_LABEL, 'FontWeight', 'bold');
ylabel('q','FontSize',FS_LABEL, 'FontWeight', 'bold');
set(gca,'FontSize',FS_AXIS, 'LineWidth', 3, 'FontWeight', 'bold');
xlim([p_min - px, p_max + px]);
ylim([q_min - qx, q_max + qx]);

% --- LAYOUT & COLORBAR ---
set(gca, 'Position', POS_AX); 
c = colorbar; colormap(cmap); c.Position = POS_CB;

c.Label.String = '\gamma'; 
c.Label.FontSize = FS_LABEL;
c.Label.Rotation = 0; 
c.Label.Units = 'normalized';
c.Label.Position = [-1.0, 0.5, 0]; 
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

figure('Name','Contraction p-q','Position',[250 250 900 650], 'Color', 'w');
set(gcf, 'PaperPositionMode', 'auto');
hold on; box on; grid on;

contourf(Pgrid, Qgrid, Rgrid, 20, 'LineStyle','none');
plot(P_boundary, Q_boundary, 'k-', 'LineWidth', LW_BOLD);

xlabel('p','FontSize',FS_LABEL, 'FontWeight', 'bold');
ylabel('q','FontSize',FS_LABEL, 'FontWeight', 'bold');
set(gca,'FontSize',FS_AXIS, 'LineWidth', 3, 'FontWeight', 'bold');
xlim([p_min - px, p_max + px]);
ylim([q_min - qx, q_max + qx]);

% --- LAYOUT & COLORBAR ---
set(gca, 'Position', POS_AX);
c = colorbar; colormap(parula); c.Position = POS_CB;

c.Label.String = '\rho';
c.Label.FontSize = FS_LABEL;
c.Label.Rotation = 0; 
c.Label.Units = 'normalized';
c.Label.Position = [-1.0, 0.5, 0]; 
c.Label.VerticalAlignment = 'middle';

% Force 3 decimal places
rng = caxis;
ticks = linspace(rng(1), rng(2), 6); 
c.Ticks = ticks;
c.TickLabels = num2str(ticks', '%.3f');


%% Figure 5: Contraction in Parameter Space
figure('Name','Contraction gamma-nu','Position',[300 300 900 650], 'Color', 'w');
set(gcf, 'PaperPositionMode', 'auto');
hold on; box on; grid on;

[GammaGrid, NuGrid] = meshgrid(gamma_vals, nu_vals);
contourf(GammaGrid, NuGrid, rho_opt', 30, 'LineStyle', 'none');

% Create the label and capture the handle 'h'
h = xlabel('\gamma', 'FontSize', FS_LABEL, 'FontWeight', 'bold');

% Switch units to normalized so (0,0) is bottom-left and (1,1) is top-right of the axis
set(h, 'Units', 'normalized');

% Set position: [x, y, z]
% x = 0.5   (Center horizontally)
% y = -0.08 (Vertical position. Standard is approx -0.12. 
%            Increase this number to move it UP towards the axis line.)
set(h, 'Position', [0.5, -0.11, 0]);

ylabel('\nu',    'FontSize', FS_LABEL, 'FontWeight', 'bold');
set(gca, 'FontSize', FS_AXIS, 'LineWidth', 3, 'FontWeight', 'bold');
set(gca, 'XScale', 'log', 'YScale', 'log');

% --- LAYOUT & COLORBAR ---
set(gca, 'Position', POS_AX);
c = colorbar; colormap(parula); c.Position = POS_CB;

c.Label.String = '\rho';
c.Label.FontSize = FS_LABEL;
c.Label.Rotation = 0; 
c.Label.Units = 'normalized';
c.Label.Position = [-1.0, 0.5, 0]; 
c.Label.VerticalAlignment = 'middle';

% Force 3 decimal places
rng = caxis;
ticks = linspace(rng(1), rng(2), 6);
c.Ticks = ticks;
c.TickLabels = num2str(ticks', '%.3f');
