%% CLEAN UP AND SETUP
% Clear workspace, close figures, clear command window
clear; close all; clc;

% Add utility functions directory to MATLAB path
addpath("utils\")

%% COMMON PARAMETERS AND SIMULATION SETUP
% Load simulation parameters from get_sim_params() function
P = get_sim_params();

% Extract parameters from structure
N  = P.N;      % Number of grid points or system size
a  = P.a;      % System parameter
M  = P.M;      % Another system parameter
T  = P.T;      % Total simulation time
c  = P.c;      % Wave speed or similar constant
dt = P.dt;     % Time step size
J  = P.J;      % Parameter for objective function

% Optimization parameter

%% OPTIMIZATION CONFIGURATION
tol = 1e-9;        % Tolerance for optimization convergence
padding = 0.05;    % Padding for plot axis limits (5% extra space)

%% FINE GRID DEFINITION FOR PARAMETER SPACE EXPLORATION
% Create logarithmic grids for gamma and nu parameters
Ng = 81;        % Number of grid points for gamma
Nn = 81;        % Number of grid points for nu

% Log-spaced parameter ranges:
% gamma from 10^(-1) to 10^1 = [0.1, 10]
% nu from 10^(-1) to 10^0 = [0.1, 1]
gamma_vals = logspace(-1, 1, Ng);   
nu_vals    = logspace(-1, 0, Nn);   

% Preallocate arrays to store optimization results
p_opt   = zeros(Ng, Nn);   % Optimal p values
q_opt   = zeros(Ng, Nn);   % Optimal q values
rho_opt = zeros(Ng, Nn);   % Optimal contraction rates (objective values)

%% OPTIMIZATION SETUP
% Configure MATLAB's optimization options
optim_options = optimset('Display', 'off', ...    % Suppress iteration output
                         'TolX', tol, ...         % Tolerance on parameter changes
                         'TolFun', tol, ...       % Tolerance on function value
                         'MaxFunEvals', 2e4, ...  % Maximum function evaluations
                         'MaxIter', 2e4);         % Maximum iterations

% Fallback initial guess if other initialization methods fail
x_fallback = [1/c, 0];  % Reasonable default: [p0, q0]

%% MAIN OPTIMIZATION LOOP OVER PARAMETER GRID
fprintf('=== Computing optimal (p,q) on fine (gamma,nu) grid ===\n');

% Loop over nu values (outer loop)
for j = 1:Nn
    nu = nu_vals(j);
    fprintf('  Row %d/%d (nu = %.3f)\n', j, Nn, nu);
    
    % Loop over gamma values (inner loop)
    for i = 1:Ng
        gamma = gamma_vals(i);
        
        %% SMART INITIALIZATION STRATEGY
        % Use previously optimized values as initial guess for next point
        % This exploits continuity in parameter space for faster convergence
        if i > 1
            % Use previous gamma value (same nu) as initial guess
            x0 = [p_opt(i-1,j), q_opt(i-1,j)];
        elseif j > 1
            % Use previous nu value (same gamma) as initial guess
            x0 = [p_opt(i,j-1), q_opt(i,j-1)];
        else
            % First point in grid: use fallback initial guess
            x0 = x_fallback;
        end
        
        % Fallback to default if initialization failed
        if any(~isfinite(x0))
            x0 = x_fallback;
        end
        
        %% OPTIMIZATION EXECUTION
        % Define objective function with current parameters
        objfun = @(x) obj_Linf(N, T, dt, J, c, gamma, nu, a, M, x(1), x(2));
        
        % Run Nelder-Mead simplex optimization (fminsearch)
        [x_opt, fval] = fminsearch(objfun, x0, optim_options);
        
        %% ERROR HANDLING FOR FAILED OPTIMIZATIONS
        % If optimization failed (non-finite results), use fallback values
        if ~isfinite(fval) || any(~isfinite(x_opt))
            x_opt = x_fallback;
            fval  = objfun(x_opt);  % Re-evaluate objective at fallback point
        end
        
        %% STORE RESULTS
        p_opt(i,j)   = x_opt(1);    % Store optimal p
        q_opt(i,j)   = x_opt(2);    % Store optimal q
        rho_opt(i,j) = fval;        % Store optimal contraction rate
    end
end

%% BOUNDARY ENVELOPE CONSTRUCTION
% Extract boundary points from the optimized grid to create envelope polygon

% Find indices of minimum and maximum gamma/nu values
[~, i_gmin] = min(gamma_vals);  % Index of smallest gamma
[~, i_gmax] = max(gamma_vals);  % Index of largest gamma
[~, j_nmin] = min(nu_vals);     % Index of smallest nu
[~, j_nmax] = max(nu_vals);     % Index of largest nu

% Extract boundary curves:
% Bottom edge: fixed gamma_min, varying nu
P_gmin = p_opt(i_gmin,:);
Q_gmin = q_opt(i_gmin,:);

% Top edge: fixed gamma_max, varying nu (reversed for proper polygon order)
P_gmax = fliplr(p_opt(i_gmax,:));
Q_gmax = fliplr(q_opt(i_gmax,:));

% Left edge: fixed nu_min, varying gamma (reversed for proper polygon order)
P_nmin = flipud(p_opt(:,j_nmin));
Q_nmin = flipud(q_opt(:,j_nmin));

% Right edge: fixed nu_max, varying gamma
P_nmax = p_opt(:,j_nmax);
Q_nmax = q_opt(:,j_nmax);

% Concatenate boundary points in clockwise order to form closed polygon
P_boundary = [p_opt(i_gmin,:), ...      % Bottom edge (left to right)
              p_opt(:,j_nmax)', ...     % Right edge (bottom to top)
              fliplr(p_opt(i_gmax,:)), ...  % Top edge (right to left)
              fliplr(p_opt(:,j_nmin)')];    % Left edge (top to bottom)

Q_boundary = [q_opt(i_gmin,:), ...      % Bottom edge
              q_opt(:,j_nmax)', ...     % Right edge
              fliplr(q_opt(i_gmax,:)), ...  % Top edge
              fliplr(q_opt(:,j_nmin)')];    % Left edge

% Calculate axis limits with padding for better visualization
p_min = min(P_boundary); p_max = max(P_boundary);
q_min = min(Q_boundary); q_max = max(Q_boundary);
px = padding * (p_max - p_min);  % Horizontal padding
qx = padding * (q_max - q_min);  % Vertical padding

%% VISUALIZATION SETTINGS (POSTER STYLE)
% Define consistent styling for all figures
FS_AXIS  = 26;      % Font size for axis tick labels
FS_LABEL = 28;      % Font size for axis labels
LW_BOLD  = 4.0;     % Line width for bold lines (boundaries)
MS_DOTS  = 10;      % Marker size for scatter plot dots

% Fixed layout positions to prevent automatic resizing
POS_AX  = [0.15, 0.15, 0.64, 0.78];  % Main axis position [left, bottom, width, height]
POS_CB  = [0.84, 0.15, 0.04, 0.78];  % Colorbar position

%% FIGURE 1: ENVELOPE AND SCATTER PLOT OF OPTIMAL POINTS
% Shows all optimal (p,q) points with boundary envelope
figure('Name','Envelope','Position',[100 100 800 700], 'Color', 'w');
set(gcf, 'PaperPositionMode', 'auto');  % Maintain size when printing
hold on; box on; grid on;

% Identify valid (finite) data points
valid = isfinite(p_opt) & isfinite(q_opt);

% Scatter plot of all optimal points (gray, semi-transparent)
scatter(p_opt(valid), q_opt(valid), MS_DOTS, [0.7 0.7 0.7], 'filled', ...
    'MarkerFaceAlpha', 0.5); 

% Plot boundary envelope as bold black line
plot(P_boundary, Q_boundary, 'k-', 'LineWidth', LW_BOLD);

% Labels and formatting
xlabel('p','FontSize',FS_LABEL, 'FontWeight', 'bold');
ylabel('q','FontSize',FS_LABEL, 'FontWeight', 'bold');
set(gca,'FontSize',FS_AXIS, 'LineWidth', 2, 'FontWeight', 'bold');
xlim([p_min - px, p_max + px]);
ylim([q_min - qx, q_max + qx]);

% Apply fixed layout
set(gca, 'Position', POS_AX);

%% FIGURE 2: ISOLINES AT FIXED NU VALUES
% Shows curves of optimal (p,q) for constant nu values
nu_levels = logspace(-1, 0, 10);  % 10 log-spaced nu levels

% Find indices in nu_vals closest to the desired levels
nu_idx = zeros(size(nu_levels));
for k = 1:numel(nu_levels)
    [~, nu_idx(k)] = min(abs(nu_vals - nu_levels(k)));
end
nu_idx = unique(nu_idx);  % Remove duplicates

figure('Name','Isolines: Fixed Nu','Position',[150 150 800 700], 'Color', 'w');
set(gcf, 'PaperPositionMode', 'auto'); 
hold on; box on; grid on;

% Plot boundary envelope as reference
plot(P_boundary, Q_boundary, 'k-', 'LineWidth', LW_BOLD);

% Color map for different nu levels
num_curves = numel(nu_idx);
cmap = parula(num_curves); 

% Plot each constant-nu curve
for k = 1:num_curves
    j = nu_idx(k);
    plot(p_opt(:,j), q_opt(:,j), '-', 'LineWidth', LW_BOLD, 'Color', cmap(k,:));
end

% Labels and formatting
xlabel('p','FontSize',FS_LABEL, 'FontWeight', 'bold');
ylabel('q','FontSize',FS_LABEL, 'FontWeight', 'bold');
set(gca,'FontSize',FS_AXIS, 'LineWidth', 2, 'FontWeight', 'bold');
xlim([p_min - px, p_max + px]);
ylim([q_min - qx, q_max + qx]);

% Apply fixed layout
set(gca, 'Position', POS_AX); 

% Add colorbar for nu values
c = colorbar;
colormap(cmap);
c.Position = POS_CB;

c.Label.String = '\nu';  % Greek letter nu
c.Label.FontSize = FS_LABEL;
c.Label.FontWeight = 'bold';
c.Label.Rotation = 0; 
c.Label.Units = 'normalized';
c.Label.Position = [-0.6, 0.5, 0];  % Position left of colorbar
c.Label.VerticalAlignment = 'middle';
caxis([min(nu_levels) max(nu_levels)]);  % Set colorbar limits
set(gca, 'ColorScale', 'log');  % Use log scaling for color mapping

%% FIGURE 3: ISOLINES AT FIXED GAMMA VALUES
% Shows curves of optimal (p,q) for constant gamma values
gamma_levels = logspace(-1, 1, 10);  % 10 log-spaced gamma levels

% Find indices in gamma_vals closest to the desired levels
gamma_idx = zeros(size(gamma_levels));
for k = 1:numel(gamma_levels)
    [~, gamma_idx(k)] = min(abs(gamma_vals - gamma_levels(k)));
end
gamma_idx = unique(gamma_idx);  % Remove duplicates

figure('Name','Isolines: Fixed Gamma','Position',[200 200 800 700], 'Color', 'w');
set(gcf, 'PaperPositionMode', 'auto');
hold on; box on; grid on;

% Plot boundary envelope as reference
plot(P_boundary, Q_boundary, 'k-', 'LineWidth', LW_BOLD);

% Color map for different gamma levels
num_curves = numel(gamma_idx);
cmap = parula(num_curves);

% Plot each constant-gamma curve
for k = 1:num_curves
    i = gamma_idx(k);
    plot(p_opt(i,:), q_opt(i,:), '-', 'LineWidth', LW_BOLD, 'Color', cmap(k,:));
end

% Labels and formatting
xlabel('p','FontSize',FS_LABEL, 'FontWeight', 'bold');
ylabel('q','FontSize',FS_LABEL, 'FontWeight', 'bold');
set(gca,'FontSize',FS_AXIS, 'LineWidth', 2, 'FontWeight', 'bold');
xlim([p_min - px, p_max + px]);
ylim([q_min - qx, q_max + qx]);

% Apply fixed layout
set(gca, 'Position', POS_AX); 

% Add colorbar for gamma values
c = colorbar;
colormap(cmap);
c.Position = POS_CB;

c.Label.String = '\gamma';  % Greek letter gamma
c.Label.FontSize = FS_LABEL;
c.Label.FontWeight = 'bold';
c.Label.Rotation = 0; 
c.Label.Units = 'normalized';
c.Label.Position = [-0.6, 0.5, 0];  % Position left of colorbar
c.Label.VerticalAlignment = 'middle';
caxis([min(gamma_levels) max(gamma_levels)]);  % Set colorbar limits
set(gca, 'ColorScale', 'log');  % Use log scaling for color mapping

%% FIGURE 4: CONTOUR PLOT OF CONTRACTION RATE IN (p,q) SPACE
% Shows the contraction rate rho as a function of p and q

% Create interpolation function for scattered (p,q,rho) data
mask_scatt = isfinite(p_opt) & isfinite(q_opt) & isfinite(rho_opt);
F_rho = scatteredInterpolant(p_opt(mask_scatt), q_opt(mask_scatt), ...
                             rho_opt(mask_scatt), 'natural', 'none');

% Create fine grid for interpolation
Np = 200; Nq = 200;
[Pgrid, Qgrid] = meshgrid(linspace(p_min, p_max, Np), ...
                          linspace(q_min, q_max, Nq));

% Interpolate rho values on the grid
Rgrid = F_rho(Pgrid, Qgrid);

% Mask values outside the boundary envelope
[in_poly, on_poly] = inpolygon(Pgrid, Qgrid, P_boundary, Q_boundary);
Rgrid(~(in_poly | on_poly)) = NaN;  % Set outside values to NaN

% Create contour plot
figure('Name','Contraction p-q','Position',[250 250 900 750], 'Color', 'w');
set(gcf, 'PaperPositionMode', 'auto');
hold on; box on; grid on;

% Plot filled contours (20 levels)
contourf(Pgrid, Qgrid, Rgrid, 20, 'LineStyle','none');

% Overlay boundary envelope
plot(P_boundary, Q_boundary, 'k-', 'LineWidth', LW_BOLD);

% Labels and formatting
xlabel('p','FontSize',FS_LABEL, 'FontWeight', 'bold');
ylabel('q','FontSize',FS_LABEL, 'FontWeight', 'bold');
set(gca,'FontSize',FS_AXIS, 'LineWidth', 2, 'FontWeight', 'bold');
xlim([p_min - px, p_max + px]);
ylim([q_min - qx, q_max + qx]);

% Apply fixed layout
set(gca, 'Position', POS_AX);

% Add colorbar for contraction rate rho
c = colorbar;
colormap(parula);
c.Position = POS_CB;

c.Label.String = '$\hat{\rho}$';  % Greek letter rho (contraction rate)
c.Label.Interpreter = 'latex';
c.Label.FontSize = FS_LABEL;
c.Label.FontWeight = 'bold';
c.Label.Rotation = 0; 
c.Label.Units = 'normalized';
c.Label.Position = [-0.6, 0.5, 0];  % Position left of colorbar
c.Label.VerticalAlignment = 'middle';

%% FIGURE 5: CONTRACTION RATE IN ORIGINAL (gamma,nu) PARAMETER SPACE
% Shows how contraction rate varies with the original parameters
figure('Name','Contraction gamma-nu','Position',[300 300 900 750], 'Color', 'w');
set(gcf, 'PaperPositionMode', 'auto');
hold on; box on; grid on;

% Create meshgrid for contour plot
[GammaGrid, NuGrid] = meshgrid(gamma_vals, nu_vals);

% Plot filled contours (30 levels)
contourf(GammaGrid, NuGrid, rho_opt', 30, 'LineStyle', 'none');

% Labels and formatting
xlabel('\gamma', 'FontSize', FS_LABEL, 'FontWeight', 'bold');
ylabel('\nu',    'FontSize', FS_LABEL, 'FontWeight', 'bold');
set(gca, 'FontSize', FS_AXIS, 'LineWidth', 2, 'FontWeight', 'bold');
set(gca, 'XScale', 'log', 'YScale', 'log');  % Log-log scale

% Apply fixed layout
set(gca, 'Position', POS_AX);

% Add colorbar for contraction rate rho
c = colorbar;
colormap(parula);
c.Position = POS_CB;

c.Label.String = '$\hat{\rho}$';  % Greek letter rho (contraction rate)
c.Label.Interpreter = 'latex';
c.Label.FontSize = FS_LABEL;
c.Label.FontWeight = 'bold';
c.Label.Rotation = 0; 
c.Label.Units = 'normalized';
c.Label.Position = [-0.6, 0.5, 0];  % Position left of colorbar
c.Label.VerticalAlignment = 'middle';