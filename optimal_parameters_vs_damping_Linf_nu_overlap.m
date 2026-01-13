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
% M  = P.M;      % Another system parameter
M = 1e-4;
% M = 0;
T  = P.T;      % Total simulation time
c  = P.c;      % Wave speed or similar constant
dt = P.dt;     % Time step size
J  = P.J;      % Parameter for objective function

%% OPTIMIZATION CONFIGURATION
tol = 1e-9;        % Tolerance for optimization convergence
padding = 0.05;    % Padding for plot axis limits (5% extra space)

%% FINE GRID DEFINITION FOR PARAMETER SPACE EXPLORATION
% We set gamma to 0 and vary nu widely
Ng = 1;         % Gamma is fixed
Nn = 100;       % High resolution for nu to make the color gradient smooth
% Parameter ranges:
gamma_vals = 0;                     % Fixed Gamma = 0
nu_vals    = logspace(-5, -3, Nn);   % Vary Nu from 1e-5 to 1e-2

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
fprintf('=== Computing optimal (p,q) for Gamma=0, varying Nu ===\n');
% Loop over nu values (outer loop)
for j = 1:Nn
    nu = nu_vals(j);
    % Loop over gamma values (inner loop - size 1)
    for i = 1:Ng
        gamma = gamma_vals(i);
        
        %% SMART INITIALIZATION STRATEGY
        if j > 1
            % Use previous nu value (same gamma) as initial guess
            x0 = [p_opt(i,j-1), q_opt(i,j-1)];
        else
            % First point: use fallback
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
        
        %% ERROR HANDLING
        if ~isfinite(fval) || any(~isfinite(x_opt))
            x_opt = x_fallback;
            fval  = objfun(x_opt);
        end
        
        %% STORE RESULTS
        p_opt(i,j)   = x_opt(1);
        q_opt(i,j)   = x_opt(2);
        rho_opt(i,j) = fval;
    end
    
    % Simple progress indicator
    if mod(j, 10) == 0
        fprintf('  Computed %d/%d (nu = %.2e)\n', j, Nn, nu);
    end
end

%% VISUALIZATION SETTINGS (POSTER STYLE)
% Define consistent styling for all figures
FS_AXIS  = 26;      % Font size for axis tick labels
FS_LABEL = 28;      % Font size for axis labels
LW_BOLD  = 4.0;     % Line width for bold lines (boundaries)
MS_DOTS  = 20;      % Marker size for scatter plot dots

% Fixed layout positions
POS_AX  = [0.15, 0.15, 0.64, 0.78];  % Main axis position
POS_CB  = [0.84, 0.15, 0.04, 0.78];  % Colorbar position

%% --- ANALYTICAL APPROXIMATION CALCULATION ---
% Using the Grand Unified Formulas derived from Minimax theory:

% 1. Compute Simulation Frequencies
w_max = pi / dt;
w_min = pi / T;

% 2. Compute q Approximation (Geometric Mean Damping)
% Formula: q ~ 0.5 * nu * (w_min * w_max) / c^3
q_approx_curve = 0.5 * (nu_vals / c^3) * (w_min * w_max);

% 3. Compute p Approximation (Minimax Phase Correction)
% Formula: p ~ 1/c - 0.25 * nu^2 * w_max^2 / c^5
% p_approx_curve = (1/c) - (nu_vals.^2 / c^5) * ( 5/16 * (w_max^2 + w_min^2) - 1/4 * w_max * w_min + M^2/(6*c^2) * (w_max^2 + w_min^2) * (w_max - w_min)^2 );
p_approx_curve = (1/c) - (nu_vals.^2 / c^5) * ( 4/16 * (w_max^2 + w_min^2) - 1/4 * w_max * w_min + M^2/(6*c^2) * (w_max^2 + w_min^2) * (w_max - w_min)^2 );

% 4. Compute Rho Approximation (Contraction Rate)
% Formula: rho ~ 0.25 * nu * (w_max - w_min) / c^2
rho_approx_curve = 0.25 * (nu_vals / c^2) * (w_max - w_min);


% Determine axis limits dynamically based on NUMERICAL data range
% (We keep the plot focused on the real data, even if approx diverges)
p_data = p_opt(1,:);
q_data = q_opt(1,:);
p_min = min(p_data); p_max = max(p_data);
q_min = min(q_data); q_max = max(q_data);

% Add minimal padding
px = padding * (p_max - p_min); if px==0, px=0.1; end
qx = padding * (q_max - q_min); if qx==0, qx=0.1; end

%% FIGURE 2: CONTINUOUSLY COLORED TRAJECTORY (Varying Nu)
% Plots (p,q) curve where color represents the value of nu
figure('Name','Optimal Parameters vs Nu','Position',[150 150 800 700], 'Color', 'w');
set(gcf, 'PaperPositionMode', 'auto'); 
hold on; box on; grid on;

% Extract data for the single gamma=0 case
x = p_data;
y = q_data;
z = zeros(size(x)); % Z-coordinate for 2D plot
c_data = log10(nu_vals); % Use log10(nu) for color so gradients are visible across orders of magnitude

% --- PLOT ANALYTICAL CURVE FIRST (so it is behind the markers if they overlap) ---
plot(p_approx_curve, q_approx_curve, 'r--', 'LineWidth', 3, 'DisplayName', 'Analytical Approx');

% --- PLOT NUMERICAL SURFACE ---
% Use surface trick to draw a line with changing color
% "FaceColor",'no' makes it a hollow mesh (just edges)
% "EdgeColor",'interp' interpolates color between vertices
% We use a dummy plot for the legend first
h_dummy = plot(nan, nan, 'b-', 'LineWidth', LW_BOLD, 'DisplayName', 'Numerical Opt.');
surface([x;x],[y;y],[z;z],[c_data;c_data], ...
        'FaceColor','no', ...
        'EdgeColor','interp', ...
        'LineWidth', LW_BOLD, ...
        'HandleVisibility', 'off'); % Hide surface from legend

% Mark the start (low nu) and end (high nu)
plot(x(1), y(1), 'ko', 'MarkerSize', 10, 'MarkerFaceColor', 'k', 'HandleVisibility', 'off'); 
text(x(1), y(1), '  Low \nu', 'FontSize', 14, 'FontWeight', 'bold');
plot(x(end), y(end), 'ks', 'MarkerSize', 10, 'MarkerFaceColor', 'k', 'HandleVisibility', 'off'); 
text(x(end), y(end), '  High \nu', 'FontSize', 14, 'FontWeight', 'bold');

% Labels and formatting
xlabel('p','FontSize',FS_LABEL, 'FontWeight', 'bold');
ylabel('q','FontSize',FS_LABEL, 'FontWeight', 'bold');
set(gca,'FontSize',FS_AXIS, 'LineWidth', 2, 'FontWeight', 'bold');
xlim([p_min - px, p_max + px]);
ylim([q_min - qx, q_max + qx]);
title(['Optimal Parameters (\gamma = ' num2str(gamma_vals) ')'], 'FontSize', 20);
legend('show', 'Location', 'best', 'FontSize', 20); % Show legend

% Apply fixed layout
set(gca, 'Position', POS_AX); 

% Add colorbar for nu values
c_bar = colorbar;
colormap(parula); % Or 'jet', 'turbo', etc.
c_bar.Position = POS_CB;

% Configure colorbar to show 10^x labels
c_bar.Label.String = 'log_{10}(\nu)'; 
c_bar.Label.FontSize = FS_LABEL;
c_bar.Label.FontWeight = 'bold';
c_bar.Label.Rotation = 0; 
c_bar.Label.Units = 'normalized';
c_bar.Label.Position = [-0.6, 0.5, 0];
c_bar.Label.VerticalAlignment = 'middle';

% Ensure color axis matches data
caxis([min(c_data) max(c_data)]); 

%% FIGURE 4: CONTRACTION RATE VS NU
% Simple 2D plot showing how performance improves/degrades with nu
figure('Name','Contraction vs Nu','Position',[250 250 800 600], 'Color', 'w');
set(gcf, 'PaperPositionMode', 'auto');
hold on; box on; grid on;

% 1. Plot Numerical Results (Blue Line)
semilogx(nu_vals, rho_opt(1,:), 'b-', 'LineWidth', LW_BOLD, 'DisplayName', 'Numerical Opt.');

% 2. Plot Analytical Prediction (Red Dashed Line)
semilogx(nu_vals, rho_approx_curve, 'r--', 'LineWidth', 3, 'DisplayName', 'Analytical Approx.');

% Labels and Formatting
xlabel('\nu (Viscoelastic Damping)', 'FontSize', FS_LABEL, 'FontWeight', 'bold');
ylabel('$\hat{\rho}$', 'Interpreter', 'latex', 'FontSize', FS_LABEL, 'FontWeight', 'bold');
set(gca,'FontSize',FS_AXIS, 'LineWidth', 2, 'FontWeight', 'bold');
title('Contraction Rate vs \nu', 'FontSize', 20);
legend('show', 'Location', 'best', 'FontSize', 20);