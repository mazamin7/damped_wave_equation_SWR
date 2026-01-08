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

%% OPTIMIZATION CONFIGURATION
tol = 1e-9;        % Tolerance for optimization convergence
padding = 0.05;    % Padding for plot axis limits (5% extra space)

%% FINE GRID DEFINITION FOR PARAMETER SPACE EXPLORATION
% We set nu to 0 and vary gamma widely
Nn = 1;         % Nu is fixed
Ng = 100;       % High resolution for gamma
% Parameter ranges:
nu_vals = 0;                        % Fixed Nu = 0
gamma_vals = logspace(-6, 0, Ng);   % Vary Gamma from 0.01 to 10

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

% Fallback initial guess
x_fallback = [1/c, 0];  % Reasonable default: [p0, q0]

%% MAIN OPTIMIZATION LOOP OVER PARAMETER GRID
fprintf('=== Computing optimal (p,q) for Nu=0, varying Gamma ===\n');
% Loop over nu (fixed size 1)
for j = 1:Nn
    nu = nu_vals(j);
    % Loop over gamma (varying)
    for i = 1:Ng
        gamma = gamma_vals(i);
        
        %% SMART INITIALIZATION STRATEGY
        if i > 1
            % Use previous gamma value as initial guess
            x0 = [p_opt(i-1,j), q_opt(i-1,j)];
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
    
    if mod(i, 10) == 0
        fprintf('  Computed %d/%d (gamma = %.2e)\n', i, Ng, gamma);
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

%% --- ANALYTICAL APPROXIMATION CALCULATION (TELEGRAPHER LIMIT) ---
% 1. Compute Simulation Frequencies
w_max = pi / dt;
w_min = pi / T;

% 2. Compute q Approximation (Constant Damping Match)
% Formula: q ~ gamma / (2 * c)
q_approx_curve = gamma_vals / (2 * c);

% 3. Compute p Approximation (Low Frequency Singular Correction)
% Formula: p ~ 1/c + gamma^2 / (16 * c * w_min^2)
% Note: The error is dominated by w_min in this regime
p_approx_curve = (1/c) + (gamma_vals.^2) / (16 * c * w_min^2);

% 4. Compute Rho Approximation (Contraction Rate)
% Formula: rho ~ 1 / (32 * w_min^2) (Linear improvement with damping)
% Note: This is valid for small gamma. For large gamma, it saturates.
rho_approx_curve = gamma_vals.^2 / (32 * w_min^2);
% rho_approx_curve = gamma_vals.^2 .* sqrt(9 * gamma_vals.^2 + 4 * w_min^2) / (64 * w_min^3);


% Determine axis limits dynamically based on NUMERICAL data range
p_data = p_opt(:,1)'; % Transpose to row vector
q_data = q_opt(:,1)';
p_min = min(p_data); p_max = max(p_data);
q_min = min(q_data); q_max = max(q_data);

% Add minimal padding
px = padding * (p_max - p_min); if px==0, px=0.1; end
qx = padding * (q_max - q_min); if qx==0, qx=0.1; end

%% FIGURE 2: CONTINUOUSLY COLORED TRAJECTORY (Varying Gamma)
% Plots (p,q) curve where color represents the value of gamma
figure('Name','Optimal Parameters vs Gamma','Position',[150 150 800 700], 'Color', 'w');
set(gcf, 'PaperPositionMode', 'auto'); 
hold on; box on; grid on;

% Extract data for the single nu=0 case
x = p_data;
y = q_data;
z = zeros(size(x)); % Z-coordinate for 2D plot
c_data = log10(gamma_vals); % Use log10(gamma) for color

% --- PLOT ANALYTICAL CURVE FIRST ---
plot(p_approx_curve, q_approx_curve, 'r--', 'LineWidth', 3, 'DisplayName', 'Analytical Approx');

% --- PLOT NUMERICAL SURFACE ---
h_dummy = plot(nan, nan, 'b-', 'LineWidth', LW_BOLD, 'DisplayName', 'Numerical Opt.');
surface([x;x],[y;y],[z;z],[c_data;c_data], ...
        'FaceColor','no', ...
        'EdgeColor','interp', ...
        'LineWidth', LW_BOLD, ...
        'HandleVisibility', 'off'); 

% Mark start and end
plot(x(1), y(1), 'ko', 'MarkerSize', 10, 'MarkerFaceColor', 'k', 'HandleVisibility', 'off'); 
text(x(1), y(1), '  Low \gamma', 'FontSize', 14, 'FontWeight', 'bold');
plot(x(end), y(end), 'ks', 'MarkerSize', 10, 'MarkerFaceColor', 'k', 'HandleVisibility', 'off'); 
text(x(end), y(end), '  High \gamma', 'FontSize', 14, 'FontWeight', 'bold');

% Labels and formatting
xlabel('p','FontSize',FS_LABEL, 'FontWeight', 'bold');
ylabel('q','FontSize',FS_LABEL, 'FontWeight', 'bold');
set(gca,'FontSize',FS_AXIS, 'LineWidth', 2, 'FontWeight', 'bold');
xlim([p_min - px, p_max + px]);
ylim([q_min - qx, q_max + qx]);
title(['Optimal Parameters (\nu = ' num2str(nu_vals) ')'], 'FontSize', 20);
legend('show', 'Location', 'best', 'FontSize', 20); 

% Apply fixed layout
set(gca, 'Position', POS_AX); 

% Add colorbar
c_bar = colorbar;
colormap(parula); 
c_bar.Position = POS_CB;
c_bar.Label.String = 'log_{10}(\gamma)'; 
c_bar.Label.FontSize = FS_LABEL;
c_bar.Label.FontWeight = 'bold';
c_bar.Label.Rotation = 0; 
c_bar.Label.Units = 'normalized';
c_bar.Label.Position = [-0.6, 0.5, 0];
c_bar.Label.VerticalAlignment = 'middle';
caxis([min(c_data) max(c_data)]); 

%% FIGURE 4: CONTRACTION RATE VS GAMMA
figure('Name','Contraction vs Gamma','Position',[250 250 800 600], 'Color', 'w');
set(gcf, 'PaperPositionMode', 'auto');
hold on; box on; grid on;

% 1. Plot Numerical Results (Blue Line)
semilogx(gamma_vals, rho_opt(:,1)', 'b-', 'LineWidth', LW_BOLD, 'DisplayName', 'Numerical Opt.');

% 2. Plot Analytical Prediction (Red Dashed Line)
semilogx(gamma_vals, rho_approx_curve, 'r--', 'LineWidth', 3, 'DisplayName', 'Analytical Approx.');

% Labels and Formatting
xlabel('\gamma (Telegrapher Damping)', 'FontSize', FS_LABEL, 'FontWeight', 'bold');
ylabel('$\hat{\rho}$', 'Interpreter', 'latex', 'FontSize', FS_LABEL, 'FontWeight', 'bold');
set(gca,'FontSize',FS_AXIS, 'LineWidth', 2, 'FontWeight', 'bold');
title('Contraction Rate vs \gamma', 'FontSize', 20);
legend('show', 'Location', 'best', 'FontSize', 20);