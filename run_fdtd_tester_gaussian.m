% FDTD_GAUSSIAN_PROPAGATION.m
% Simulates a Gaussian pulse centered at 2/3 Lx.
% Time T is chosen short enough so waves do not reach the boundaries (Infinite Domain approx).
clear all; close all; clc;
addpath("utils\") 

fprintf('FDTD Gaussian Propagation (Infinite Domain Approximation).\n');

%% 1. Domain and Numerics Configuration
P = get_sim_params();

N  = P.N; % ignore
a  = P.a; % ignore
M  = P.M; % ignore
b  = P.b; % ignore
Lx = P.Lx; % ignore
Lx = 10;
% T  = P.T;

% c  = P.c;
c = 1;
gamma = P.gamma;
nu = P.nu;

dh = P.dh;
dt = P.dt;
J  = P.J; % ignore


% "Infinite Domain" Condition:
% The pulse is at 2/3 Lx. The distance to the nearest boundary (Right) is 1/3 Lx.
% We need c*T < 1/3 Lx to avoid reflections.
dist_to_bound = (1/3) * Lx;
T = 0.9 * (dist_to_bound / c); % Set T slightly less than travel time to boundary

%% 2. Initial Conditions (Gaussian)
A0    = 1.0;
x_c   = (2/3) * Lx;       % Center location
sigma = Lx / 40;          % Pulse width (narrow enough to be distinct)

% Gaussian displacement, zero velocity
u0 = @(x) A0 * exp( - (x - x_c).^2 ./ (2*sigma^2) );
v0 = @(x) zeros(size(x)); 

% Visual check of Initial Conditions
x_ic = linspace(0, Lx, 1000);
figure('Name','Initial conditions');
plot(x_ic, u0(x_ic), 'LineWidth', 2); 
grid on; xlabel('x'); ylabel('u_0');
title('Initial Condition: Gaussian at x = 2/3 L_x');
xlim([0 Lx]);

%% 3. FDTD Solution
fprintf('Running FDTD...\n');
% Utilizing the existing run_fdtd utility
u_fdtd = run_fdtd(u0, v0, Lx, T, c, dh, dt, gamma, nu);

%% 4. Visualization
% Generate grid vectors for plotting
Nx_sol = size(u_fdtd, 1);
Nt_sol = size(u_fdtd, 2);
x = linspace(0, Lx, Nx_sol);
t = linspace(0, T, Nt_sol);

figure('Name', 'Space-Time Evolution');
mesh(x, t, u_fdtd'); 
xlabel('x'); ylabel('t'); zlabel('u');
title('FDTD Solution (No Boundary Reflections)');
view(0, 90); % Top-down view to see the spread
colorbar;
axis tight;

fprintf('Simulation complete.\n');