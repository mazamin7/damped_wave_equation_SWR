% FDTD_vs_ANALYTIC_1MODE_1D — single mode, Dirichlet, damping consistent with u_tt+gamma u_t = c^2 u_xx + nu u_txx
clear all; close all; clc;
addpath("utils\")

fprintf('FDTD vs analytical, single mode, Dirichlet BCs.\n');

%% Domain and numerics
N = 2; a = 0.3; M = 0.1; Lx = N*a + M;
T = 1.4; c = 1.0;

% target steps
dh = 0.002;
dt = 0.002;

CFL = c*dt/dh;

%% Physics (PDE: u_tt + gamma u_t = c^2 u_xx + nu u_txx)
gamma = 10;
nu    = 0;

%% Single eigenmode
n0    = 1;
A0    = 1.0;
v0amp = 0.0;

k0     = n0*pi/Lx;        % wavenumber
omega0 = c*k0;            % natural frequency

%% Initial conditions
u0 = @(x) A0    * sin(k0*x);
v0 = @(x) v0amp * sin(k0*x);

% Visual check
x_ic = linspace(0,Lx,1000);
figure('Name','Initial conditions');
subplot(2,1,1); plot(x_ic,u0(x_ic),'LineWidth',2); grid on; xlabel('x'); ylabel('u_0');
subplot(2,1,2); plot(x_ic,v0(x_ic),'LineWidth',2); grid on; xlabel('x'); ylabel('v_0');
sgtitle('Initial conditions (single sine mode)');

%% FDTD solution (util enforces Dirichlet and damping terms)
fprintf('Running FDTD...\n');
u_fdtd = run_fdtd_1D(u0, v0, Lx, T, c, dh, dt, gamma, nu);

Nx = size(u_fdtd,1);
Nt = size(u_fdtd,2);
x  = linspace(0,Lx,Nx).';
t  = linspace(0,T,Nt);

u_an = analytic_solution_single_mode(x, t, c, gamma, nu, k0, A0, v0amp);

%% Final-time errors
u_f   = u_fdtd(:,end);
u_a   = u_an(:,end);
diff  = u_f - u_a;

abs_Linf_final = max(abs(diff));
rel_Linf_final = abs_Linf_final / max(max(abs(u_a)), 1e-15);

dx_eff   = Lx/(Nx-1);
abs_L2_final = sqrt(dx_eff * sum(diff.^2));
rel_L2_final = abs_L2_final / max(sqrt(dx_eff * sum(u_a.^2)), 1e-15);

fprintf('Final absolute L2  error: %.3e\n', abs_L2_final);
fprintf('Final relative L2  error: %.3e\n', rel_L2_final);
fprintf('Final absolute Linf error: %.3e\n', abs_Linf_final);
fprintf('Final relative Linf error: %.3e\n', rel_Linf_final);

%% Space–time fields
figure('Name','Space–time fields');
subplot(1,2,1)
mesh(x, t, u_fdtd'); xlabel('x'); ylabel('t'); zlabel('u');
title('FDTD (numerical)'); grid on; view(120,30); colorbar;
subplot(1,2,2)
mesh(x, t, u_an');   xlabel('x'); ylabel('t'); zlabel('u');
title('Analytical'); grid on; view(120,30); colorbar;
sgtitle('Space–time solutions: numerical vs analytical');
