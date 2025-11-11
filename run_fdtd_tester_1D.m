% FDTD_vs_ANALYTIC_1MODE_1D — single mode, Dirichlet, damping consistent with u_tt+gamma u_t = c^2 u_xx + nu u_txx
clear all; close all; clc;
addpath("utils\")

fprintf('FDTD vs analytical, single mode, Dirichlet BCs.\n');

%% Domain and numerics
N = 2; a = 0.3; M = 0.1; Lx = N*a + M;
T = 1.4; c = 1.0;

% target steps
dh_t = 1e-3;
dt_t = 1e-3;

% snap to grid-exact values to satisfy util asserts
Nx = max(3, round(Lx/dh_t) + 1);
Nt = max(3, round(T /dt_t) + 1);
dh = Lx/(Nx-1);
dt = T /(Nt-1);
CFL = c*dt/dh;
fprintf('Snapped: dh_t=%.3e -> dh=%.3e, dt_t=%.3e -> dt=%.3e, CFL=%.3g\n', dh_t, dh, dt_t, dt, CFL);

%% Physics (PDE: u_tt + gamma u_t = c^2 u_xx + nu u_txx)
gamma = 0.1;
nu    = 0.1;

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

%% Analytical solution consistent with u_tt + gamma u_t = c^2 u_xx + nu u_txx (no forcing)
% Modal ODE: q'' + (gamma + nu*k0^2) q' + omega0^2 q = 0
A0x  = A0;  V0x = v0amp;
ge   = 0.5*(gamma + nu*k0^2);    % effective half-damping  (NOTE: k0^2, not omega0^2)
disc = ge^2 - omega0^2;

if disc < -1e-14
    omegad = sqrt(omega0^2 - ge^2);
    Ct = cos(omegad*t); St = sin(omegad*t);
    q  = exp(-ge*t) .* ( A0x*Ct + ((V0x + ge*A0x)/omegad) * St );
elseif abs(disc) <= 1e-14
    q  = exp(-ge*t) .* ( A0x + (V0x + ge*A0x)*t );
else
    s  = sqrt(disc);
    r1 = -ge + s;  r2 = -ge - s;
    C1 = (V0x - r2*A0x)/(r1 - r2);
    C2 = (r1*A0x - V0x)/(r1 - r2);
    q  = C1*exp(r1*t) + C2*exp(r2*t);
end

u_an = sin(k0*x) * q;

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
