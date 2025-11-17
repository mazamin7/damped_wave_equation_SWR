clear all; close all; clc;
addpath("utils\")

fprintf('Scanning nu values for final relative Linf error (gamma = 0).\n');

%% Domain and numerics
N = 2; a = 0.3; M = 0.1; Lx = N*a + M;
T = 1; c = 1.0;

dh = 0.002;
dt = 0.002;

%% Physics (gamma = 0 fixed)
gamma = 0;

%% Single eigenmode parameters
n0    = 1;
A0    = 1.0;
v0amp = 0.0;

k0     = n0*pi/Lx;
omega0 = c*k0;

%% Nu range
nu_vals = 0:0.01:1;
Nnu = numel(nu_vals);

rel_Linf_final = zeros(Nnu,1);

%% Loop over nu values
for inu = 1:Nnu
    nu = nu_vals(inu);
    fprintf('nu = %.2f\n', nu);

    %% FDTD simulation
    u_fdtd = run_fdtd_1D(@(x) A0*sin(k0*x), ...
                         @(x) v0amp*sin(k0*x), ...
                         Lx, T, c, dh, dt, gamma, nu);

    Nx = size(u_fdtd,1);
    Nt = size(u_fdtd,2);

    x = linspace(0,Lx,Nx).';
    t = linspace(0,T,Nt);

    %% Analytical solution
    u_an = analytic_solution_single_mode(x, t, c, gamma, nu, k0, A0, v0amp);

    %% Final-time error
    diff  = u_fdtd(:,end) - u_an(:,end);
    abs_Linf_final = max(abs(diff));
    abs_Linf_ref   = max(abs(u_an(:,end)));
    rel_Linf_final(inu) = abs_Linf_final / max(abs_Linf_ref, 1e-15);
end

%% Plot results
figure;
semilogy(nu_vals, rel_Linf_final, 'LineWidth', 2);
xlabel('\nu'); ylabel('Final relative L_\infty error');
title('Final relative Linf error vs viscoelastic coefficient \nu (with \gamma = 0)');
grid on;
