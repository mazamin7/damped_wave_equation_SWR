clear; close all; clc;
addpath("utils\")

%% Parameters (edit here)
c     = 1.0;        % wave speed
gamma = 0;        % telegraph damping
nu    = 0;        % viscoelastic damping

dx = 0.002;
dt = 0.002;

fprintf('Dispersion analysis for c=%.3f, gamma=%.3f, nu=%.3f, dx=%.4f, dt=%.4f\n', ...
        c, gamma, nu, dx, dt);

%% Non-dimensional parameters
r2   = (c*dt/dx)^2;           % (c dt / dx)^2
sigma = nu*dt/(2*dx^2);       % nu dt / (2 dx^2)
G     = gamma*dt/2;           % gamma dt / 2

%% Wavenumber range
% Avoid k=0 and Nyquist to skip degenerate cases / branch issues
Nk = 400;
kmin = 1e-6;                  % small positive
kmax = pi/dx - 1e-6;          % just below Nyquist
k    = linspace(kmin, kmax, Nk);
kdx  = k*dx;                  % nondimensional wave number

%% Storage for numerical frequency
omega_num = zeros(1, Nk);

for j = 1:Nk
    kappa = kdx(j);                  % nondimensional wavenumber

    % Discrete Laplacian symbol: s = (e^{ikdx} -2 + e^{-ikdx})/dx^2
    % but we factor out 1/dx^2, so here s is dimensionless:
    s = 2*cos(kappa) - 2;            % = -4 sin^2(kappa/2), no 1/dx^2

    % From the scheme:
    % (u^{n+1}-2u^n+u^{n-1})/dt^2 + (gamma/(2dt))(u^{n+1}-u^{n-1}) =
    %   c^2/dx^2 * s * u^n + (nu/(2dt dx^2)) * s * (u^{n+1}-u^{n-1})
    %
    % Assume modal ansatz u^n ~ lambda^n:
    %
    % -> [lambda - 2 + lambda^{-1}] + G (lambda - lambda^{-1})
    %    = r2 * s + sigma * s (lambda - lambda^{-1})
    %
    % Multiply by lambda to get a quadratic:
    %
    %   A*lambda^2 + B*lambda + C = 0
    %
    A = 1 + G - sigma * s;
    B = -2 - r2 * s;
    C = 1 - G + sigma * s;

    roots_lambda = roots([A, B, C]);

    % Select the physically relevant root: |lambda| <= 1 (decaying / stable mode)
    [~, idx_min] = min(abs(roots_lambda));  % could also pick the one with |lambda|<=1
    lambda = roots_lambda(idx_min);

    % Numerical frequency from lambda = exp(-i omega dt)
    % -> omega = -i/dt * log(lambda)
    omega_num(j) = -1i/dt * log(lambda);
end

%% Physical dispersion relation
% PDE: u_tt + gamma u_t = c^2 u_xx + nu u_txx
% Try u(x,t) = exp(i(kx - omega t)):
%   -omega^2 + i gamma omega = -c^2 k^2 - i nu k^2 omega
% => omega^2 + i (gamma + nu k^2) omega - c^2 k^2 = 0
%
% Solve quadratic in omega:
%   omega^2 + i alpha omega - c^2 k^2 = 0,  alpha = gamma + nu k^2
%
% omega_phys = [ -i alpha +/- sqrt(-alpha^2 + 4 c^2 k^2) ] / 2
alpha = gamma + nu*k.^2;
disc  = -alpha.^2 + 4*c^2.*k.^2;
sqrt_disc = sqrt(disc);

omega_phys_1 = (-1i*alpha + sqrt_disc)/2;
omega_phys_2 = (-1i*alpha - sqrt_disc)/2;

% Choose one branch (the one with positive real frequency magnitude)
% This is somewhat arbitrary; adjust if you want the other mode.
choose1 = abs(real(omega_phys_1)) >= abs(real(omega_phys_2));
omega_phys = omega_phys_1;
omega_phys(~choose1) = omega_phys_2(~choose1);

%% Phase velocity and damping
% Phase velocity: Re(omega)/k
phase_num  = real(omega_num)./k;
phase_phys = real(omega_phys)./k;

% Damping rate: -Im(omega) (should be >= 0 if solution decays)
damp_num  = -imag(omega_num);
damp_phys = -imag(omega_phys);

%% Plot: phase velocity vs k
figure;
subplot(1,2,1);
plot(k, phase_phys, 'k-', 'LineWidth', 2); hold on;
plot(k, phase_num,  'r--', 'LineWidth', 2);
xlabel('wavenumber k');
ylabel('phase velocity v_p = \Re(\omega)/k');
title('Phase velocity: physical vs numerical');
legend('physical','numerical','Location','best');
grid on;

%% Plot: damping rate vs k
subplot(1,2,2);
plot(k, damp_phys, 'k-', 'LineWidth', 2); hold on;
plot(k, damp_num,  'r--', 'LineWidth', 2);
xlabel('wavenumber k');
ylabel('damping rate -\Im(\omega)');
title('Damping: physical vs numerical');
legend('physical','numerical','Location','best');
grid on;

sgtitle(sprintf('Dispersion and damping: c=%.2f, \\gamma=%.2f, \\nu=%.2f, dt=%.4f, dx=%.4f', ...
    c, gamma, nu, dt, dx));
