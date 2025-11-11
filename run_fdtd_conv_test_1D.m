%% convergence_test_damped_wave_using_util.m
clearvars; close all; clc;
addpath("utils\")

% Physical parameters
c = 1;  L = 1.0;  T = 1.0;

k = 1;  A0 = 1.0;  V0 = 0.0;

% Old damping -> new PDE: u_tt + gamma u_t = c^2 u_xx + nu u_txx
alpha1 = 0.05 * c;
alpha2 = 0.05 * c;
gamma  = 2*alpha1;
nu     = 2*alpha2*c^2;

% target dx list
dx_target = [0.1, 0.1/2, 0.1/4, 0.1/8, 0.1/16, 0.1/32, 0.1/64];

num_cases = numel(dx_target);
L2_err  = zeros(num_cases,1);
Linf_err = zeros(num_cases,1);
dx_used = zeros(num_cases,1);
dt_used = zeros(num_cases,1);

fprintf('PDE: u_tt + gamma u_t = c^2 u_xx + nu u_txx,  gamma=%.4g, nu=%.4g\n', gamma, nu);

u0_fun = @(x) A0 * sin(k*pi*x/L);
v0_fun = @(x) V0 * sin(k*pi*x/L);

for ic = 1:num_cases
    % --- snap steps so util asserts pass ---
    dx_t   = dx_target(ic);
    Nx     = max(3, round(L/dx_t) + 1);
    dx     = L/(Nx-1);
    dt_t   = dx/c;                      % target CFL=1
    Nt     = max(3, round(T/dt_t) + 1);
    dt     = T/(Nt-1);
    CFL    = c*dt/dx;

    fprintf('Case %d: dx_t=%.3e -> dx=%.3e, dt=%.3e, CFL=%.3g, Nx=%d, Nt=%d\n', ...
            ic, dx_t, dx, dt, CFL, Nx, Nt);

    % --- numerical solution on snapped grid ---
    u_fdtd = run_fdtd_1D(u0_fun, v0_fun, L, T, c, dx, dt, gamma, nu);
    % grid that matches u_fdtd
    x = linspace(0,L,size(u_fdtd,1)).';
    t = linspace(0,T,size(u_fdtd,2));

    % --- analytical single-mode (no forcing) consistent with PDE ---
    lambda_k = (k*pi/L)^2;
    omega0   = c*sqrt(lambda_k);
    ge       = 0.5*(gamma + nu*lambda_k);
    disc     = ge^2 - omega0^2;

    if disc < -1e-14
        omegad = sqrt(omega0^2 - ge^2);
        q = exp(-ge*t) .* ( A0*cos(omegad*t) + ((V0 + ge*A0)/omegad)*sin(omegad*t) );
    elseif abs(disc) <= 1e-14
        q = exp(-ge*t) .* ( A0 + (V0 + ge*A0)*t );
    else
        s  = sqrt(disc);
        r1 = -ge + s;  r2 = -ge - s;
        C1 = (V0 - r2*A0)/(r1 - r2);
        C2 = (r1*A0 - V0)/(r1 - r2);
        q  = C1*exp(r1*t) + C2*exp(r2*t);
    end
    u_an = sin(k*pi*x/L) * q;

    % --- final-time errors on snapped grid ---
    diff = u_fdtd(:,end) - u_an(:,end);
    dx_eff = L/(numel(x)-1);
    L2_err(ic)   = sqrt(dx_eff * sum( diff.^2 ));
    Linf_err(ic) = max(abs(diff));

    dx_used(ic) = dx_eff;
    dt_used(ic) = dt;
end

% rates
p_L2   = polyfit(log(dx_used),   log(max(L2_err,  eps)), 1);
p_Linf = polyfit(log(dx_used),   log(max(Linf_err,eps)), 1);
rate_L2   = p_L2(1);
rate_Linf = p_Linf(1);

fprintf('\nConvergence results at T = %.3f\n', T);
fprintf(' dx        dt          L2_err       rate_L2    Linf_err    rate_Linf\n');
for i=1:num_cases
    if i==1, rL2=NaN; rLinf=NaN;
    else
        rL2   = log(L2_err(i)/L2_err(i-1))/log(dx_used(i)/dx_used(i-1));
        rLinf = log(Linf_err(i)/Linf_err(i-1))/log(dx_used(i)/dx_used(i-1));
    end
    fprintf('%.5g  %.5g  %.3e  %.3f     %.3e  %.3f\n', dx_used(i), dt_used(i), L2_err(i), rL2, Linf_err(i), rLinf);
end
fprintf('\nEstimated global rates (linear fit): rate_L2 = %.3f, rate_Linf = %.3f\n', rate_L2, rate_Linf);

% plot with h^2 ref
figure('Name','Final-time error vs dx');
loglog(dx_used, L2_err,   'o-', 'LineWidth',1.8, 'DisplayName','L2'); hold on;
loglog(dx_used, Linf_err, 's--','LineWidth',1.8, 'DisplayName','L_\infty');
[dx_min, idx_min] = min(dx_used);
ref_h2 = L2_err(idx_min) * (dx_used/dx_min).^2;
loglog(dx_used, ref_h2, ':', 'LineWidth',1.5, 'DisplayName','h^2 ref');
set(gca,'XDir','reverse'); grid on; legend('Location','southwest');
xlabel('h = \Delta x'); ylabel('Error at T');
title(sprintf('Convergence: rate L2=%.2f, rate L_\\infty=%.2f', rate_L2, rate_Linf));

% sanity check vs old formula
A_old = modalA_old(t, A0, V0, alpha1, alpha2, c, k, L);
fprintf('Max |A_old - A_new|: %.3e\n', max(abs(A_old - q)));

function At = modalA_old(t, A0, V0, alpha1, alpha2, c, k, L)
    lambda_k = (k*pi/L)^2;  beta = alpha1 + alpha2*c^2*lambda_k;  omega2 = c^2*lambda_k;
    t = t(:)'; 
    if beta^2 < omega2
        od = sqrt(omega2 - beta^2);
        At = exp(-beta*t).*( A0*cos(od*t) + ((V0 + beta*A0)/od).*sin(od*t) );
    elseif abs(beta^2-omega2) < 1e-12
        At = exp(-beta*t).*( A0 + (V0 + beta*A0).*t );
    else
        s = sqrt(beta^2 - omega2); r1=-beta+s; r2=-beta-s;
        C1=(V0 - r2*A0)/(r1-r2); C2=A0-C1;
        At = C1*exp(r1*t) + C2*exp(r2*t);
    end
end
