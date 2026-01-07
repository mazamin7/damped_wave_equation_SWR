%% SCRIPT 4: TELEGRAPHER COEFFICIENT FITTING (ROBUST)
% Fits Kp, Kq and the non-linear integers A, B for the refined rho formula
clear; close all; clc; addpath("utils\");

c = 1.0; nu = 0; N = 2; a = 0.5; M = 0; J = 1000;

% Tight tolerances are crucial as p-correction scales with gamma^2
options = optimset('Display','off', 'TolX',1e-14, 'TolFun',1e-14);

fprintf('--- TELEGRAPHER REFINED FITTING ---\n');
NumPoints = 100;

% 1. Logarithmic Random Sampling
% gamma: 10^-3 to 10^-0.5 (~0.0001 to 0.01)
% We avoid extremely small gamma (<1e-4) to ensure p-correction signal > noise
gamma_rand = 10.^( -6 + 2 * rand(NumPoints, 1) );

% w_min: 0.1 to 10 (Simulation time scales)
w_min_rand = 10.^( -1 + 1.5 * rand(NumPoints, 1) );

% w_max: Fixed high enough to be asymptotic
w_max = 100; 
dt = pi/w_max;

Y_q = zeros(NumPoints,1); 
Y_p = zeros(NumPoints,1); 
Y_rho = zeros(NumPoints,1);

fprintf('Generating data...\n');
for i=1:NumPoints
    g = gamma_rand(i); 
    w = w_min_rand(i);
    T = pi/w;
    
    % Smart Initial Guess
    % p ~ 1/c + g^2/(16 c w^2)
    p_guess = 1/c + g^2/(16*c*w^2);
    q_guess = g/(2*c);
    
    obj = @(x) obj_Linf(N,T,dt,J,c,g,nu,a,M,x(1),x(2));
    [res, fval] = fminsearch(obj, [p_guess, q_guess], options);
    
    Y_p(i) = res(1) - 1/c; % Positive correction
    Y_q(i) = res(2);
    Y_rho(i) = fval;
    
    if mod(i, 20) == 0, fprintf('Processed %d/%d...\n', i, NumPoints); end
end

%% STEP 1: LINEAR COEFFICIENT FITTING (Robust)

% 1. Fit q (Linear in gamma -> Strong signal)
X_q = gamma_rand / c;
K_q = X_q \ Y_q;

% 2. Fit p (Quadratic in gamma -> Needs filtering)
X_p = (gamma_rand.^2) ./ (c * w_min_rand.^2);

% Filter Outliers (Noise Floor Rejection)
% Rejects points where optimizer didn't move enough or hit numerical noise
valid_p_idx = (Y_p > 1e-13) & (Y_p < 10.0);

if sum(valid_p_idx) < 5
    warning('Too few valid points for p-fit.');
    K_p = 0;
else
    K_p = X_p(valid_p_idx) \ Y_p(valid_p_idx);
end

fprintf('\nLINEAR RESULTS:\n');
fprintf('  K_q (Expect 0.50):   %.4f\n', K_q);
fprintf('  K_p (Expect 0.0625): %.4f  (Fitted on %d valid points)\n', K_p, sum(valid_p_idx));

%% STEP 2: NON-LINEAR RHO FITTING (Integers A, B)
fprintf('\nNON-LINEAR RESULTS (Refined Rho):\n');

% Model: rho = gamma^2 * sqrt(A*g^2 + B*w^2) / (64 * w^3)
C_fixed = 64;
func_rho = @(p, x) (x(:,1).^2) .* sqrt(p(1)*x(:,1).^2 + p(2)*x(:,2).^2) ./ (C_fixed * x(:,2).^3);

X_data = [gamma_rand, w_min_rand];

% We use all data for rho fitting as rho signal is usually strong
% (Unless gamma is tiny, but lsqcurvefit handles noise well)
[p_fit, resnorm] = lsqcurvefit(func_rho, [9, 4], X_data, Y_rho, [0,0], [100,100], options);

fprintf('  A (Expect 9): %.4f\n', p_fit(1));
fprintf('  B (Expect 4): %.4f\n', p_fit(2));

%% VISUALIZATION
figure('Name','Tele Fit Analysis','Color','w', 'Position', [100 100 1200 400]);

% Plot q Fit
subplot(1,3,1); 
plot(X_q, Y_q, 'bo', 'MarkerFaceColor','b'); hold on;
plot(X_q, K_q*X_q, 'r-', 'LineWidth', 1.5);
title(['q Fit: K=',num2str(K_q,'%.4f')]); 
xlabel('Model Term (\gamma/c)'); ylabel('Actual q');
grid on; axis square;

% Plot p Fit (Log-Log for better visibility of outliers)
subplot(1,3,2); 
loglog(X_p(valid_p_idx), Y_p(valid_p_idx), 'bo', 'MarkerFaceColor','b', 'DisplayName','Valid');
hold on;
if any(~valid_p_idx)
    loglog(X_p(~valid_p_idx), Y_p(~valid_p_idx), 'kx', 'DisplayName','Noise');
end
% Theory line
x_sort = sort(X_p);
loglog(x_sort, K_p*x_sort, 'r-', 'LineWidth', 1.5, 'DisplayName','Fit');

title(['p Fit: K=',num2str(K_p,'%.4f')]); 
xlabel('Model Term (\gamma^2 / c\omega_{min}^2)'); ylabel('Actual p correction');
grid on; axis square; legend('Location','best');

% Plot Refined Rho Fit
subplot(1,3,3); 
Y_pred = func_rho(p_fit, X_data);
plot(Y_rho, Y_pred, 'ko', 'MarkerFaceColor','k'); 
hold on; 
plot([0 max(Y_rho)], [0 max(Y_rho)], 'r-', 'LineWidth', 1.5);

title(['Refined \rho Model (A=',num2str(p_fit(1),'%.1f'),', B=',num2str(p_fit(2),'%.1f'),')']); 
xlabel('Actual \rho'); ylabel('Formula Prediction');
grid on; axis square;