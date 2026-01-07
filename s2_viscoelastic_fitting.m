%% SCRIPT 2: VISCOELASTIC COEFFICIENT FITTING
% Fits the constants Kp, Kq, and K_rho with Robust Outlier Rejection
clear; close all; clc; addpath("utils\");

c = 1.0; gamma = 0; N = 2; a = 0.5; M = 0; J = 1000;
% Tight tolerances required because p-correction is very small (nu^2)
options = optimset('Display','off', 'TolX',1e-14, 'TolFun',1e-14);

fprintf('--- VISCOELASTIC COEFFICIENT FITTING ---\n');
NumPoints = 100;

% 1. Logarithmic Random Sampling
% nu: 10^-7 to 10^-5
nu_rand    = 10.^( -7 + 2 * rand(NumPoints, 1) ); 
% w_max: 10 to 40
w_max_rand = 10.^( 1 + (log10(40)-1) * rand(NumPoints, 1) );
% w_min: 0.1 to 1.0
w_min_rand = 10.^( -1 + 1 * rand(NumPoints, 1) );

Y_q = zeros(NumPoints,1);
Y_p = zeros(NumPoints,1);
Y_rho = zeros(NumPoints,1);

for i=1:NumPoints
    nu = nu_rand(i); w_max = w_max_rand(i); w_min = w_min_rand(i);
    T = pi/w_min; dt = pi/w_max;
    
    % Smart Initial Guess (Theoretical)
    q_guess = 0.5 * nu * w_min * w_max / c^3;
    p_guess = 1/c - 0.25 * nu^2 * w_max^2 / c^5;
    
    obj = @(x) obj_Linf(N,T,dt,J,c,gamma,nu,a,M,x(1),x(2));
    [res, fval] = fminsearch(obj, [p_guess, q_guess], options);
    
    Y_p(i) = 1/c - res(1); % The correction drop
    Y_q(i) = res(2);
    Y_rho(i) = fval;
    
    if mod(i, 20) == 0, fprintf('Processed %d/%d...\n', i, NumPoints); end
end

% --- FIT MODELS ---

% 1. Fit q (Linear in nu -> Strong signal)
X_q = (nu_rand .* w_min_rand .* w_max_rand) / c^3;
K_q = X_q \ Y_q;

% 2. Fit p (Quadratic in nu -> Weak signal, needs filtering)
X_p = (nu_rand.^2 .* w_max_rand.^2) / c^5;

% Filter Outliers:
% Remove points where optimization failed to move (noise) 
% or produced non-physical negative drops.
% Threshold 1e-13 is chosen just above typical double-precision noise for these operations.
valid_p_idx = (Y_p > 1e-13) & (Y_p < 1.0); 

if sum(valid_p_idx) < 5
    warning('Too few valid points for p-fit. Try increasing nu range.');
    K_p = 0;
else
    K_p = X_p(valid_p_idx) \ Y_p(valid_p_idx);
end

% 3. Fit rho (Linear in nu -> Strong signal)
X_rho = (nu_rand .* (w_max_rand - w_min_rand)) / c^2;
K_rho = X_rho \ Y_rho;

fprintf('\nFINAL RESULTS:\n');
fprintf('  K_q   (Expect 0.50): %.4f\n', K_q);
fprintf('  K_p   (Expect 0.25): %.4f  (Fitted on %d valid points)\n', K_p, sum(valid_p_idx));
fprintf('  K_rho (Expect 0.25): %.4f\n', K_rho);

% --- VISUALIZATION ---
figure('Name','Visco Fit','Color','w', 'Position', [100 100 1200 400]);

% Plot q
subplot(1,3,1); 
plot(X_q, Y_q, 'bo', 'MarkerFaceColor','b'); hold on; 
plot(X_q, K_q*X_q, 'r-', 'LineWidth', 1.5);
title(['q Fit: K=',num2str(K_q,'%.4f')]); xlabel('Model Term'); ylabel('Actual q');
grid on; axis square;

% Plot p (Highlighting Valid vs Outliers)
subplot(1,3,2); 
plot(X_p(valid_p_idx), Y_p(valid_p_idx), 'bo', 'MarkerFaceColor','b', 'DisplayName', 'Valid Data'); 
hold on; 
if any(~valid_p_idx)
    plot(X_p(~valid_p_idx), Y_p(~valid_p_idx), 'kx', 'DisplayName', 'Outliers (Noise)');
end
plot(X_p, K_p*X_p, 'r-', 'LineWidth', 1.5, 'DisplayName', 'Fit');
title(['p Fit: K=',num2str(K_p,'%.4f')]); xlabel('Model Term'); ylabel('Actual p drop');
legend('Location','best'); grid on; axis square;

% Plot rho
subplot(1,3,3); 
plot(X_rho, Y_rho, 'ko', 'MarkerFaceColor','k'); hold on; 
plot(X_rho, K_rho*X_rho, 'r-', 'LineWidth', 1.5);
title(['\rho Fit: K=',num2str(K_rho,'%.4f')]); xlabel('Model Term'); ylabel('Actual \rho');
grid on; axis square;