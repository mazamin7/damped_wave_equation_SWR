function G = matrix_G(N, s, theta1, theta2, c, gamma, nu, a, b)
% COMPUTE_TRANSFER_MATRIX Compute the transfer matrix G for given parameters
% red black SWR algorithm for N subdomains of length b = a+M and overlap M
%
% Inputs:
%   N       - Number of segments
%   s       - Laplace variable
%   theta1  - First interface parameter
%   theta2  - Second interface parameter  
%   c       - Wave speed
%   gamma   - First damping coefficient
%   nu      - Second damping coefficient
%   a       - Length of the first portion of a subdomain
%   b       - Total length of a subdomain
%
% Output:
%   G       - Transfer matrix (2N x 2N)

% Calculate M from a and b
M = b - a;

% Define position functions
aj = @(j) a*(j-1);
bj = @(j) aj(j+1) + M;

Lx = bj(N);

% Calculate kx
kx = -1i * sqrt((s.^2+gamma.*s) ./ (c^2+nu.*s));

interface_term = s * theta1 + theta2;

% Define matrix elements
s_111 = exp(1i * kx * bj(1)) * (1i * kx + interface_term) + exp(-1i * kx * bj(1)) * (1i * kx - interface_term);
s_122 = -s_111;

s_j11 = @(j) exp(1i * kx * aj(j)) * (1i * kx - interface_term);
s_j12 = @(j) -exp(-1i * kx * aj(j)) * (1i * kx + interface_term);
s_j21 = @(j) exp(1i * kx * bj(j)) * (1i * kx + interface_term);
s_j22 = @(j) -exp(-1i * kx * bj(j)) * (1i * kx - interface_term);

s_N11 = exp(1i * kx * aj(N)) * (1i * kx - interface_term) + exp(-1i * kx * aj(N)) * exp(2 * 1i * kx * Lx) * (1i * kx + interface_term);
s_N22 = -s_N11 * exp(-2 * 1i * kx * Lx);

S_block_1 = [s_111, 0;
             0, s_122];

S_block_j = @(j) [s_j11(j), s_j12(j);
                  s_j21(j), s_j22(j)];

S_block_N = [s_N11, 0;
             0, s_N22];

v_111 = exp(1i * kx * bj(1)) * (1i * kx + interface_term);
v_112 = -exp(-1i * kx * bj(1)) * (1i * kx - interface_term);
v_121 = v_111;
v_122 = v_112;

v_j11 = @(j) s_j11(j);
v_j12 = @(j) s_j12(j);
v_j21 = @(j) s_j21(j);
v_j22 = @(j) s_j22(j);

v_N11 = exp(1i * kx * aj(N)) * (1i * kx - interface_term);
v_N12 = -exp(-1i * kx * aj(N)) * (1i * kx + interface_term);
v_N21 = v_N11;
v_N22 = v_N12;

V_block_1 = [v_111, v_112;
             v_121, v_122];

V_block_jL = @(j) [v_j11(j), v_j12(j);
                   0, 0];

V_block_jR = @(j) [0, 0;
                   v_j21(j), v_j22(j)];

V_block_N = [v_N11, v_N12;
             v_N21, v_N22];

% Build S matrix
S = zeros(2*N);

S(1:2,1:2) = S_block_1;

for j = 2:N-1
    if mod(j,2) == 0
        S(2*j-1:2*j,2*j-3:2*j-2) = -V_block_jL(j);
        S(2*j-1:2*j,2*j+1:2*j+2) = -V_block_jR(j);
    end
    S(2*j-1:2*j,2*j-1:2*j) = S_block_j(j);
end

if mod(N,2) == 0
    S(end-1:end,end-3:end-2) = -V_block_N;
end
S(end-1:end,end-1:end) = S_block_N;

% Build V matrix
V = zeros(2*N);

V(1:2,3:4) = V_block_1;

for j = 2:N-1
    if mod(j,2) == 1
        V(2*j-1:2*j,2*j-3:2*j-2) = V_block_jL(j);
        V(2*j-1:2*j,2*j+1:2*j+2) = V_block_jR(j);
    end
end

if mod(N,2) == 1
    V(end-1:end,end-3:end-2) = V_block_N;
end

% Compute transfer matrix G
warning('off', 'MATLAB:nearlySingularMatrix');
G = S \ V;
warning('on', 'MATLAB:nearlySingularMatrix');

end