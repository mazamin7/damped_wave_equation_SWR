% Precompute contraction factor function with finite-domain factors
function r = rho(N, s, theta1, theta2, c, gamma, nu, a, b)
    if(N == 2)
        % Spectral radius with finite-domain factors a, b
        ikappa = sqrt((s.^2+gamma.*s) ./ (c^2+nu.*s));

        exp_pa = 1;% exp(ikappa * a);
        exp_ma = 1;% exp(-ikappa * a);
        exp_pb = 1;% exp(ikappa * b);
        exp_mb = 1;% exp(-ikappa * b);
    
        num = ikappa .* (exp_pa + exp_ma) - (theta1*s + theta2) .* (exp_pa - exp_ma);
        den = ikappa .* (exp_pb + exp_mb) + (theta1*s + theta2) .* (exp_pb - exp_mb);
    
        r = num ./ den * exp(-ikappa * (b-a));
    else
        G = matrix_G(N, s, theta1, theta2, c, gamma, nu, a, b);
        r = max(abs(eig(G)));
    end
end