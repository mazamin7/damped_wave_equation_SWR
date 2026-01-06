% Precompute contraction factor function with finite-domain factors
function r = rho(N, s, theta1, theta2, c, gamma, nu, a, b)
    if(N == 2)
        % Spectral radius with finite-domain factors a, b
        ikappa = sqrt((s.^2+gamma.*s) ./ (c^2+nu.*s));
    
        num = ikappa - (theta1*s + theta2);
        den = ikappa + (theta1*s + theta2);
    
        r = num ./ den * exp(-ikappa * (b-a));
    else
        G = matrix_G(N, s, theta1, theta2, c, gamma, nu, a, b);
        r = max(abs(eig(G)));
    end
end