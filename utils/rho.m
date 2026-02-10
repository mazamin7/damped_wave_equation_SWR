% Precompute contraction factor function with finite-domain factors
function r = rho(N, s, theta1, theta2, c, gamma, nu, a, b)
    if(N == 2)
        % Spectral radius with finite-domain factors a, b
        % ikappa = sqrt((s.^2+gamma.*s) ./ (c^2+nu.*s));
        ikappa = s/c .* sqrt((1+gamma./s) ./ (1+nu.*s/c^2));
    
        num = ikappa - (theta1*s + theta2);
        den = ikappa + (theta1*s + theta2);
    
        % r = num ./ den * exp(-ikappa * (b-a));
        r = num ./ den;
    else
        G = matrix_G(N, s, theta1, theta2, c, gamma, nu, a, b);
        r = max(abs(eig(G)));
    end
end