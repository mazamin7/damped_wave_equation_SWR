% Precompute contraction factor function with finite-domain factors
function r = rho(N, s, theta1, theta2, c, gamma, nu, a, b, ky)
    if isnan(ky)
        ky = 0;
        disp('Error: ky was NaN')
    end

    if(N == 2)
        % Spectral radius with finite-domain factors a, b
        ikappa = sqrt((s.^2+gamma.*s) ./ (c^2+nu.*s) + ky^2);

        exp_pa = exp(ikappa * a);
        exp_ma = exp(-ikappa * a);
        exp_pb = exp(ikappa * b);
        exp_mb = exp(-ikappa * b);
    
        num = ikappa .* (exp_pa + exp_ma) - (theta1*s + theta2) .* (exp_pa - exp_ma);
        den = ikappa .* (exp_pb + exp_mb) + (theta1*s + theta2) .* (exp_pb - exp_mb);
    
        r = num ./ den;
    else
        % error('Not implemented for N != 2');
        ky = 0;
        G = matrix_G(N, s, theta1, theta2, c, gamma, nu, a, b, ky);
        r = max(abs(eig(G)));
    end
end