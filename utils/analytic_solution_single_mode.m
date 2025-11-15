function u_an = analytic_solution_single_mode(x, t, c, gamma, nu, k0, A0, v0amp)

    % natural angular frequency
    omega0 = c * k0;

    % effective damping
    ge   = 0.5*(gamma + nu*k0^2);
    disc = ge^2 - omega0^2;

    if disc < -1e-14
        omegad = sqrt(omega0^2 - ge^2);
        Ct = cos(omegad*t);
        St = sin(omegad*t);
        q  = exp(-ge*t) .* ( A0*Ct + ((v0amp + ge*A0)/omegad) * St );
    elseif abs(disc) <= 1e-14
        q  = exp(-ge*t) .* ( A0 + (v0amp + ge*A0)*t );
    else
        s  = sqrt(disc);
        r1 = -ge + s;
        r2 = -ge - s;
        C1 = (v0amp - r2*A0)/(r1 - r2);
        C2 = (r1*A0 - v0amp)/(r1 - r2);
        q  = C1*exp(r1*t) + C2*exp(r2*t);
    end

    u_an = sin(k0 * x(:)) * q;
end
