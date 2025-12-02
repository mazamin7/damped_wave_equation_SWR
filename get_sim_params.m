function P = get_sim_params()
    % Global simulation parameters for the 1D SWR project

    % Domain / geometry / decomposition
    P.N  = 2;
    P.a  = 0.3;
    P.M  = 0.1;
    P.b  = P.a + P.M;
	P.Lx = P.N * P.a + P.M;
	P.T  = 5.0;

    % Physical parameters
    P.c  = 1.0;
	P.gamma = 0.0;
	P.nu = 0.0;

    % Time / space / frequency discretization
    P.dh = 0.002;
    P.dt = 0.002;
    P.J  = 1000;
end
