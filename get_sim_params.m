function P = get_sim_params()
    % Global simulation parameters for the 1D SWR project

    % Domain / geometry / decomposition
    P.N  = 2;
    % P.a  = 0.3;
    P.a = 2;
	P.T  = 5.0;

    % Physical parameters
    P.c  = 1.0;
	P.gamma = 0.0;
	P.nu = 0.0;

    % Time / space / frequency discretization
    P.dh = 0.02;
    P.dt = 0.02;
    P.J  = 1000;

    P.M = 10*P.dh;
    P.b = P.a + P.M;
    P.Lx = P.N * P.a + P.M;
end
