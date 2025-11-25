function interactive_rho_all_params()
    % INTERACTIVE_RHO_ALL_PARAMS
    % Interactive GUI to plot Contraction Factor rho(omega)
    % allowing real-time adjustment of physics parameters (gamma, nu)
    % and transmission parameters (p, q).
    % 
    % FEATURES:
    % - Real-time curve updates.
    % - Automatic calculation of "Best" p,q (Red dashed line).
    % - "Dots" (MajorTicks) on sliders indicating the optimal values.

    close all; clc;
    addpath("utils\")

    %% 1. SETUP CONSTANTS
    c  = 1.0;
    T  = 5.0;
    dt = 0.002;
    Lx = 1.0;

    % Geometry
    delta = 0.1;
    a = (Lx - delta)/2;
    M = delta;        
    b = a + M;        

    % Solver params
    N  = 2;
    ky = 0;
    
    % Frequency Grid (High Res for Plotting)
    omega_min = pi / T;
    omega_max = pi / dt;
    
    J_plot = 500;
    omegas = logspace(log10(omega_min), log10(omega_max), J_plot);
    s_vals = 1i * omegas;

    %% 2. INITIAL VALUES
    init_gamma = 0;
    init_nu    = 0;
    init_p     = 1.0/c;
    init_q     = 0.0;

    %% 3. BUILD UI
    fig = uifigure('Name', 'SWR Parameter Explorer', 'Position', [100, 100, 1000, 700]);

    % Layout: Plot (Top), Controls (Bottom)
    gl = uigridlayout(fig, [2, 1]);
    % INCREASED HEIGHT: Changed from 200 to 350 to fix label overlapping
    gl.RowHeight = {'1x', 350}; 

    % --- AXES ---
    ax = uiaxes(gl);
    ax.XScale = 'log';
    ax.XGrid = 'on';
    ax.YGrid = 'on';
    ax.FontSize = 12;
    title(ax, 'Contraction Factor \rho(\omega)');
    xlabel(ax, 'Frequency \omega [rad/s]');
    ylabel(ax, '|\rho|');
    xlim(ax, [omega_min, omega_max]); 
    ylim(ax, [0, 1.1]);

    % --- PLOT OBJECTS ---
    hold(ax, 'on');
    
    % 1. Optimal Reference (Red Dashed)
    % Calculated initially for 0,0
    r_opt_init = compute_rho_curve(init_gamma, init_nu, init_p, init_q); 
    h_opt_line = plot(ax, omegas, r_opt_init, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Optimal (Auto)');
    
    % 2. Current Manual (Blue Solid)
    r_init = compute_rho_curve(init_gamma, init_nu, init_p, init_q);
    h_line = plot(ax, omegas, r_init, 'b-', 'LineWidth', 2, 'DisplayName', 'Current (Manual)');
    
    legend(ax, 'Location', 'northeast');
    
    % 3. Horizontal Max Line
    max_init = max(r_init);
    h_max_line = yline(ax, max_init, 'b:', 'LineWidth', 1.5, ...
        'Label', sprintf('Max: %.4f', max_init), 'LabelHorizontalAlignment', 'left');

    hold(ax, 'off');

    % --- CONTROLS PANEL ---
    panel = uipanel(gl, 'Title', 'Parameters');
    % Grid: 4 Rows (Gamma, Nu, P, Q), 3 Cols (Label, Slider, Value)
    pgrid = uigridlayout(panel, [4, 3]);
    % Use fixed width (80px) for the value labels to prevent resizing/jumping
    pgrid.ColumnWidth = {'fit', '1x', 80};
    
    % Add spacing between rows
    pgrid.RowSpacing = 15; 

    % Slider 1: Gamma
    uilabel(pgrid, 'Text', 'Gamma (Friction):');
    sld_gamma = uislider(pgrid, 'Limits', [0, 20], 'Value', init_gamma);
    lbl_gamma = uilabel(pgrid, 'Text', sprintf('%.2f', init_gamma));
    sld_gamma.MajorTicks = []; % Clean ticks

    % Slider 2: Nu
    uilabel(pgrid, 'Text', 'Nu (Viscosity):');
    sld_nu = uislider(pgrid, 'Limits', [0, 0.1], 'Value', init_nu);
    lbl_nu = uilabel(pgrid, 'Text', sprintf('%.4f', init_nu));
    sld_nu.MajorTicks = []; % Clean ticks

    % Slider 3: p
    uilabel(pgrid, 'Text', 'p (Time Deriv):');
    sld_p = uislider(pgrid, 'Limits', [0, 2], 'Value', init_p);
    lbl_p = uilabel(pgrid, 'Text', sprintf('%.4f', init_p));

    % Slider 4: q
    uilabel(pgrid, 'Text', 'q (Constant):');
    sld_q = uislider(pgrid, 'Limits', [-5, 20], 'Value', init_q);
    lbl_q = uilabel(pgrid, 'Text', sprintf('%.4f', init_q));

    %% 4. UPDATE LOGIC
    % Store last optimal p/q to avoid re-optimizing if only p/q change
    last_g = -999; 
    last_n = -999;
    
    function updatePlot(src, event)
        % Get current values (handle dragging)
        g_val = sld_gamma.Value;
        n_val = sld_nu.Value;
        p_val = sld_p.Value;
        q_val = sld_q.Value;

        if strcmp(event.EventName, 'ValueChanging')
            new_val = event.Value;
            if src == sld_gamma, g_val = new_val;
            elseif src == sld_nu, n_val = new_val;
            elseif src == sld_p, p_val = new_val;
            elseif src == sld_q, q_val = new_val;
            end
        end

        % Update Labels
        lbl_gamma.Text = sprintf('%.2f', g_val);
        lbl_nu.Text    = sprintf('%.4f', n_val);
        lbl_p.Text     = sprintf('%.4f', p_val);
        lbl_q.Text     = sprintf('%.4f', q_val);

        % --- OPTIMIZATION LOGIC ---
        % STRICT CHECK: Only re-optimize if PHYSICS sliders are the source
        % This prevents P/Q dragging from triggering expensive re-optimization or tick artifacts
        physics_changed = (src == sld_gamma) || (src == sld_nu) || strcmp(event.EventName, 'Init');
        
        if physics_changed && ((abs(g_val - last_g) > 1e-4) || (abs(n_val - last_n) > 1e-5))
            
            [opt_p, opt_q] = get_optimum(g_val, n_val);
            
            % Update Red Curve
            rho_opt = compute_rho_curve(g_val, n_val, opt_p, opt_q);
            h_opt_line.YData = rho_opt;
            
            % Update Slider "Dots" (MajorTicks)
            % Only show tick if it's within the slider's range to avoid visual bugs
            if opt_p >= sld_p.Limits(1) && opt_p <= sld_p.Limits(2)
                sld_p.MajorTicks = opt_p;
                sld_p.MajorTickLabels = {sprintf('★ %.2f', opt_p)};
            else
                sld_p.MajorTicks = [];
            end
            
            if opt_q >= sld_q.Limits(1) && opt_q <= sld_q.Limits(2)
                sld_q.MajorTicks = opt_q;
                sld_q.MajorTickLabels = {sprintf('★ %.2f', opt_q)};
            else
                sld_q.MajorTicks = [];
            end
            
            last_g = g_val;
            last_n = n_val;
        end

        % --- CURRENT CURVE ---
        % Recompute and Update Blue Plot (Always happens)
        new_rho = compute_rho_curve(g_val, n_val, p_val, q_val);
        h_line.YData = new_rho;
        
        % Update Horizontal Line
        max_rho = max(new_rho);
        h_max_line.Value = max_rho;
        h_max_line.Label = sprintf('Max: %.4f', max_rho);
        
        drawnow limitrate;
    end

    % Bind Callbacks
    sld_gamma.ValueChangingFcn = @updatePlot; sld_gamma.ValueChangedFcn = @updatePlot;
    sld_nu.ValueChangingFcn    = @updatePlot; sld_nu.ValueChangedFcn    = @updatePlot;
    sld_p.ValueChangingFcn     = @updatePlot; sld_p.ValueChangedFcn     = @updatePlot;
    sld_q.ValueChangingFcn     = @updatePlot; sld_q.ValueChangedFcn     = @updatePlot;

    % Trigger initial update to set ticks
    updatePlot(sld_gamma, struct('EventName','Init'));

    %% 5. PHYSICS HELPER
    function r_vec = compute_rho_curve(g, n, p, q)
        r_vec = zeros(1, J_plot);
        % Use library rho function
        for k = 1:J_plot
            val = rho(N, s_vals(k), p, q, c, g, n, a, b, ky);
            r_vec(k) = abs(val);
        end
    end

    %% 6. OPTIMIZATION HELPER (Lightweight)
    function [op, oq] = get_optimum(g, n)
        % Use a coarse grid for real-time responsiveness
        J_fast = 50; 
        
        % Initial guess based on asymptotics
        p0 = 1.0/c;
        q0 = g/(2*c); 
        
        % Use library obj_Linf (Minimax) with fminsearch
        % We assume obj_Linf is in path. 
        % Note: obj_Linf usually takes J. We pass J_fast.
        
        fun = @(x) obj_Linf(N, T, dt, J_fast, c, g, n, a, M, x(1), x(2), ky);
        
        % Fast options
        opts_fast = optimset('Display','off', 'TolX',1e-3, 'TolFun',1e-3);
        
        x_res = fminsearch(fun, [p0, q0], opts_fast);
        op = x_res(1);
        oq = x_res(2);
    end
end