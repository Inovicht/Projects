function dydt = dp_linear_ODE(t, y, sys, params)
    % Check if K_lqr is provided
    if isfield(params, 'K_lqr')
        K_lqr = params.K_lqr;
        desired_angles = params.desired_angles;
        
        % Control input (LQR)
        u = -K_lqr * (y - [desired_angles; 0; 0]);
    else
        % Default control input (no LQR)
        u = [0;0];
    end

    % State-space equation
    dydt = sys.A * y + sys.B * u;
end