function dydt = double_pendulum_ode(t, y, sys, K_lqr, desired_angles)
    % State vector
    x = y;

    % Control input (LQR)
    u = -K_lqr * (x - [desired_angles; 0; 0]);  % Incorporate desired angles into control law

    % State-space equation
    dydt = sys.A * x + sys.B * u;
end