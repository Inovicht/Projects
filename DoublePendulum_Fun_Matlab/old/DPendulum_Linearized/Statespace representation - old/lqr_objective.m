function cost = lqr_objective(x, sys, init, tspan, desired_angles)
    % Extract Q and R from optimization variables
    Q = diag(x(1:4));
    R = x(5);
    
    % Ensure Q and R are positive definite
    if min(eig(Q)) <= 0 || R <= 0
        disp('Q or R not positive definite.');
        cost = Inf;
        return;
    end
    
    % Regularization term to ensure numerical stability
    Q = Q + 1e-5 * eye(size(Q));
    R = R + 1e-5;
    
    % Compute LQR gain matrix
    try
        K_lqr = lqr(sys.A, sys.B, Q, R);
    catch
        disp('LQR computation failed.');
        cost = Inf;
        return;
    end
    
    % Simulate the system with LQR control
    try
        [t, y] = ode45(@(t, y) double_pendulum_ode(t, y, sys, K_lqr, desired_angles), tspan, init);
    catch
        disp('ODE simulation failed.');
        cost = Inf;
        return;
    end
    
    % Calculate the tracking error
    tracking_error = y(:, 1:2) - desired_angles';  % Error in theta_1 and theta_2
    
    % Define performance metrics (e.g., integrated square error)
    cost = trapz(t, sum(tracking_error.^2, 2));
end
