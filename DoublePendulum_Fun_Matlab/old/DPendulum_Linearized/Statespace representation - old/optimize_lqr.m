function [K_lqr_opt, Q_opt, R_opt] = optimize_lqr(sys, init, tspan, desired_angles)
    % Initial guesses for Q and R
    Q0 = diag([10, 10, 1, 1]);
    R0 = 1;

    epsilon = 1e-5;

    % Optimization options
    options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'sqp', 'MaxIterations', 1000);

    % Optimization problem
    x0 = [diag(Q0); R0]; % Initial guess for optimization variables
    lb = epsilon * ones(size(x0)); % Lower bounds to ensure positive definiteness
    ub = Inf(size(x0));   % Upper bounds

    % Objective function
    obj_fun = @(x) lqr_objective(x, sys, init, tspan, desired_angles);

    % Display initial values
    disp('Initial guess for optimization variables:');
    disp(x0);

    % Check the initial guess for the objective function
    initial_cost = obj_fun(x0);
    if isinf(initial_cost)
        error('Initial guess leads to an infeasible solution. Adjust Q0 and R0.');
    end

    % Optimize
    x_opt = fmincon(obj_fun, x0, [], [], [], [], lb, ub, [], options);

    % Extract optimized Q and R
    Q_opt = diag(x_opt(1:4));
    R_opt = x_opt(5);

    % Compute optimized LQR gain matrix
    K_lqr_opt = lqr(sys.A, sys.B, Q_opt, R_opt);
end
