function [K_lqr_opt, Q_opt, R_opt] = optimize_lqr(sys, params)
    % Initial guesses for Q and R
    Q0 = params.Q0;
    R0 = params.R0;
    epsilon = params.epsilon;
    options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'sqp', 'MaxIterations', 50);

    % Initial guess and boundaries
    x0 = [diag(Q0); R0];            % Initial guess for optimization variables
    lb = epsilon * ones(size(x0));  % Lower bounds to ensure positive definiteness
    ub = Inf(size(x0));             % Upper bounds

    % Objective function
    obj_fun = @(x) lqr_objective(x, sys, params);
    
    % Check the initial guess for the objective function
    initial_cost = obj_fun(x0)
    if isinf(initial_cost)
        error('Initial guess leads to an infeasible solution. Adjust Q0 and R0.');
    end

    % Optimize
    [x_opt] = fmincon(obj_fun, x0, [], [], [], [], lb, ub, [], options);

    % Extract optimized Q and R
    Q_opt = diag(x_opt(1:4));
    R_opt = x_opt(5);

    % Compute optimized LQR gain matrix
    K_lqr_opt = lqr(sys.A, sys.B, Q_opt, R_opt);
end
