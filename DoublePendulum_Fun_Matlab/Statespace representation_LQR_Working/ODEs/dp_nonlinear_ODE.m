function dydt = dp_nonlinear_ODE(t, y,sys ,params)
    % Unpack the state vector
    theta_1 = y(1);
    theta_2 = y(2);
    theta_dot_1 = y(3);
    theta_dot_2 = y(4);
    delta = theta_1-theta_2;
    m = params.m;
    L = params.L;
    g = params.g;
    m1 = m(1);
    m2 = m(2);
    l1 = L(1);
    l2 = L(2);


    if isfield(params, 'K_lqr')
    K_lqr = params.K_lqr;
    desired_angles = params.desired_angles;
    u = -K_lqr * (y - [desired_angles; 0; 0]);
    else
        u = [0;0];
    end
    % Matrices
    M = [(m1+m2)*l1     , m2*l2*cos(delta); l1*cos(delta)  , l2]; % Mass matrix
    C = [0              , m2*l2*sin(delta); l1*sin(delta)  , 0];  % Coupling matrix
    G = [(m1+m2)*g*sin(theta_1); g*sin(theta_2)];                 % Gravitational matrix
    
    ddtheta = M\(-C*y(3:4).^2-G+u);                                 % Acceleration
    
    dydt = [theta_dot_1; theta_dot_2; ddtheta];
end 