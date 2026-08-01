function params = getParams()
    %% Physical Parameters
    params.m = [1, 1];
    params.L = [1; 1];
    params.g = 9.81; % Gravitational acceleration

    %% State space parameters
    params.q = [1, 2];                           % Force input

    %% LQR Parameters
    params.epsilon = 1e-5;
    params.Q0 = diag([10, 10, 1, 1]);
    params.R0 = 1;
end