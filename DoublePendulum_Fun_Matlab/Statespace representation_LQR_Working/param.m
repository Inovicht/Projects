function params = param
%% Physical Parameters
params.m = [1,1];
params.L = [1;1];
params.g = 9.81; % Gravitational acceleration

%% State space parameters
% Mapping the output vector (Which node/s)
params.d       =   [1];                     % Displacement output
params.ddot    =   [];                      % Velocity output
params.dddot   =   [];                      % Acceleration output
% Mapping the input vector
params.q = [1,2];                           % Force input

%% LQR Parameters
params.epsilon = 1e-5;
params.Q0 = diag([10, 10, 1, 1]);
params.R0 = 1;
end

