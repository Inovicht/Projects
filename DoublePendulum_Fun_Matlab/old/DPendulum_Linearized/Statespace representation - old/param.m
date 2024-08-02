function params = param
% Parameters
params.m = [5,1];
params.L = [3;5];
params.g = 9.81; % Gravitational acceleration

%% State space parameters
% Mapping the output vector (Which node/s)
params.d       =   [1,2];   % Displacement output
params.ddot    =   [];                      % Velocity output
params.dddot   =   [];                      % Acceleration output
% Mapping the input vector
params.q = [1];                             % Force input

L = params.L;
m = params.m;

end

