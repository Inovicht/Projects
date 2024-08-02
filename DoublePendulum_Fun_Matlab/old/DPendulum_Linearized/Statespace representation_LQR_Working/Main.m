clear all; close all; clc;
addpath([cd,'/State_Space_Matrix/'])
addpath([cd,'/LQR/'])
addpath([cd,'/ODEs/'])
addpath([cd,'/Plotting/'])

% Define parameters
params = param;
params.init = [pi/4, -pi/4, 0, 0];
params.tspan = [0, 1];
params.desired_angles = [-pi/4 ; pi/4];

% Retrieving state space matrices
sys = Sys(params);

% Optimize the LQR controller
[params.K_lqr_opt, Q_opt, R_opt] = optimize_lqr(sys, params);
params.K_lqr = params.K_lqr_opt;

% Simulate with optimized LQR controller
[t, y] = Dpendulum_simulation(sys, params);
plot_dynamics(t, y,1);
%animate_double_pendulum(t, y, params,10);

function [t, y] = Dpendulum_simulation(sys, params)
  % options = odeset('RelTol', 1e-13, 'AbsTol', 1e-13);
    [t, y] = ode45(@(t, y) dp_nonlinear_ODE(t, y, sys, params), params.tspan, params.init);
end

