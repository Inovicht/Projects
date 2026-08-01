clear all; close all; clc;
addpath([cd,'/State_Space_Matrix/'])
params = param;
params.init = [0, 0, 0, 0];
params.tspan = [0 1];
% Desired angles
params.desired_angles = [pi/4; pi/4];  % Target angles for theta_1 and theta_2

sys = Sys(params);

[params.K_lqr_opt, params.Q_opt, params.R_opt] = optimize_lqr(sys, params.init, params.tspan, params.desired_angles);

% Simulate with optimized LQR controller
[t,y] = Dpendulum_simulation(sys, params);
plot_dynamics(t,y)
% AmpFrec_analysis(t, y);
% compute_omega_n(M, K); % Optional: Compute natural frequencies
animate_double_pendulum(t, y, params.L); % Optional: Animate double pendulum
% plot_dynamics(t,y)

function [t,y] = Dpendulum_simulation(sys, params)
    L = params.L;
    init = params.init;
    tspan = params.tspan;
    K_lqr = params.K_lqr_opt;
    desired_angles = params.desired_angles;
    options = odeset('RelTol', 1e-8, 'AbsTol', 1e-10); % Relaxing the tolerances
    [t, y] = ode45(@(t, y) double_pendulum_ode(t, y, sys, K_lqr, desired_angles), tspan, init);
end