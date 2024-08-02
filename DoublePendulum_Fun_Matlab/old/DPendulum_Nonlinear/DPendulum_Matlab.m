clear all; close all; clc;
double_pendulum_simulation

function double_pendulum_simulation
    % Parameters
    m1 = 1;  % Mass of the first pendulum
    m2 = 1;  % Mass of the second pendulum
    l1 = 1;  % Length of the first pendulum
    l2 = 1;  % Length of the second pendulum
    g = 9.81; % Gravitational acceleration

    % Initial conditions [theta1, theta2, theta1_dot, theta2_dot]
    initial_conditions = [pi, pi/4, 0, 0];

    % Time span for the simulation
    tspan = [0 5];

    % Solve the ODEs
    options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);  % Relaxing the tolerances
    [t, y] = ode45(@(t, y) double_pendulum_ode(t, y, m1, m2, l1, l2, g), tspan, initial_conditions);

    % Plot the results
    subplot(2,1,1);
    plot(t, y(:, 1), 'r', t, y(:, 2), 'b');
    title('Double Pendulum Simulation');
    xlabel('Time (s)');
    ylabel('Angle (rad)');
    legend('$\theta_1$', '$\theta_2$', 'Interpreter', 'latex');
    
    subplot(2,1,2);
    plot(t, y(:, 3), 'r', t, y(:, 4), 'b');
    xlabel('Time (s)');
    ylabel('Angular Velocity (rad/s)');
    legend('$\dot{\theta}_1$', '$\dot{\theta}_2$', 'Interpreter', 'latex');

    % Animate the results
    animate_double_pendulum(t, y, l1, l2);
end

function dydt = double_pendulum_ode(t, y, m1, m2, l1, l2, g)
    % Unpack the state vector
    theta_1 = y(1);
    theta_2 = y(2);
    theta_dot_1 = y(3);
    theta_dot_2 = y(4);

    % Calculate the derivatives
    alpha_1 = m2/(m1+m2) * l1/l2 * cos(theta_1-theta_2);
    alpha_2 = l1/l2 * cos(theta_1-theta_2);
    A = [1, alpha_1; alpha_2, 1];

    beta_1 = - m2/(m1+m2) * l1/l2 * theta_dot_2^2 * sin(theta_1-theta_2) - g/l1 * sin(theta_1);

    beta_2 = theta_dot_1^2 * l1/l2 * sin(theta_1-theta_2) - g/l2 * sin(theta_2);

    B = [beta_1;beta_2];
    T = [0;0];
    B = B+T;
    acceleration = A\B;
    
    dydt = [theta_dot_1; theta_dot_2; acceleration];
end

function animate_double_pendulum(t, y, l1, l2)
    % Create a figure for the animation
    figure;
    hold on;
    axis equal;
    axis([-2*(l1+l2) 2*(l1+l2) -2*(l1+l2) 2*(l1+l2)]);
    plot([-2*(l1+l2) 2*(l1+l2)], [0 0], 'k'); % ground line
    
    % Loop through the time steps to create the animation
    for i = 1:length(t)
        theta1 = y(i, 1);
        theta2 = y(i, 2);
        
        % Positions of the pendulums
        x1 = l1 * sin(theta1);
        y1 = -l1 * cos(theta1);
        x2 = x1 + l2 * sin(theta2);
        y2 = y1 - l2 * cos(theta2);
        
        % Plot the pendulum
        plot([0, x1], [0, y1], 'r', 'LineWidth', 2); % first rod
        plot([x1, x2], [y1, y2], 'b', 'LineWidth', 2); % second rod
        plot(x1, y1, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r'); % first mass
        plot(x2, y2, 'bo', 'MarkerSize', 10, 'MarkerFaceColor', 'b'); % second mass
        
        drawnow;
        
        % Clear the figure except the ground line for the next frame
        if i < length(t)
            clf;
            hold on;
            axis equal;
            axis([-2*(l1+l2) 2*(l1+l2) -2*(l1+l2) 2*(l1+l2)]);
            plot([-2*(l1+l2) 2*(l1+l2)], [0 0], 'k'); % ground line
        end
    end
end