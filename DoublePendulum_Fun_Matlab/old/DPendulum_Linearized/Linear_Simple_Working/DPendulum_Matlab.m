clear all; close all; clc;
% Parameters
m = [5,1];
L = [3;5];
g = 9.81; % Gravitational acceleration
init = [pi/4, pi/4, 0, 0];
tspan = [0 100];

% Mass and stiffness matrix
M = [L(1)^2*m(1) + L(1)^2*m(2) , L(1)*L(2)*m(2);
     L(1)*L(2)*m(2)          , L(2)^2*m(2)];
K = g*[L(1)*m(2)+L(1)*m(2)     , 0
       0               , L(2)*m(2)];

Dpendulum_simulation(M,K,L,init,tspan)

function double_pendulum_simulation = Dpendulum_simulation(M,K,L,init,tspan)
    % Solve the ODEs
    omega_n_rad_F = [0,0];
    options = odeset('RelTol', 1e-17, 'AbsTol', 1e-17);  % Relaxing the tolerances
    [t, y] = ode45(@(t, y) double_pendulum_ode(t, y, M,K,omega_n_rad_F), tspan, init);
    AmpFrec_analysis(t,y);
    compute_omega_n(M,K);
   % animate_double_pendulum(t, y, L);
end

function [omega_n_Hz, eigenvectors] = compute_omega_n(M,K)
    % Solve the generalized eigenvalue problem
    [eigenvectors, eigenvalues] = eig(K, M);
    
    % Extract the eigenvalues from the diagonal matrix
    eigenvalues = diag(eigenvalues);
    
    % Compute the natural frequencies (in rad/s)
    natural_frequencies_rad_per_s = sqrt(eigenvalues);
    
    % Convert the natural frequencies to Hertz (Hz)
    omega_n_Hz = natural_frequencies_rad_per_s / (2 * pi);
    
    % Display the natural frequencies in Hz
    disp('The natural frequencies (Hz) of the system are:');
    disp(omega_n_Hz);
    
    % Display the eigenvectors (mode shapes)
    disp('The eigenvectors (mode shapes) of the system are:');
    disp(eigenvectors);
end

function AmpFrec_analysis(t,y)
    [frec, amp] = amplitudeSpectrumOneSided(t,y);
    semilogy(frec,amp);
end

function dydt = double_pendulum_ode(t, y, M,K,frequency_Force_rad)
    % Unpack the state vector
    theta_1 = y(1);
    theta_2 = y(2);
    theta_dot_1 = y(3);
    theta_dot_2 = y(4);

    % Calculate the derivatives
    theta = [theta_1;theta_2];
    T = [sin(frequency_Force_rad(1) * t); sin(frequency_Force_rad(2) * t)];
    ddtheta= M\(-K*theta+T);

    dydt = [theta_dot_1; theta_dot_2; ddtheta];
end

function animate_double_pendulum(t, y, L)
    l1 = L(1);
    l2 = L(2);
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

function plot_dynamics(t,y)
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
end



