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
