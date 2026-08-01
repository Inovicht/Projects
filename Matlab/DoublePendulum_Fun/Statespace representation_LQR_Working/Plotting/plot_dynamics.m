function total_time = plot_dynamics(t, y, end_time)
    % Find the indices of t that are within the range from t(1) to end_time
    idx = t <= end_time;
    
    % Filter t and y to only include data up to end_time
    t_filtered = t(idx);
    y_filtered = y(idx, :);
    
    % Calculate total simulation time
    total_time = end_time - t(1);
    
    % Plot the results
    subplot(2,1,1);
    plot(t_filtered, y_filtered(:, 1), 'r', t_filtered, y_filtered(:, 2), 'b');
    ylim([-pi-0.2,pi+0.2])
    grid on
    title('Double Pendulum Simulation');
    xlabel('Time (s)');
    ylabel('Angle (rad)');
    legend('$\theta_1$', '$\theta_2$', 'Interpreter', 'latex');
    
    subplot(2,1,2);
    plot(t_filtered, y_filtered(:, 3), 'r', t_filtered, y_filtered(:, 4), 'b');
    grid on
    xlabel('Time (s)');
    ylabel('Angular Velocity (rad/s)');
    legend('$\dot{\theta}_1$', '$\dot{\theta}_2$', 'Interpreter', 'latex');

    % Display the total simulation time
    disp(['Total Simulation Time: ', num2str(total_time), ' seconds']);
end