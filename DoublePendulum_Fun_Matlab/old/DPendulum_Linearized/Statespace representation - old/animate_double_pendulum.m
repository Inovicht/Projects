function animate_double_pendulum(t, y, L)
    l1 = L(1);
    l2 = L(2);
    % Create a figure for the animation
    figure;
    hold on;
    axis equal;
    axis([-2*(l1+l2) 2*(l1+l2) -2*(l1+l2) 2*(l1+l2)]);
    plot([-2*(l1+l2) 2*(l1+l2)], [0 0], 'k', 'LineWidth', 1); % ground line
    
    % Create time annotation
    time_text = text(-1.5*(l1+l2), 1.8*(l1+l2), '', 'FontSize', 12, 'FontWeight', 'bold');
    
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
        
        % Update the time annotation
        set(time_text, 'String', sprintf('Time: %.2f s', t(i)));
        
        drawnow;
        
        % Clear the figure except the ground line and time annotation for the next frame
        if i < length(t)
            clf;
            hold on;
            axis equal;
            axis([-2*(l1+l2) 2*(l1+l2) -2*(l1+l2) 2*(l1+l2)]);
            plot([-2*(l1+l2) 2*(l1+l2)], [0 0], 'k', 'LineWidth', 1); % ground line
            time_text = text(-1.5*(l1+l2), 1.8*(l1+l2), sprintf('Time: %.2f s', t(i)), 'FontSize', 12, 'FontWeight', 'bold');
        end
    end
end