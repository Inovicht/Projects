function [w,w_dot,w_ddot] = NewtonRaphston(t,q)
N = length(t);      % Length of time
es = 0.00000001;     % Tolerance
maxit = 10000;       % Maximum iteration


w = zeros(length(q),N); % Saving positions based on time-step
w_dot = zeros(length(q),N);         % Saving velocities based on time-step
w_ddot = zeros(length(q),N);        % Saving accelerations based on time-step
for i = 1:N
    % Newton Raphton for each time-step (Position and orientation)
    iter = 0; % Defining iteration
    % Newton-Raphton for specific time-step
    while(1)
        dq = Phi_q(q)\Phi(t(i),q); % Defining the change
        q = q-dq; % Removing change from previous guess
        iter = iter+1; % Counting Iteration
        ea = 100*max(abs(dq./q)); % Looking at error-value
        if(iter>=maxit || ea<=es)  % See if conditions are satisfied
            w(:,i) = q; % Saving data if satisfied
            break
        end
    end
    w_dot(:,i) = Phi_q(q)\Nu;                % Defining Velocity
    w_ddot(:,i) = Phi_q(q)\Gamma(q,w_dot);   % Defining Acceleration
end
