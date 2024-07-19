function w = NewtonRaphston(t,q)
es = 0.000000000001;     % Tolerance
maxit = 100000;       % Maximum iteration
iter = 0; % Defining iteration
driving = zeros(1,6);
driving(3) = 1;
while(1)
    Phi1 = [Phi(params,t,q);q(3)*1*t];
    Jac = [Phi_q(q,params); driving];
    dq = Jac\Phi1; % Defining the change
    q = q-dq; % Removing change from previous guess
    iter = iter+1; % Counting Iteration
    ea = 100*max(abs(dq./q)); % Looking at error-value
    if(iter>=maxit || ea<=es)  % See if conditions are satisfied
        w = q; % Saving data if satisfied
        break
    end
end