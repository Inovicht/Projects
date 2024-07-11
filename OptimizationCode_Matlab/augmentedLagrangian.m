function [xk,vk]=augmentedLagrangian(x0,f,h,df,dh,pars)
mu = pars.mu0;  % Initial mu penalty paramter
x_pre = x0';    % Guessed start value
i = 0;          % Amount of iterations
epsilon = pars.epsilon; % Accuracy of 
k = 0;  % Amount of iterations
v_k = pars.v0;
if isempty(v_k)
    v_k = 0;
    h = @(x) [0];
    dh = @(x) zeros(length(x_pre),1);
end

while mu < pars.muMax && (k == 0 || norm(x_pre-x_pre_pre) > pars.tau)   % Checking the conditions
    lags = @(x) f(x) - v_k'*h(x) + sum((mu/2)*h(x).^2); % Accumulated Lagrangian
    dlag = @(x) df(x)' - (v_k-mu*h(x))'*dh(x)'; % Gradient of each variables
    [x,q] = ConjugateGradient(lags,dlag,x_pre,1e3,epsilon,1e-6); % Minimizing lagrangian with conjugated gradient and line search
    i = i+q;    % Accumulating the cost function
    v_k = v_k-mu*h(x);  % Updating my lagrangian multipliers
    mu = pars.factor*mu;    % Updating my penalty parameter
    epsilon = epsilon/pars.eps_factor; % Updatign accuracy
    x_pre_pre = x_pre;  % Saving x^i-2
    x_pre = x;  % Saving previous x values
    k = k+1    % counting iteration
end
xk = x'; % Returning final x-values
vk = v_k; % Returning Final lagrangian multiplier
end