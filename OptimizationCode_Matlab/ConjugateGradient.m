function [x,i] = ConjugateGradient(f,g,x0,nmax,tol,lstol)
    function meh = FindBeta(c_k,c_k_pre) % Using rules to determine direction
    y_k = c_k-c_k_pre;
    
    Beta_FR = sum(c_k*c_k.')/sum(c_k_pre*c_k_pre.');
    Beta_PR = sum(c_k*y_k.')/sum(c_k_pre*c_k_pre.');
    
    if 0 <= Beta_PR && Beta_PR <= Beta_FR
        meh = Beta_PR;
    elseif Beta_PR > Beta_FR
        meh = Beta_FR;
    else
        meh = 0;
    end
    end
n = length(x0); % Getting amount of variables
k = 0; % Counting iterations
d_k_pre = 0; % setting direction to 0
c_k_pre = 0; % Setting previous gradient as 0
c_k = 0; % Setting gradient as 0
x_pre = zeros(n);   % Initialization previous x-value as 0
x = x0; % Importing x-values
i = 0; %Counting iterations of cost function

while (abs(norm(g(x))) >= tol && abs(norm(x-x_pre)) >= tol) && k < nmax % Checking conditions
    c_k = g(x); % Getting current gradient
    if mod(k,n+1) == 0  % Remaineder of division based on the number of variables
        beta_k = 0;
    else
        beta_k = FindBeta(c_k,c_k_pre);
        i = i+3;
    end % Getting the beta value to determine the direction
    d_k = -c_k+beta_k*d_k_pre;  % Getting the new direction
    qq=@(a)(f(x+a*d_k)); % % Making step size function based on the current x-values and direction
    [alpha,w] = LineSearch(qq,1e-1,nmax,lstol); % Line Searching based on the direction
    i = i+w; % Accumulating cost function iterations
    x_pre = x; % Saving previous x-value data
    x = x+alpha*d_k; % Getting new x-values
    d_k_pre = d_k; % Saving previous direction
    c_k_pre = c_k; % Saving previous gradient
    k = k+1; % Iteration the iteration
end
end