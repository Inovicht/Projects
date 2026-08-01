function [x,i] = SteepestDescent(f,g,x0,xl,xu,nmax,tol,lstol)
n = length(x0);                                             % Getting amount of variables
k = 0;                                                      % Counting iterations
x_pre = zeros(n);                                           % Initialization previous x-value as 0
x = x0;                                                     % Importing x-values
i = 0;                                                      % Counting iterations of cost function

while (abs(norm(g(x))) >= tol && abs(norm(x-x_pre)) >= tol) && k < nmax % Checking conditions
    d_k = -g(x);                                            % Getting current gradient
    qq= @(a) (f(p(x+a*d_k)));                               % Making step size function based on the current x-values and direction
    [alpha,w] = LineSearch(qq,1e-1,nmax,lstol);             % Line Searching based on the direction with projection
    i = i+w;                                                % Accumulating cost function iterations
    x_pre = x;                                              % Saving previous x-value data
    x = p(x+alpha*d_k);                                     % Getting new x-values with projection
    k = k+1;                                                % Iteration the iteration
end
    function  x = p(x)                                      % Projection method
        for qw = 1:length(x)                                % Length of x
            if x(qw) < xl                                   % If under lower boundary
                x(qw) = xl;                                 % Setting lower boundary value
            elseif x(qw) > xu                               % If over upper boundary
                x(qw) = xu;                                 % Setting higher boundary value
            end
        end
    end
end