function [x,i] = LineSearch(f,delta,kmax,tol)
i = 0; % Amount of cost functions
% f is the cost function
% delta is the initial step size which can be corrected
% kmax is the total amount of iterations
% tol is the tolerance
    function [step,i] = Phase0(func,i,step,kmax)
    k = 0; % Iterating until kmax
    while func(step)>func(0) && k < kmax
        i = i+2;
        step = step/2; % Decreasing step if step too high
        k = k+1;
    end
    end
    function [alpha,i] = Phase1(func,i,kmax,delta)
    a_1 = 0; % First check
    GR = 1.618; % Golden ratio
    for n = 1:kmax
        a_0 = a_1;
        a_1 = a_1+delta*GR^n;
        i = i+2; % Iterating through until finding range
        if func(a_0) <= func(a_1)
            a_0 = a_0-delta*GR^(n-1);
            alpha = [a_0,a_1]; % Getting range
            break
        end
    end
    if ~exist('alpha', 'var') % get error if something wrong
      error('Kmax is hit. You have found an unbounded minimum. Change delta');
    end
    end
    function [alpha,i] = Phase2(func,alpha,i,tol,kmax)
    k = 0;
    tau = 0.618;
    a_l = alpha(1); % Lower range
    a_u = alpha(2); % Higher range
    I = a_u-a_l;
    while I >= tol && k < kmax % See if tolerance is hit
        a_a = a_l+(1-tau)*I;
        a_b = a_l+tau*I;
        aa = func(a_a);
        ab = func(a_b);
        if aa < ab % Checking each range
            a_u = a_b;
        elseif aa > ab % Checking each range
            a_l = a_a;
        elseif aa == ab % Checking each range
            a_l = a_a;
            a_u = a_b;
        end
        I = a_u-a_l; % Looking at the range
        i = i+2;
        k = k+1;
    end
    alpha = [a_l a_u]; % Range
    end
    function [alpha,i] = Interpol(func,alpha,i)
    a_u = alpha(1); % Interpolating one time between the range
    a_l = alpha(2);
    a_i = mean(alpha); % Mean
    fau = func(a_u); % Cost function values
    fal = func(a_l);
    fai = func(a_i);
    
    a_2 = (1/(a_u-a_i))*((fau-fal)/(a_u-a_l)-(fai-fal)/(a_i-a_l));
    a_1 = (fai-fal)/(a_i-a_l)-a_2*(a_l+a_i);
    a_0 = fal-a_1*a_l-a_2*a_l^2;
    i = i+3;
    alpha = -a_1/(2*a_2);
    if isnan(alpha) % Does alpha exits?
        alpha = a_i; % Finding the final value
    end
    end
% Iterating through the functions
[step,i] = Phase0(f,i,delta,kmax);
[alpha,i] = Phase1(f,i,kmax,step);
[alpha,i] = Phase2(f,alpha,i,tol,kmax);
[x,i] = Interpol(f,alpha,i);
% outputting final x step size
% i is the final amount of use of cost function
end