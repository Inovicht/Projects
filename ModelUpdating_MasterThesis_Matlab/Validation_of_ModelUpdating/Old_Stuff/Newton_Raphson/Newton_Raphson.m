function [theta]= Newton_Raphson(params,theta,delta)
    %% Outputting parameters
    C_c = params.C_c;
    B = params.B;
    G = params.G;
    lambda_tilt = params.lambda_tilt;
    %% Newton Raphson parameters
    es = 1e-6;
    maxit=10000;
    iteration = 0;
    %% Sizes of for loop
    m = size(G,1);   % Amount of inputs
    n = size(C_c,1); % Amount of outputs
    step = 1e-3;
    %% Initializing Newton Raphson
    for p = 1:maxit
        [A_nom,~,~] = Amat(params,theta);
        lambda = zeros(2*m,1);
        J = zeros(m*2,length(theta));
        for i = 1:m
            Ag = A_nom + B*G(i,:)*C_c;
            Ji = Jacw(params,G(i,:),theta,delta);
            l = eig(Ag);
            indx = indxfind(l,lambda_tilt,i,n);
            J(n*i-n+1:n*i,:) = Ji(indx,:);                 % Saving the eigenvalue jacobian
            lambda(n*i-n+1:n*i) = l(indx);                 % Saving the eigenvalue
        end
    %% Updating the variables
    Allow_eigs = params.Allow_eigs;
    Jtemp = J(1:Allow_eigs,:);
    lambda_tilttemp = lambda_tilt(1:Allow_eigs,:);
    lambda_temp = lambda(1:Allow_eigs,:);
    Jtemp = [real(Jtemp);imag(Jtemp)];
    dtheta = ((Jtemp.')*Jtemp)\(Jtemp.')*[real(lambda_tilttemp)-real(lambda_temp);imag(lambda_tilttemp)-imag(lambda_temp)];
    theta = theta+dtheta*step;
    theta = theta;
    iteration = iteration + 1;
    ea = max(abs(dtheta./theta));

    %% Checking the rank of the jacobian
    if rank(J) < length(theta) || any(dtheta > 1e8)
         error('Jacobian rank is reduced or dtheta is too big')
    end
    %% Checking if the criterea is satisfied
     if ea <= es
        break
     end
    if mod(p,10)==1
        fprintf('Error: %.10f Iteration: %d Stepsize: %.6f\n',ea,iteration,step)
        theta
    end
   end
end