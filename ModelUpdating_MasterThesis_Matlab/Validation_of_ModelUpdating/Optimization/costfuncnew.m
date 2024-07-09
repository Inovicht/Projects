function [f,g] = costfuncnew(params,theta)

 %% Outputting parameters
sys_nom = Sys(params,theta);
lambda_tilt = params.lambdatilt;

%% Initializing Newton Raphson
lambda = findl(lambda_tilt,sys_nom,params);

%% Getting the gradient
g = Grad(params,sys_nom,params.G,theta,0.99).';

%% Updating the variables
lambda_choose = length(lambda);
lambda_tilttemp = lambda_tilt(1:lambda_choose,:);
lambda_temp = lambda(1:lambda_choose,:);
diff = lambda_tilttemp-lambda_temp;
f = abs(sum(diff.*diff));
end
