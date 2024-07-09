function [f,lambda] = costfunc(params,theta)
%This function computes the cost function that is used for the optimization
%problem

%inputs:
%structure params (Check params.m for content and description)
%vector that contains the values of the boundary springs and axial tension
%[k1,k2,k3,k4,N]

%Outputs: 
%f: cost function evaluation
%lambda: Model prediction of eigenvalues

 %% Outputting parameters
sys_nom = Sys(params,theta);
lambda_tilt = params.lambdatilt;

%% Initializing Newton Raphson
lambda = findl(lambda_tilt,sys_nom,params);
%% Updating the variables
lambda_choose = length(lambda);
lambda_tilttemp = lambda_tilt(1:lambda_choose,:);
lambda_temp = lambda(1:lambda_choose,:);
f = abs(sum((lambda_tilttemp-lambda_temp).*(lambda_tilttemp-lambda_temp)));

