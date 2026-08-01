function [lambda_tilt] = findl(lambda,sys_tilt,params)
% This function finds the eigenvalues from the pertubated state space model 
% "sys_tilt" that have the least euclidean distance from the eigenvalues from the
% nominal state space model defined in "lambda"

% Input
    % lambda
        % Contains the nomial eigenvalues
    % sys_tilt
        % Contains pertubated state space model
    % params
        % Contains "eiga" that determine which method of eigenstructure
        % assignment is used
        % The designed gain "G"
% Output
    % lambda_tilt
        % Contains the eigenvalues from pertubated state space model 
        % that have minimal euclidean distance from the eigenvalues in "lambda"

%% Using right eigenstructure assignment
if params.eiga == "r"
    n = size(params.G,2);
    m = size(params.G,1);
    lambda_tilt = zeros(n*m,1);
    for i = 1:m
        Ag_tilt = sys_tilt.A - sys_tilt.B*params.G(i,:)*sys_tilt.Cc;    % Pertubated A-matrix
        l = eig(Ag_tilt);                                               % Pertubated eigenvalues
        indx = indxfind(l,lambda(i*n-n+1:i*n));                         % Getting index of pertubated eigenvalues with least distance
        lambda_tilt(n*i-n+1:n*i) = l(indx);                             % Extracting the eigenvalues
    end
else 
%% Using left eigenstructure assignment
    lambda_tilt = zeros(length(lambda),1);
    for i = 1:length(lambda)
        Ag_tilt = sys_tilt.A - sys_tilt.B*params.G(i,:)*sys_tilt.Cc;    % Pertubated A-matrix
        l = eig(Ag_tilt);                                               % Pertubated eigenvalues
        indx = indxfind(l,lambda(i));                                   % Getting index of pertubated eigenvalues with least distance
        lambda_tilt(i) = l(indx);                                       % Extracting the eigenvalues
    end
end
