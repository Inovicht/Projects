%%  FE formulation with "in-line" Euler-Bernoulli beam elements
clear all; close all; clc;
addpath([cd,'/Covariance_identification/'])
addpath([cd,'/EigStrucAssign/'])
addpath([cd,'/fmin/'])
addpath([cd,'/Generalfunctions/'])
addpath([cd,'/Newton_Raphson/'])
addpath([cd,'/State_Space_Matrix/'])
addpath([cd,'/Statematrices/'])
addpath([cd,'/Subspace_Identification/'])
format short
nodes = 20;
max_it = 10;
%% Parameters
%theta = [10e6,1000,10e6,1000,26*1e4]';    % ltrans,rtrans,lrot, rrot
theta = [2.806*1e8,11382,2.896*1e8,12783,3.23e4].';
per = 0*[-160,-105,-400,-560,-180]';           % Pertubation
% System identication
nsr = 0.01;
t = 100;
Nsim = 1;
%% Kinematic boundary conditions and system matrices
params = params(nodes);
sys_nom = Sys(params,theta);
sys_tilt = Sys(params,theta+per);

% Place eigenvalues
sys_nom.lambdapl = eigen(params,sys_nom.lambda);
%sys_nom.lambdapl = reshape(sys_nom.lambdapl.',1,[]).';
%% State-Space Identification
sample_rates = sqrt(eig(sys_nom.K,sys_nom.M))/(2*pi);
fs = 3*sample_rates(3);
% Making filename
z = 1;
% Check if file exits
while(1)
    if isfile('MS_Simu' + string(nsr) + string(z) + '.csv')
    z = z+1;
    else
        filename = 'MS_Simu' + string(nsr) + string(z) + '.csv';
        break;
    end
end
% Simulation of Monte Carlo
% n_dof = length(sys_nom.A)/2*size(sys_nom.Cc,1);
n_dof = length(sys_nom.lambdapl);
params.G = eigassign(sys_nom,params);
lambdatilt = zeros(max_it,2*n_dof);
for i = 1:max_it
    [uk,y,ns] = virtualdata(sys_tilt,fs,t,nsr);
    n = 1:size(y,2);
    sys_id= SysID(fs,uk,y,params);
    lambda_NO_MI = findl(sys_nom.lambdapl,sys_id,params);
    % sys_id_OD= SysID(3*ns,fs,uk,y);                        % We conduct sysID with higher model order
    % sys_id = TruncateSysID(sys_nom.lambda,sys_id_OD);      % We now truncate the higher order SysID to order ns
    lambda_temp = findl(sys_nom.lambdapl,sys_id,params);
    [abs(sys_nom.lambdapl-lambda_NO_MI),abs(sys_nom.lambdapl-lambda_temp)]
    for j = 1:n_dof
        lambdatilt(i,1:n_dof) = real(lambda_temp);
        lambdatilt(i,n_dof+1:end) = imag(lambda_temp);
    end
    i
end
%%
writematrix(lambdatilt,filename,'Delimiter','tab','WriteMode','append');


