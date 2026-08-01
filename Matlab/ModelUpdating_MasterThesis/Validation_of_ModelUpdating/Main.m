% This is the main file where the code can be run. The parameters can be
% adjusted in the "params" file which contain most parameters.


%%  Importing functions
clear all; close all; clc;
addpath([cd,'/Covariance_identification/'])
addpath([cd,'/State_Space_Matrix/'])
addpath([cd,'/Optimization/'])
addpath([cd,'/Newton_Raphson/'])
addpath([cd,'/Subspace_Identification/'])
addpath([cd,'/Generalfunctions/'])
addpath([cd,'/EigStrucAssign/'])
format default

%% Parameters
rng(9)                                                              % Using same random generator
nodes = 20;                                                         % Defining amount of nodes for FEM model
nsr = 0.01;                                                         % Procent measurement noise for virtual data
t = 100;                                                             % Total time for virtual data
params = params(nodes);                                             % Defining params struc for all variables
theta = [2.806*1e8,11382,2.896*1e8,12783,3.23e4].';                 % Model parameters (trans,rot,trans,for,tension)
per = 0*theta;                                                      % Pertubated model parameters

%% Kinematic boundary conditions and system matrices
sys_nom = Sys(params,theta);                                        % Nominal system model
sys_tilt = Sys(params,theta-per);                                   % Pertubated system model
% Place eigenvalues
sys_nom.lambdapl = sort(eigen(params,sort(sys_nom.lambda)));        % Placed eigenvalues from nominal system model
%sys_nom.lambdapl = reshape(sys_nom.lambdapl.',1,[]).';             % Reshape may be required if using right eigenstructure assignemt
%% Retrieve Data
sample_rates = abs(diag(sys_tilt.A))/(2*pi);                        % Natural frequencies for pertubated system model
fs = 3*sample_rates(2*3);                                           % Chosen sample rate for generating virtual data
[uk,y,ns] = virtualdata(sys_tilt,fs,t,nsr);                         % Generated virtual data (Input, output, system order)
%% System Identification
sys_id= SysID(fs,uk,y,params);                                      % Generated state space model from subspace identification
res_id = sysID_cov(fs ,uk,y,params);                                % Generated state space model from covariance system identification
%% Eigenstructure Assignment
params.G = eigassign(sys_nom,params);                               % Computed gains from eigenstructure assignment
[params.lambdatiltideal] = findl(sys_nom.lambdapl,sys_tilt,params); % Eigenvalues from pertubated system model generated from subspace identification
[params.lambdatiltsysid] = findl(sys_nom.lambdapl,sys_id,params);   % Eigenvalues from pertubated system model generated from covariance system identification
%% Calculating the covariance
[params.cov,params.lambdatilt] = CovarianceMat(res_id,params,fs);   % Generated covariances and eigenvalues
%% Model update of stiffness
theta_new = optimization(params,theta);                             % Estimated model parameters from optimization problem
sys_new = Sys(params,theta_new);                                    % New state space matrix using estimated model parameters
lambda_new = findl(sys_nom.lambdapl,sys_new,params);                % Estimated eigenvalues using estimated model parameters
%% Compare confidensinterval with Monte Carlo Simulation
MS_simu = readmatrix("MS_Simu0.011.csv");                           % Importing Monte Carlo estimates
drawelipse(MS_simu,params,params.lambdatiltideal,lambda_new)        % Drawing uncertainty ellipse with Monte Carlo Estimates and estimated eigenvalues