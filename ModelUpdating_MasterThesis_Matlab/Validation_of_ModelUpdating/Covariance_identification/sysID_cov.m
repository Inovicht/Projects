function res_id = sysID_cov(fs,uk,y,params)
% This function generate state space matrix from covariance system
% identification

% Inputs
    % fs
        % Sample rate
    % uk
        % Input data
    % y
        % Output data
    % params
        % Parameters defined
        % sysCoNhanek
        % sysCoNcocompute
        % SysCons
% Output
    % res_id
        % Struc of state space matrices


%% Defining parameters
ns = params.SysCons;    % Model order
dof = ns/2;             % Degree of freedom
dt = 1/fs;              % Time Step
r = size(uk,1);         % Amount of input data
m = size(y,1);          % Amount of output data

yUukU=[y.' uk.'];

%% Local param struc used for covariance system identification
param.mthd='IOcov_all'; % I/O covariance-driven identification
param.pchannel=1:m;  % Projection channels, in case of many outputs
param.p=params.sysCoNhankel;    % Number of block columns in Hankel matrix
param.q=param.p+1;
param.nb = params.sysCoNcocompute; % Number of blocks for covariance computation
param.nmax=2*dof; % Model order
param.inputs=r; 
param.r=m;      
param.dt=dt;
%% Covariance system identification
res_id=system_id(yUukU,param); % Output A,B,C,D,Delta A,Delta B,Delta C,Delta D