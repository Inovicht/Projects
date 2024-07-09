function params = params(nodes)
% This function contains the parameters used in the code

%% Material, geometry, and topology
params.n_trunk = 2*3;               % Amount of nodes the model is truncated to
noElem = nodes-1;                   % Amount of elements
params.E =     2.07e11;             % E-modulus
D = (2*5.9)/1000;                   % Diameter
params.I =     pi*D^4/64;           % Area moment of Inertia
params.A =     (D/2)^2*pi;          % Area               
A = params.A;
rho =   7850;                       % Uniform density
x_in =  0;                          % Start point of beam
x_fi =  0.225;                      % End point of beam (Length)
x_b =   [x_in:x_fi/noElem:x_fi]';   % X-coordinates for nodes
y_b =   zeros(noElem+1,1);          % Y-coordinates for nodes
params.coor =  [x_b y_b];           % Coordinates for nodes
params.no   =  [[1:numel(x_b)-1]' [2:numel(x_b)]' params.I*ones(numel(x_b)-1,1)...
params.E*ones(numel(x_b)-1,1) A*ones(numel(x_b)-1,1) rho*ones(numel(x_b)-1,1)];
params.fixeddof=linspace(1,(1+noElem)*3-2,1+noElem);
%% Bolt tension estimation - Marie Brøns PhD
pro = 0.58;             % Extracted data from Marie Brøns PhD
sigma = 640*1e6*pro;
params.N = sigma*A;     % Computed bolt tension
%% EigenStructure Assignment
params.n_CL = 2;        % Number of closed-loop eigenvalue for each open-loop eigenvalue
params.n_OL = 3;        % Number of open-loop eigenvalue that closed-loop eigenvalues will be assigned
params.eiga = "l";      % Determining to use left "l" or "r" eigenstructure assignment
%% Damping parameters
params.zeta_nom = 4.75e-4*ones(2*length(x_b),1); % Computed damping ratio from Marie Brøns
params.zeta_tilt = params.zeta_nom;              % Pertubated state space matrix damping ratio
%% State space parameters
% Mapping the output vector (Which node/s)
params.d       =   [nodes*2-1,nodes*2-3];   % Displacement output
params.ddot    =   [];                      % Velocity output
params.dddot   =   [];                      % Acceleration output
% Mapping the input vector
params.q = [1];                             % Force input

%% Parameters
params.sysidNhankel = 12;                   % Parameter for hankel matrix in subspace identification
params.sysCoNhankel = 10;                   % Parameter for covariance hankel blocks in covariance system identification
params.sysCoNcocompute = 100;              % Number of blocks for covariance computation in covariance system identification
params.SysCons = params.n_trunk;            % Model order 
%% Configuration of optimization problem:
params.Ab = [];
params.b = [];
params.Aeq = [];
params.beq = [];
params.lb = [1e3,1e2,1e3,1e2,0];            % Lower limit for optimization problem for model parameters
params.ub = [4e8,8e5,4e8,8e5,params.N*1.1]; % Upper limit for optimization problem for model parameters
end