function params = params(nodes)
%% Material, geometry, and topology
params.n_trunk = 6;
noElem = nodes-1;
params.E =     2.07e11;                % E-modulus
D = (2*5.9)/1000;
params.I =     pi*D^4/64;         % Area moment of Inertia
A =     (D/2)^2*pi;                     % Area
params.A = A;
rho =   7850;                       % Densit
x_in =  0;                          % Start point of beam
x_fi =  0.225;                      % End point of beam
x_b =   [x_in:x_fi/noElem:x_fi]';   % X-coordinates for nodes
y_b =   zeros(noElem+1,1);          % Y-coordinates for nodes
params.coor =  [x_b y_b];                  % Coordinates for nodes
params.no   =  [[1:numel(x_b)-1]' [2:numel(x_b)]' params.I*ones(numel(x_b)-1,1)...
params.E*ones(numel(x_b)-1,1) A*ones(numel(x_b)-1,1) rho*ones(numel(x_b)-1,1)];
params.fixeddof=linspace(1,(1+noElem)*3-2,1+noElem);
%% Data
pro = 0.76;
sigma = 640*1e6*pro;
params.N = sigma*A;
%% EigenStructure Assignment
params.n_CL = 2;
params.n_OL = int32(params.n_trunk/2);
params.n_OL = 3;
params.eiga = "l";
%% Damping parameters
params.zeta_nom = 4.75e-4*ones(2*length(x_b),1);
%params.zeta_nom = 0*ones(2*length(x_b),1);
params.zeta_tilt = params.zeta_nom;
%% State space parameters
% Mapping the output vector
params.d       =   [];
params.ddot    =   [];
params.dddot   =   [1];
% Mapping the input vector
params.q = [1];

%% Parameters
params.sysidNhankel = 38*2;
params.sysCoNhankel = 70;
params.sysCoNcocompute = 200; % Number of blocks for covariance computation
params.SysCons = 6;
%% Configuration of optization problem:
params.Ab = [];
params.b = [];
params.Aeq = [];
params.beq = [];
params.lb = [1e3,1e2,1e3,1e2,1e3];
params.ub = [1e9,1e6,1e9,1e6,1e8];
end

