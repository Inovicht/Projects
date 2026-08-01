function [sysID] = SysID(fs,uk,y,params)
% This functions generates the state space model from subspace
% identification using input and output data.

% Inputs
    % fs
        % Contains the sample rate
    % uk
        % Contains the input data
    % y
        % Contains the output data

    % params
        % Contains the parameters
% Output
    % sysID
        % Containts the generated state space model from subspace
        % identification


ns = params.n_trunk; % Truncated model order
dt = 1/fs;           % Time step
dof = ns/2;          % Degree of freedom
Nsim = 1;            % Amount of simulations

for i=1:Nsim
    [AdID,BdID,CdID,DdID]=n4sid(uk,y,ns,params);  % Estimate state-space model using n4sid
    AID=logm(AdID)/dt;                            % Computing continuous A-matrix
    BID =AID/(AdID-eye(size(AdID)))*BdID;         % Computing continuous B-matrix
    CID = CdID;                                   % Computing continuous C-matrix
    lamIDr=eig(AID);                              % Eigvalues of state-space model
    [~,i2]=sort(imag(lamIDr));                    % Taking index of eigenvalues by sorting size
    lamIDrr=[lamIDr(i2(dof+1:end)) ; conj(lamIDr(i2(dof+1:end)))];
    lamID(:,i)=lamIDrr;                           % saving Eigenvalues 
end

% Compiling the matrices into a struc
sysID.A = AID;
sysID.B = BID;
sysID.Cc = CdID;
sysID.lambda = lamID;