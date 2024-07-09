function [uk,y,ns] = virtualdata(sys,fs,t,nsr)
% This function generates virtual data for system identification

% Inputs
    % sys
        % Contains the continuous state space matrix
    % fs
        % Contains the sample rate
    % t
        % Containts the total simulation time for virtual data
    % nsr
        % Contains the system order
% Outputs
    % uk
        % Contains the virtual data input
    % y
        % Contains the virtual data output
    % ns
        % Contains the model order
%% Extracting the state space matrices
M = sys.M;                          % Mass matrix
B = sys.B;                          % B-matrix
C_c = sys.Cc;                       % C-matrix
Ac = sys.A;                         % Continuous A-matrix
H = sys.H;                          % D-matrix in state matrix

%% Defining variables
m = size(C_c,1);                    % Amount of outputs
r = size(B,2);                      % Amount of inputs
dt = 1/fs;                          % Step size
Nsamp = int32(t/dt);                % Amount of total samples
dof = length(Ac)/2;                 % Amount of degree of freedom
t = linspace(0,t,Nsamp);            % Array of the time
zold=zeros(2*dof,1);                % Accumulated state response

% Discretizing the state space matrix
Ad=expm(Ac*dt);                     % Discrete state space
Bu=inv(Ac)*(Ad-eye(2*dof))*B;       % Discrete input-matrix
ukk=randn(r,Nsamp);                 % Array of random samples by normal distribution
sig_in = std(ukk,0,2);              % Standard deviation of input

yk = zeros(size(C_c,1),Nsamp);
for i=1:Nsamp
    zn=Ad*zold+Bu*ukk(:,i);          % Defining new output by random input
    yk(:,i)=C_c*zold-H*ukk(:,i);     % Saving new output
    zold=zn;                         % Saving accumulated state reponse
end
std_y=std(yk');                      % Standard deviation of output
ns=[2*dof];                          % Model order
noise_in=randn(r,Nsamp).*sig_in*nsr; % Generating the noise input
uk=ukk+noise_in;                     % Adding the noise on the input
y=real(yk+nsr*std_y'.*randn(m,Nsamp)); % Adding the noise on the output

% There are imaginary numbers on the output. This should not be there! This
% is probably cause by the large time step. To solve it, decrease the time
% step
end