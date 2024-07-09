function [G,Ag] = eig_assignR(A,B,C,eigs)
% This function designs the gain from left eigenstructure assignment using
% the state space matrices and the desired placed eigenvalues. It is
% required the model is observable and controllable

%% Defining matrix sizes
r = size(B,2);      % Input
m = size(C,1);      % Ouput
nn = size(A,1);     % 2*DOF
j = length(eigs);   % Eigenvalues
%% Initiating matrices
Omega = zeros(nn+r,r,j);
chi = zeros(nn,r,j);
P = zeros(r,r,j);
psi = zeros(nn,1,j);
Psi = zeros(nn,m);
Gamma = zeros(r,m);
alpha = 1;
%% Iterating the eigenvalues
for i = 1:j
    [~,~,NO] = svd([A.'-eye(length(A))*eigs(i),-B]);
    Omega(:,:,i) = NO(:,nn+1:end);
    % Omega(:,:,i) = null([A-eye(length(A))*eigs(i),-B]);
    chi(:,:,i)   = Omega(1:nn,:,i);
    P(:,:,i)     = Omega(end-r+1:end,:,i);
    psi(:,:,i)   = chi(:,:,i)*alpha;
    Psi(:,i)     = psi(:,:,i);
    Gamma(:,i)   = P(:,:,i)*alpha;
end
%% Determining the Gain and closed loop input
Theta = C*Psi;
G = Gamma*inv(Theta); % Gain
Ag = A-B*G*C;
end