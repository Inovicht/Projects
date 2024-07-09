function [G,Ag] = eig_assignL(A,B,C,eigs)
% This function designs the gain from left eigenstructure assignment using
% the state space matrices and the desired placed eigenvalues. It is
% required the model is observable and controllable

%% Defining matrix sizes
r = size(B,2);      % Input
m = size(C,1);      % Ouput
nn = size(A,1);     % 2*DOF
j = r;              % Placement Left
%% Initiating matrices
Omega = zeros(nn+m,m,j);
chi = zeros(nn,m,j);
P = zeros(m,m,j);
phi = zeros(nn,1,j);
Phi = zeros(nn,j);
Gamma = zeros(m,j);
%% Iterating the eigenvalues
for i = 1:j
    [~,~,NO] = svd([A.'-eye(length(A))*eigs(i),-C.']); % Extracting all right singular vectors
    Omega(:,:,i) = NO(:,nn+1:end); % Extracting nullspace
    % Omega(:,:,i) = null([A.'-eye(length(A))*eigs(i),-C.']); % Second method for computing null-space
    chi(:,:,i)   = Omega(1:nn,:,i); % Computing Chi
    [~,~,V] = svd(chi);
    alpha = V(:,1);  % Defining alpha constants from the first right singular vector
    P(:,:,i)     = Omega(nn+1:end,:,i);
    phi(:,:,i)   = chi(:,:,i)*alpha;
    Phi(:,i)     = phi(:,:,i);
    Gamma(:,i)   = P(:,:,i)*alpha;
end
%% Determining the Gain and closed loop input
Theta = B.'*Phi;
G_transpose = Gamma/Theta; % Gain
G = G_transpose.';
Ag = A-B*G*C;
end