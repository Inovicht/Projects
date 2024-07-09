function G = eigassign(sys,params)
% This function designs the gain based on left or right eigenstructure
% assignment. Recommended to use left eigenstructure assignment

% Inputs
    % sys
        % Contains the state space model the gain will be designed from
        % (Usually the nominal model)
    % params
        % Contains 

% Output
    % G
        % Contains the designed gains

eigs = sys.lambdapl; % Placed eigenvalues

% Computing the gains
for i = 1:size(eigs,1)
    if params.eiga == "r"
        [G(i,:),~] = eig_assignR(sys.A,sys.B,sys.Cc,eigs(i,:)); % Right eigenstructure assignment
    else
        [G(i,:),~] = eig_assignL(sys.A,sys.B,sys.Cc,eigs(i,:)); % Left eigenstructure assignment
    end
end

% Resizing the gain matrix
if params.eiga == "r"
    G(1,:) = zeros(1,size(G,2));
else
    for ol = 1:params.n_CL+1:(1+params.n_CL)*params.n_OL
        G(ol,:) = zeros(1,size(G,2));
    end
end