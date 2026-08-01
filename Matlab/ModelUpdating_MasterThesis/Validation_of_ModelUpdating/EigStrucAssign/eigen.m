function CL_eigs = eigen(params,lambda_OL)
% This function places the eigenvalues from the nominal state space model

% Inputs
    % Params
        % Contain the number of desired open-loop eigenvalues that
        % closed-loop eigenvalues will placed from
        % Contain the number of desired closed-loop eigenvalues for each
        % open-loop eigenvalue
% Output
    % CL_eigs
        % Contains all the eigenvalues 

% There is some erorr when having multiple inputs
if length(params.q) > 1
    error("Code does not work for multiple inputs");
end

n_OL = params.n_OL; % Amount of open-loop eigenvalues
n = params.n_CL;    % Amount of desired placed closed-loop eigenvalue for each open-loop eigenvalue
eigen = zeros(n,n_OL);
lambda_OL = sort(lambda_OL); % Sorted open-loop eigenvalues

% Extracting non-conjugated eigenvalues
for m = 1:n_OL
    eigen(1,m) = lambda_OL(2*m);
end

% Placing closed-loop eigenvalues
for i = 1:n         % Rows is amount of closed loops placements
    for j = 1:n_OL  % Columns is the amount of open loop eigenvalues to place
        eigen(i+1,j) = (1+1e-1*i)*lambda_OL(2*j); % Method to place closed-loop eigenvalue (This can be changed)
    end
end

% Reshaping "CL_eigs" for later use
if params.eiga == "r"
    CL_eigs = eigen;
else
    CL_eigs = reshape(eigen,1,[]).';
end
end