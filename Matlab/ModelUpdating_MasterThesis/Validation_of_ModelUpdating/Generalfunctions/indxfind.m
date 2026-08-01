function indx = indxfind(l,lambda)
% This function finds the index of the eigenvalues in "lambda" that have the
% least euclidean distance from eigenvalues in "l". If the same index is
% find multiple times, another index will be found

% Input
    % l
        % Contains all the eigenvalues
    % lambda
        % Contains the eigenvalues that we want to find in "l"
% Output
    % indx
        % Contains the indexes of the eigenvalues in "lambda" that are
        % closed to the eigenvalues in "l"

n = length(lambda); % Amount of eigenvalues we want to find in "l"
indx = zeros(n,1);
for j = 1:n
    [~,indx(j)] = min(abs(l-lambda(j))); % computing the index of the eigenvalues in "lambda" that are closest to a eigenvalue in "l"
    if sum(indx == indx(j)) > 1 % If the eigenvalue in "l" is found multiple times, another will be generated
        l_temp = l; % Computing a temporary array of the eigenvalues in "l"
        l_temp(indx(j)) = inf; % Setting the eigenvalue found multiple times to infinite
        [~,indx(j)] = min(abs(l_temp-lambda(j))); % Computing the index of the eigenvalues in "lambda" that are closest to a eigenvalue in "l" again
    end
end