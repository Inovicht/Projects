function [cov,lambdatilt] = CovarianceMat(res_id,params,fs)
%This function computes the confidence interval of each target eigenvalue

%Input: 
%res_id: Structure that contain the system quadruplet in discrete time
%domain and jacobians
%params structure that contains the augmented gain matrix
%Sample frequency fs

%Outputs:
%cell cov that cointains covarians matrix for each target eigenvalue
%vector lambdatilt that contain the target eigenvalues

if params.eiga == "r"                               %If we are using left or right eigenstructure assignment
    % Sizes of for loop
    m = size(params.G,1);   % Amount of inputs
    n = size(params.G,2); % Amount of outputs
    %% Finding confidence matrix
    lambdatilt = zeros(m*n,1);   
    for i = 1:m
        [cov_dAd] = Confidensinterval(res_id,params.G(i,:),fs);     %Computes a structure than contain a vector of closed-loop eigenvalues and a cell that contains confidence matrix for each closed-loop eigenvalues
        l = cov_dAd.e_vald;                                         %Closed-loop eigenvalues
        indx = indxfind(l,params.lambdatiltideal(n*i-n+1:n*i));     %We find the index numbers of the assigned eigevalue in the vector l
        place =n*i-n+1:n*i; 
        cov{place(1)} = cov_dAd.cov_lambda_Ab{indx(1)};             %The corresponding covariance matrices to each target eigenvalue are inserted in the cell cov
        cov{place(2)} = cov_dAd.cov_lambda_Ab{indx(2)};
        cov{place(3)} = cov_dAd.cov_lambda_Ab{indx(3)};
        lambdatilt(n*i-n+1:n*i) = l(indx);
    end
else
    lambdatilt = zeros(length(params.lambdatiltideal),1);
    for i = 1:length(params.lambdatiltideal)
        [cov_dAd] = Confidensinterval(res_id,params.G(i,:),fs);
        l = cov_dAd.e_vald;
        indx = indxfind(l,params.lambdatiltideal(i));
        cov{i} = cov_dAd.cov_lambda_Ab{indx};
        lambdatilt(i) = l(indx);
    end
end
