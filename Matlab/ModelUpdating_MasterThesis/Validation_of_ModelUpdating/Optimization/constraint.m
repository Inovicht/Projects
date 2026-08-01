function [h,heq] = constraint(params,theta)
%This function computes the inequality constraints that are used for the optimization
%problem. In this case, the inequality are based on the confidence ellipse
%that are generated from the confidence intervals of of the target
%eigenvalues

%inputs:
%structure params (Check params.m for content and description)
%vector that contains the values of the boundary springs and axial tension
%[k1,k2,k3,k4,N]


sys_opd = Sys(params,theta);
p=0.99;                                                         %Confidence level
s=-2*log(1-p);
lambdatilt = params.lambdatilt;
h = zeros(length(lambdatilt),1);
lambda_temp = findl(lambdatilt,sys_opd,params);
for i = 1:length(lambdatilt)
    lambda_nom = [real(lambdatilt(i));imag(lambdatilt(i))];    %We augment the predicted eigenvalues into the real and imaginary part
    lambda_tilt =[real(lambda_temp(i));imag(lambda_temp(i))];  %We augment the target eigenvalues into the real and imaginary part
    diff = lambda_tilt-lambda_nom;                      
    h(i) = real(diff'/params.cov{i}*diff) - s;                 %inequality constraint the need to be lower or equal to zero
end
heq = [];