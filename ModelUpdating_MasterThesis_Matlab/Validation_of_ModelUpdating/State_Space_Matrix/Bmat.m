function [B,B2] = Bmat(M,b)
%Function that computes the input matrix B

%Inputs:
% Global mass matrix M
% vector b that contains the DoF that can undergo external excitation

%Outputs:
%input matrix B in state space
%Input matrix B2 in second order formulation

n = length(M);
r = length(b);
B2 = zeros(n,r);
k = 1;

for i = 1:length(b)
    B2(b(i),k) = 1;
    k = k + 1;
end

B = [zeros(n,r);inv(M)*B2];
end
