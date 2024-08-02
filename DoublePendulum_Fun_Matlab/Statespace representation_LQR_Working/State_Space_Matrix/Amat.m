function [A,K,D,M] = Amat(params)
L = params.L;
m = params.m;
g = params.g;


M = [L(1)^2*m(1) + L(1)^2*m(2), L(1)*L(2)*m(2);
     L(1)*L(2)*m(2), L(2)^2*m(2)];
K = g * [L(1)*m(1) + L(1)*m(2), 0;
         0, L(2)*m(2)];
D = zeros(size(M));
A = [zeros(length(M)), eye(length(M));
     -M\K, zeros(length(M))];
end