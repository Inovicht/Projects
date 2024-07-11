clear all; close all; clc;

f1 = @(x) (x-1)^2;
f2 = @(x) (x-1.75)^2;
f3 = @(x) (x-0.01)^2;
f4 = @(x) (x-5)^2-x;
f5 = @(x) -(x-5)^3+5*x;
%f1 = @(x)((x^2-0.001)^4)

f = f3
delta = 0.5;
kmax = 1e4;
tol = 1e-10;
[x,n] = LineSearch(f,delta,kmax,tol);
sprintf("x-value is: %f",x)
sprintf("Iteration of cost function is: %.f",n)
% sprintf("Function value is: %.5f",f(x)