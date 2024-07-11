clear all; close all; clc;
f1 = @(x) 1000*x(1)^2+100*x(2)^2+10*x(3)^2+x(4)^2+100*x(1)+10*x(2)+x(3)+0.1*x(4);
g1 = @(x) [2000*x(1)+100,200*x(2)+10,20*x(3)+1,2*x(4)+0.1];


clear; close all; clc
f=@(x)(x(1)^2+x(2)^2-x(1)*x(2));
g=@(x)[2*x(1)-x(2),2*x(2)-x(1)];
x0=[3,5];
[x,i]=ConjugateGradient(f,g,x0,500,1e-4,1e-4);
sprintf('x-value is: %f\n',x)
sprintf("Iteration of cost function is: %.f",i)
sprintf("Function value is: %.5f",f(x))

%%
clear; close all; clc
f=@(x)(3*x(1)^2+2*x(2)^2+2*x(1)+7);
x0=[2,1];
g=@(x) [6*x(1)+2,4*x(2)];
[x,i]=ConjugateGradient(f,g,x0,500,1e-4,1e-4);
sprintf('x-value is: %f\n',x)
sprintf("Iteration of cost function is: %.f",i)
sprintf("Function value is: %.5f",f(x))

%%
clear; close all; clc
f=@(x)(x(1)^2+x(2)^2-2*x(1)-2*x(2)+4);
g=@(x) [2*x(1)-2,2*x(2)-2];
x0=[45,32];
[x]=ConjugateGradient(f,g,x0,500,1e-4,1e-4)
sprintf('x-value is: %f\n',x)
sprintf("Iteration of cost function is: %.f",i)
sprintf("Function value is: %.5f",f(x))
%%
clear; close all; clc
f=@(x)((x(1)-2)^2+(x(2)-1)^2);
g=@(x) [2*x(1)-4,2*x(2)-2];
x0=[45,32];
[x]=ConjugateGradient(f,g,x0,500,1e-4,1e-4)
sprintf('x-value is: %f\n',x)
sprintf("Iteration of cost function is: %.f",i)
sprintf("Function value is: %.5f",f(x))
%% 
clear; close all; clc
f=@(x)((x(1)+x(2))^2+(x(2)+x(3))^2);
g=@(x) [2*x(1)+2*x(2),2*x(1)+4*x(2)+2*x(3),2*x(2)+2*x(3)];
x0=[45,32,0];
[x]=ConjugateGradient(f,g,x0,500,1e-4,1e-4)
sprintf('x-value is: %f\n',x)
sprintf("Iteration of cost function is: %.f",i)
sprintf("Function value is: %.5f",f(x))


% f = f;
% g = g;
% x0 = [2,2,2,2];
% nmax = 500;
% tol = 1e-4;
% lstol = 1e-3;
% [x,i] = ConjugatedGradient(f,g,x0,nmax,tol,lstol);
sprintf('x-value is: %f\n',x)
sprintf("Iteration of cost function is: %.f",i)
sprintf("Function value is: %.5f",f(x))