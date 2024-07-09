%%  FE formulation with "in-line" Euler-Bernoulli beam elements
clear all; close all; clc;
format short
addpath([cd,'/Scripts/'])
%% Parameters
nodes = 20;
params = params(nodes);
%% Making filename
z = 1;
% Check if file exits
while(1)
    if isfile('test' + string(z) + '.csv')
    z = z+1;
    else
        filename = 'test' + string(z) + '.csv';
        break;
    end
end
%% Initiation of simulation
% Iterate
n = 1+5-2+3; % Length of storage
i = 1;
force = params.N;
iteration = 500 ;

% Stiffnesses to check the frequency over
transit = logspace(-2,10,iteration);
%transit = 1e12;
rotit = logspace(-2,10,iteration);
rotit = 1e15;
forceit = linspace(0,1e6,iteration);
% Defining size of storage matrix
storage = zeros(length(rotit)+length(transit),n);
size(storage);
%% Simulation
storage = zeros(length(rotit)+length(transit),n);
counter = 0;
for trans1 = transit
    for rot1 = rotit
        for force = forceit
            theta = [1e15,trans1,1e15,1e15,force];
            sys = Sys(params,theta);
            omega = sqrt(eigs(sys.K,sys.M,3,"smallestabs"))'; % Take the smallest 3 frequencies
            storage(i,:) = [i,omega,theta([1 2 5])];
            i = i+1;
        end
      counter = counter+1 % See counter size of iteration to see how many
    end
    counter = counter+1 % See counter size of iteration to see how many
end
%% Plotting
clc; close all
figure(1)
scatter3(storage(:,6),storage(:,5),storage(:,3))
set(gca,'xscale','log','YScale','log')
title("Normalized frequency 1 over stiffnesses")
xlabel("Rotational Stiffness")
ylabel("Translational Stiffness")
zlabel("Frequency")
%% Save the data
Name = "Data_N" + string(nodes) + "It" + string(iteration) + "F" + string(force) + "K2INF" + ".csv";
writematrix(storage,Name,'Delimiter','tab');