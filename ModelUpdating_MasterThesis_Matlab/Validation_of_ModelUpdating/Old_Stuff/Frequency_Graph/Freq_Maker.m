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
%force = 0;
iteration = 100;

% Stiffnesses to check the frequency over
transit = logspace(-2,10,iteration);
rotit = logspace(-2,10,iteration);

% Defining size of storage matrix
storage = zeros(length(rotit)+length(transit),n);
size(storage);
%% Simulation
storage = zeros(length(rotit)+length(transit),n);
counter = 1;
for trans1 = transit
    for rot1 = rotit
            theta = [trans1,rot1,trans1,rot1,force];
            sys = Sys(params,theta);
            omega = sqrt(eigs(sys.K,sys.M,3,"smallestabs"))'; % Take the smallest 3 frequencies
            storage(i,:) = [i,omega,theta([1 2 5])];
            i = i+1;
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
Name = "Data_N" + string(nodes) + "It" + string(iteration) + "F" + string(force) +".csv";
writematrix(storage,Name,'Delimiter','tab');