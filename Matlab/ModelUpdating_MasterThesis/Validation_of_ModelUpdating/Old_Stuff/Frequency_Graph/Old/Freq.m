%%  FE formulation with "in-line" Euler-Bernoulli beam elements
clear all; close all; clc;
format short
%% Parameters
nodes = 25;
params = params(nodes);

%% Extracting parameteres from beam without boundaries
[K_beam,~,M,~]= beam2D_freg(params);
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
n = 1+5+3; % Length of storage
i = 1;
force = params.N;
iteration = 200;

% Stiffnesses to check the frequency over
transit = logspace(2,10,iteration);
rotit = logspace(1,7,iteration);

% Defining size of storage matrix
storage = zeros(length(rotit)+length(transit),n);
size(storage);
%% Simulation
storage = zeros(length(rotit)+length(transit),n);
counter = 1;
for trans1 = transit
    for rot1 = rotit
            theta = [trans1,rot1,trans1,rot1,force];
            K = K_Matrix_freg(K_beam,theta);
            omega = sqrt(eigs(K,M,3,"smallestabs"))'; % Take the smallest 3 frequencies
            storage(i,:) = [i,omega,theta];
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
writematrix(storage,"Good_Data_200.csv",'Delimiter','tab');