clc; close all; clear all;
format short
addpath([cd,'/Scripts/'])
Data_p0 = readmatrix("Data_N20It500F0K4INF.csv");
Data_p1 = readmatrix("Data_N20It500F1000000K4INF.csv");
iteration = Data_p0(end,1);
n = sqrt(iteration);
for i = 1:3
    omega_p0{i} = Data_p0(:,1+i);
    omega_p1{i} = Data_p1(:,1+i);
end
for i = 1:3
    theta_p0{i} = Data_p0(:,4+i);
    theta_p1{i} = Data_p1(:,4+i);
end
force = Data_p1(:,end);
% Rearrange Data for 3D plot
S = [sqrt(iteration),sqrt(iteration)];
for i = 1:2
    theta_p1_surf{i} = reshape(theta_p1{i},S);
end
force_reshaped = reshape(force(1:length(theta_p0{1})),S);
for i = 1:3
    omega_p0_surf{i} = reshape(omega_p0{i},S);
    omega_p1_surf{i} = reshape(omega_p1{i},S);
end
omegak2 = omega_p1_surf{1}./omega_p0_surf{1};
%% Show 3D Plot
close all
figure(1)
hold on
contourf(force_reshaped,theta_p1_surf{1},omegak2,20)
set(gca,'YScale','log','GridLineWidthMode','auto')
colormap(flipud(hot))
c=colorbar;
title("\textbf{Contour of the first normalized natural frequency (rad/s)}",'FontName','Times New Roman','FontSize',12,'Interpreter','latex');
mehx = xlabel("Bolt Tension (N)",'FontName','Times New Roman','FontSize',12,'Interpreter','latex');
mehy = ylabel("Translational Stiffness - $K_3$ (N/m)",'FontName','Times New Roman','FontSize',12,'Interpreter','latex');
zlabel("Natural frequency (rad/s)",'FontName','Times New Roman','FontSize',12,'Interpreter','latex')
hold off
figure(2)
hold on
contourf(force_reshaped,theta_p1_surf{2},omegak2,25)
set(gca,'YScale','log')
colormap(flipud(hot))
c=colorbar;
title("\textbf{Contour of the first normalized natural frequency (rad/s)}",'FontName','Times New Roman','FontSize',12,'Interpreter','latex');
mehx = xlabel("Bolt Tension (N)",'FontName','Times New Roman','FontSize',12,'Interpreter','latex');
mehy = ylabel("Rotational Stiffness - $K_4$ (Nm/rad)",'FontName','Times New Roman','FontSize',12,'Interpreter','latex');
zlabel("Natural frequency (rad/s)",'FontName','Times New Roman','FontSize',12,'Interpreter','latex')
hold off