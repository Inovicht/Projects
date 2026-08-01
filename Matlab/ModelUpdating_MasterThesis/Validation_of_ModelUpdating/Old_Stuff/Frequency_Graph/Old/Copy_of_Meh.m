clc; close all; clear all;
format short
Data = readmatrix("Good_Data_200.csv");
iteration = Data(end,1);
n = sqrt(iteration);
for i = 1:3
    omega{i} = Data(:,1+i);
    omega_temp = omega{i};
end
for i = 1:5
    theta{i} = Data(:,4+i);
end

% Rearrange Data for 3D plot
close all; clc
theta_q{1} = unique(theta{1});
theta_q{2} = unique(theta{2});
S = [numel(theta_q{1}),numel(theta_q{2})];
theta_surf{1} = reshape(theta{1},S);
theta_surf{2} = reshape(theta{2},S);
theta2d{1} = max(theta_surf{1},[],1);
theta2d{2} = max(theta_surf{2},[],2);
for i = 1:3
    omega_surf{i} = reshape(omega{i},S);
    omega2d_trans{i} = max(omega_surf{i},[],1);
    omega2d_rot{i} = max(omega_surf{i},[],2)
end
uf = [0 1e9];
%% Show 2D Plot
close all;
figure(1)
hold on
plot(theta2d{1},omega2d_rot{1},'LineWidth',3)
ylim([0 7000]);
xlim(uf)
set(gca,'xscale','log','FontName','cmr15')
title("\textbf{First natural frequency}",'FontName','Times New Roman','FontSize',12,'Interpreter','latex')
grid on
xlabel("Rotational Stiffness (Nm/rad)",'FontName','Times New Roman','FontSize',12,'Interpreter','latex')
ylabel("Natural frequency (Hz)",'FontName','Times New Roman','FontSize',12,'Interpreter','latex')
hold off
figure(2)
hold on
xlim(uf)
ylim([0 2e4]);
plot(theta2d{1},omega2d_rot{2},'LineWidth',3)
set(gca,'xscale','log','FontName','cmr15')
title("\textbf{Second natural frequency}",'FontName','Times New Roman','FontSize',12,'Interpreter','latex')
grid on
xlabel("Rotational Stiffness (Nm/rad)",'FontName','Times New Roman','FontSize',12,'Interpreter','latex')
ylabel("Natural frequency (Hz)",'FontName','Times New Roman','FontSize',12,'Interpreter','latex')
hold off
figure(3)
hold on
plot(theta2d{1},omega2d_trans{1},'LineWidth',3)
ylim([0 7000]);
xlim(uf)
set(gca,'xscale','log','FontName','cmr15')
title("\textbf{First natural frequency}",'FontName','Times New Roman','FontSize',12,'Interpreter','latex')
grid on
xlabel("Translational Stiffness (N/m)",'FontName','Times New Roman','FontSize',12,'Interpreter','latex')
ylabel("Natural frequency (Hz)",'FontName','Times New Roman','FontSize',12,'Interpreter','latex')
hold off
figure(4)
hold on
plot(theta2d{1},omega2d_trans{2},'LineWidth',3)
ylim([0 2e4]);
xlim(uf)
set(gca,'xscale','log','FontName','cmr15')
title("\textbf{Second natural frequency}",'FontName','Times New Roman','FontSize',12,'Interpreter','latex')
grid on
xlabel("Translational Stiffness (N/m)",'FontName','Times New Roman','FontSize',12,'Interpreter','latex')
ylabel("Natural frequency (Hz)",'FontName','Times New Roman','FontSize',12,'Interpreter','latex')
hold off
%% Show 3D Plot
figure(5)
x = 2;
y = 1;
z = 2;
surf(theta_surf{x},theta_surf{y},omega_surf{z},"LineStyle","-")
colormap([1 1 1])
set(gca,'xscale','log','YScale','log')
title("\textbf{Behavior of first natural frequency}",'FontName','Times New Roman','FontSize',12,'Interpreter','latex');
mehx = xlabel("Rotational Stiffness (Nm/rad)",'FontName','Times New Roman','FontSize',12,'Interpreter','latex');
%mehx.Position = mehx.Position + [30000 0 0];
mehy = ylabel("Translational Stiffness (N/m)",'FontName','Times New Roman','FontSize',12,'Interpreter','latex');
mehy.Position = mehy.Position - [-0 -200000 -0.14];
zlabel("Natural frequency (Hz)",'FontName','Times New Roman','FontSize',12,'Interpreter','latex')