close all; clear all
% Defining start values
t = 0:0.01:10;   % Time-span
N = length(t);          % Length of time
params = params;
q = [1; 1; pi/6; 0.5; 0.5; pi/3;];  % Guess Values

% Finding Position, velocity and acceleration
[q,q_dot,q_ddot] = NewtonRaphston(t,q);
pp = 100; % Percentage into the simulation
p = round(N*(pp/100)); % Point in simulation
% Geometry of arms
A = @(theta) [cos(theta) -sin(theta); sin(theta) cos(theta)];   %   A-Matrix
B = @(theta) [-sin(theta) -cos(theta); cos(theta) -sin(theta)]; %   B-Matrix
% Finding points
for i = 1:N
    P1(:,i) = [q(1,i);q(2,i)] - A(q(3,i))*[2;0]; % Mid point of circle
    P2(:,i) = [q(1,i);q(2,i)] + A(q(3,i))*[2;0]; % Circle perimeter
    P3(:,i) = [q(4,i);q(5,i)] - A(q(6,i))*[8;0]; % Start point long arm
    P4(:,i) = [q(4,i);q(5,i)] + A(q(6,i))*[8;0]; % End point long arm
    P5(:,i) = P2(:,i)+A(q(6,i))*[0;2];           % Trans point on Long arm
end
for i = 1:N
    part1(:,i) = q_ddot(1:2,i).';
    part2(:,i) = q_dot(6,i)^2*A(q(6,i))*[8;0];
    part3(:,i) = q_ddot(6,i)*B(q(6,i))*[8;0];
    acc(:,i) = part1(:,i)-part2(:,i)+part3(:,i);
end
L_2 = params.L_2;
s2p = [L_2/2;0];
i = 1;
meh = zeros(2,N);
blyat = zeros(2,N);
for i = 1:N
    meh(:,i) = q_ddot(4:5,i)-q_dot(6,i)^2*A(q(6,i))*s2p+q_ddot(6,i)*B(q(6,i))*s2p;
    blyat(:,i) = norm(meh(:,i));
    ewq(i) = L_2/2*q_ddot(6,i);
end
blyat = blyat(1,:)
max(blyat)
max(ewq)
figure(8)
plot(t,blyat,t,ewq)
hold on
legend('Resultant','Angular')
title('Acceleration for point E')
xlabel('Time (s)')
ylabel('Acceleration (m/s^2)')
hold off

% Plotting Trajectory
figure(1)
plot(P1(1,:),P1(2,:),'r')
hold on
plot(P2(1,:),P2(2,:),'r')
plot(P3(1,:),P3(2,:),'m')
plot(P4(1,:),P4(2,:),'b')
plot(P1(1,p),P1(2,p),'ro-')
plot(P3(1,p),P3(2,p),'bo-')
plot(P5(1,p),P5(2,p),'mo-')
plot(P2(1,p),P2(2,p),'mo-')
plot(P4(1,p),P4(2,p),'bo-')
Body1 = line([P1(1,p) P2(1,p)], [P1(2,p) P2(2,p)]);
Body1.Color = 'r';
Body2 = line([P2(1,p) P5(1,p)], [P2(2,p) P5(2,p)]);
Body2.Color = 'm';
Body3 = line([P3(1,p) P4(1,p)], [P3(2,p) P4(2,p)]);
Body3.Color = 'b';
xlim([-9 9])
ylim([-5 13])
title('Trajectory')
xlabel('X')
ylabel('Y')
hold off
% 
% % Live Trajectory
% figure(8)
% plot(P1(1,:),P1(2,:),'r');
% hold on
% plot(P2(1,:),P2(2,:),'r');
% plot(P3(1,:),P3(2,:),'m');
% plot(P4(1,:),P4(2,:),'b');
% hold on
% title('Live Trajectory')
% xlabel('X-axis')
% ylabel('Y-axis')
% xlim([-10 10])
% ylim([-5 13])
% for p = 1:p
%     Body1 = line([P1(1,p) P2(1,p)], [P1(2,p) P2(2,p)]);
%     Body1.Color = 'r';
%     Body2 = line([P2(1,p) P5(1,p)], [P2(2,p) P5(2,p)]);
%     Body2.Color = 'm';
%     Body3 = line([P3(1,p) P4(1,p)], [P3(2,p) P4(2,p)]);
%     Body3.Color = 'b';
%     pause(0.01)
%     delete(Body1);
%     delete(Body2);
%     delete(Body3);
% end
% plot(P1(1,p),P1(2,p),'ro-')
% plot(P3(1,p),P3(2,p),'bo-')
% plot(P5(1,p),P5(2,p),'mo-')
% plot(P2(1,p),P2(2,p),'mo-')
% plot(P4(1,p),P4(2,p),'bo-')
% Body1 = line([P1(1,p) P2(1,p)], [P1(2,p) P2(2,p)]);
% Body1.Color = 'r';
% Body2 = line([P2(1,p) P5(1,p)], [P2(2,p) P5(2,p)]);
% Body2.Color = 'm';
% Body3 = line([P3(1,p) P4(1,p)], [P3(2,p) P4(2,p)]);
% Body3.Color = 'b';
% hold off
% 
% % Position - Arm 1
figure(2)
sgtitle('Arm 1 - Position')
subplot(2,2,1)
plot(t,q(1,:))
title('X')
xlabel('Time (s)')
ylabel('Position (m)')
subplot(2,2,2)
plot(t,q(2,:))
title('Y')
xlabel('Time (s)')
ylabel('Position (m)')
subplot(2,2,[3 4])
plot(t,q(3,:))
title('\theta_1')
xlabel('Time (s)')
ylabel('Orientation (Rad)')
% 
% % Position - Arm 2
figure(3)
sgtitle('Arm 2 - Position')
subplot(2,2,1)
plot(t,q(4,:))
title('Position')
xlabel('Time (s)')
ylabel('X (m)')
subplot(2,2,2)
plot(t,q(5,:))
title('Position')
xlabel('Time (s)')
ylabel('Y (m)')
subplot(2,2,[3 4])
plot(t,q(6,:))
title('\theta_1')
xlabel('Time (s)')
ylabel('Orientation (Rad)')

% Velocity - Arm 1
figure(4)
sgtitle('Arm 1 - Velocity')
subplot(2,2,1)
plot(t,q_dot(1,:))
title('dX')
xlabel('Time (s)')
ylabel('Velocity (m/s)')
subplot(2,2,2)
plot(t,q_dot(2,:))
title('dY')
xlabel('Time (s)')
ylabel('Velocity (m/s)')
subplot(2,2,[3 4])
plot(t,q_dot(3,:))
title('\theta_1')
xlabel('Time (s)')
ylabel('Angular Speed (Rad/s)')

% Velocity - Arm 2
figure(5)
sgtitle('Arm 2 - Velocity')
subplot(2,2,1)
plot(t,q_dot(4,:))
title('X')
xlabel('Time (s)')
ylabel('Velocity (m/s)')
subplot(2,2,2)
plot(t,q_dot(5,:))
title('Y')
xlabel('Time (s)')
ylabel('Velocity (m/s)')
subplot(2,2,[3 4])
plot(t,q_dot(6,:))
title('d\theta_1')
xlabel('Time (s)')
ylabel('Angular Speed (Rad/s)')

% Acceleration - Arm 1
figure(6)
sgtitle('Arm 1 - Acceleration')
subplot(2,2,1)
plot(t,q_ddot(1,:))
title('ddX')
xlabel('Time (s)')
ylabel('Acceleration (m/s^2)')
subplot(2,2,2)
plot(t,q_ddot(2,:))
title('ddY')
xlabel('Time (s)')
ylabel('Acceleration (m/s^2)')
subplot(2,2,[3 4])
plot(t,q_ddot(3,:))
title('dd\theta_1')
xlabel('Time (s)')
ylabel('Angular Acceleration (Rad/s^2)')

% Acceleration - Arm 2
figure(7)
sgtitle('Arm 2 - Acceleration')
subplot(2,2,1)
plot(t,q_ddot(4,:))
title('ddX')
xlabel('Time (s)')
ylabel('Acceleration (m/s^2)')
subplot(2,2,2)
plot(t,q_ddot(5,:))
title('ddY')
xlabel('Time (s)')
ylabel('Acceleration (m/s^2)')
subplot(2,2,[3 4])
plot(t,q_ddot(6,:))
title('dd\theta_1')
xlabel('Time (s)')
ylabel('Angular Acceleration (Rad/s^2)')