close all; clear all;
tspan = 0:0.001:5;
params = params;
L_2 = params.L2;
L_3 = params.L3;
y0 = zeros(1,17);
%q_guess = [-2; 0.5; pi/4; 1; 1; pi/10];

grashof = [params.L1,params.L2,params.L3,params.G];
if max(grashof)+min(grashof) < sum(grashof)-max(grashof)+min(grashof) % See if Grashof is satisfied
    %q_guess = NewtonRaphston(tspan(1),q_guess)
else
    error("Grashof is not satisfied")
end
q_guess =  [1.0036; 2.2393; 0; 5.0018; 1.1197; 5.4405];
y0 = [q_guess;zeros(11,1)];
% Runge Kutta with correction
opts = odeset; opts=odeset(opts,'RelTol',1e-13,'AbsTol',1e-14);
[t y] = ode45(@(t,y) EOM(t,y,params),tspan,y0',opts);
% Position, velocity, acceleration
q = y(:,1:6);
q_dot = y(:,6:12);
B = @(theta) [-sin(theta) -sin(theta); cos(theta) -sin(theta)];
% Acceeleration, Constraint Forces, Constraint Torso
dy = zeros(length(tspan),17);
lambda = zeros(5,length(tspan));
for i = 1:length(tspan)
    dy(i,:) = EOM(tspan(i),y(i,:)',params);
    lambda(:,i) = dy(i,13:17);
    CF(i,:) = -Phi_q(params,q(i,:)).'*lambda(:,i);
    phik = Phi_q(params,q(i,:));
    CT(i,1) = ([L_2/2;0].'*B(q(i,3)).'*phik(2:3,1:2).'-phik(2:3,3).')*lambda(2:3,i);
    CT(i,2) = ([-L_3/2;0].'*B(q(i,6)).'*phik(2:3,4:5).'-phik(2:3,6).')*lambda(2:3,i);
end
q_ddot = dy(:,7:12);

% Constraint forces - Arm 1
figure(14)
plot(t,CT(:,1),'b',t,CT(:,2),'r')
hold on
axis padded
grid on
legend('Arm_1','Arm_2')
title('Constraint Moment')
ylabel('Moment (mN)')
xlabel('Time (s)')
hold off


% Spring-Damper Forces
K0 = params.K0;
C0 = params.C0;
Spring_Force = zeros(length(tspan));
Damper_Force_Force = zeros(length(tspan));
for i = 1:length(tspan)
    Spring_Force(i) = K0*(q(i,6)-y0(6));
    Damper_Force(i) = C0*q_dot(i,6);
end

% Getting points
A = @(theta) [cos(theta) -sin(theta); sin(theta) cos(theta)];
for i = 1:length(tspan)
    P1(i,:) = [q(i,1);q(i,2)] + A(q(i,3))*[-L_2/2;0]; % Second arm left point
    P2(i,:) = [q(i,1);q(i,2)] + A(q(i,3))*[L_2/2;0];  % Second arm right point
    P3(i,:) = [q(i,4);q(i,5)] + A(q(i,6))*[-L_3/2;0]; % Third arm left point
    P4(i,:) = [q(i,4);q(i,5)] + A(q(i,6))*[L_3/2;0];  % Third arm right point
end

figure(1)% Live Trajectory
hold on
grid on
title('Live Trajectory')
xlabel('X-axis')
ylabel('Y-axis')
xlim([-4 8])
ylim([-4 5])
tic
for p = 1:length(tspan)
    Body1 = line([P1(p,1) P2(p,1)], [P1(p,2) P2(p,2)]);
    Body1.Color = 'r';
    Body2 = line([P3(p,1) P4(p,1)], [P3(p,2) P4(p,2)]);
    Body2.Color = 'b';
    Body3 = line([0 P1(p,1)], [0 P1(p,2)]);
    Body3.Color = 'm';
    pause(0.000001)
    delete(Body1);
    delete(Body2);
    delete(Body3);
end
toc
plot(q(p,1),q(p,2),'ro-');
plot(q(p,4),q(p,5),'bo-');
Body1 = line([P1(end,1) P2(end,1)], [P1(end,2) P2(end,2)]);
Body1.Color = 'r';
Body2 = line([P3(end,1) P4(end,1)], [P3(end,2) P4(end,2)]);
Body2.Color = 'b';
Body3 = line([0 P1(end,1)], [0 P1(end,2)]);
Body3.Color = 'm';
hold off
% 
% Trajectory
% figure(2)
% plot(P1(:,1),P1(:,2),'.','Color','b');
% hold on
% plot(P2(:,1),P2(:,2),'b');
% plot(P3(:,1),P3(:,2),'.','Color','r');
% plot(P4(:,1),P4(:,2),'r');
% plot(q(end,1),q(end,2),'bo-');
% plot(q(end,4),q(end,5),'ro-');
% axis equal
% grid on
% axis padded
% title('Trajectory')
% xlabel('X-axis')
% ylabel('Y-axis')
% Body3 = line([0 P1(end,1)], [0 P1(end,2)]);
% Body3.Color = 'm';
% Body1 = line([P1(end,1) P2(end,1)], [P1(end,2) P2(end,2)]);
% Body1.Color = 'b';
% Body2 = line([P3(end,1) P4(end,1)], [P3(end,2) P4(end,2)]);
% Body2.Color = 'r';
% hold off

% Constraint forces - Arm 1
% figure(9)
% plot(t,CF(:,1),'b',t,CF(:,2),'r',t,CF(:,3),'m')
% hold on
% axis padded
% grid on
% legend('x','y','\theta')
% title('Constraint Forces - Arm 1')
% ylabel('Forces (N)')
% xlabel('Time (s)')
% hold off
% 
% % Constraint forces - Arm 2
% figure(10)
% plot(t,CF(:,4),'b',t,CF(:,5),'r',t,CF(:,6),'m')
% hold on
% axis padded
% grid on
% legend('x','y','\theta')
% title('Constraint Forces - Arm 2')
% ylabel('Forces (N)')
% xlabel('Time (s)')
% hold off

% Damper Forces
figure(15)
plot(t,Damper_Force,'r')
hold on
axis padded
grid on
legend('Damper')
title('Damper - Arm 2')
ylabel('Forces (N)')
xlabel('Time (s)')
hold off

% Spring Force
% figure(11)
% plot(t,Spring_Force,'b')
% hold on
% axis padded
% grid on
% legend('Spring')
% title('Spring - Arm 2')
% ylabel('Forces (N)')
% xlabel('Time (s)')
% hold off

% % Position - Arm 1
% figure(3)
% sgtitle('Arm 1 - Position')
% subplot(2,2,1,'Color','r')
% plot(t,q(:,1))
% grid on
% axis padded
% title('Position (m)')
% xlabel('Time (s)')
% ylabel('X-axis')
% subplot(2,2,2)
% plot(t,q(:,2))
% grid on
% axis padded
% title('Position (m)')
% xlabel('Time (s)')
% ylabel('Y-axis')
% subplot(2,2,[3 4],'Color','b')
% plot(t,q(:,3))
% grid on
% axis padded
% title('\theta_1')
% xlabel('Time (s)')
% ylabel('Orientation (Rad)')
% 
% % Velocity - Arm 1
% figure(4)
% sgtitle('Arm 1 - Velocity')
% grid on
% axis padded
% subplot(2,2,1)
% plot(t,q_dot(:,1))
% grid on
% axis padded
% title('dX')
% xlabel('Time (s)')
% ylabel('Velocity (m/s)')
% subplot(2,2,2)
% plot(t,q_dot(:,2))
% grid on
% axis padded
% title('dY')
% xlabel('Time (s)')
% ylabel('Velocity (m/s)')
% subplot(2,2,[3 4])
% plot(t,q_dot(:,3))
% grid on
% axis padded
% title('\theta_1')
% xlabel('Time (s)')
% ylabel('Angular Speed (Rad/s)')
% 
% % Acceleration - Arm 1
% figure(5)
% sgtitle('Arm 1 - Acceleration')
% grid on
% axis padded
% subplot(2,2,1)
% plot(t,q_ddot(:,1))
% grid on
% axis padded
% title('ddX')
% xlabel('Time (s)')
% ylabel('Acceleration (m/s^2)')
% subplot(2,2,2)
% plot(t,q_ddot(:,2))
% grid on
% axis padded
% title('ddY')
% xlabel('Time (s)')
% ylabel('Acceleration (m/s^2)')
% subplot(2,2,[3 4])
% plot(t,q_ddot(:,3))
% grid on
% axis padded
% title('dd\theta_1')
% xlabel('Time (s)')
% ylabel('Angular Acceleration (Rad/s^2)')
% 
% % Position - Arm 2
% figure(6)
% sgtitle('Arm 2 - Position')
% subplot(2,2,1,'Color','r')
% plot(t,q(:,4))
% grid on
% axis padded
% title('Position (m)')
% xlabel('Time (s)')
% ylabel('X-axis')
% subplot(2,2,2)
% plot(t,q(:,5))
% grid on
% axis padded
% title('Position (m)')
% xlabel('Time (s)')
% ylabel('Y-axis')
% subplot(2,2,[3 4],'Color','b')
% plot(t,q(:,6))
% grid on
% axis padded
% title('\theta_1')
% xlabel('Time (s)')
% ylabel('Orientation (Rad)')
% 
% % Velocity - Arm 2
% figure(7)
% sgtitle('Arm 2 - Velocity')
% grid on
% axis padded
% subplot(2,2,1)
% plot(t,q_dot(:,4))
% grid on
% axis padded
% title('dX')
% xlabel('Time (s)')
% ylabel('Velocity (m/s)')
% subplot(2,2,2)
% plot(t,q_dot(:,5))
% grid on
% axis padded
% title('dY')
% xlabel('Time (s)')
% ylabel('Velocity (m/s)')
% subplot(2,2,[3 4])
% plot(t,q_dot(:,6))
% grid on
% axis padded
% title('\theta_1')
% xlabel('Time (s)')
% ylabel('Angular Speed (Rad/s)')
% 
% % Acceleration - Arm 2
% figure(8)
% sgtitle('Arm 2 - Acceleration')
% grid on
% axis padded
% subplot(2,2,1)
% plot(t,q_ddot(:,4))
% grid on
% axis padded
% title('ddX')
% xlabel('Time (s)')
% ylabel('Acceleration (m/s^2)')
% subplot(2,2,2)
% plot(t,q_ddot(:,5))
% grid on
% axis padded
% title('ddY')
% xlabel('Time (s)')
% ylabel('Acceleration (m/s^2)')
% subplot(2,2,[3 4])
% plot(t,q_ddot(:,6))
% grid on
% axis padded
% title('dd\theta_1')
% xlabel('Time (s)')
% ylabel('Angular Acceleration (Rad/s^2)')

