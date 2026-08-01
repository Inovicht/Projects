function param = params
%double_pendulum_simulation
m = [1,1];
L = [1;1];
g = 9.81; % Gravitational acceleration
init = [0, 0, pi, 0];
tspan = [0 5];

M = [L(1)^2*m(1) + L(1)^2*m(2) , L(1)*L(2)*m(2);
     L(1)*L(2)*m(2)          , L(2)^2*m(2)];
K = g*[L(1)*m(2)+L(1)*m(2)     , 0
       0               , L(2)*m(2)];

end

