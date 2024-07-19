function dy = EOM(t,y,params)
m2 = params.m2;
m3 = params.m3;
J2 = params.J2;
J3 = params.J3;
g  = params.g;
K0 = params.K0;
C0 = params.C0;
persistent y0;
q = y(1:6);
q_dot = y(7:12);
lambda = y(13:17);


if t == 0
    y0 = q(6);
end
n_1 = K0*(q(6)-y0) + C0*q_dot(6);

M = diag([m2 m2 J2 m3 m3 J3]);                                     % Mass Matrix
Q_A = [0 -m2*g 0 0 -m3*g -n_1].';                                     % Applied Force
A = [M Phi_q(params,q).'; Phi_q(params,q) zeros(5)];               % Combined Matrix
B = [Q_A;Gamma(params,q,q_dot)];                                   % Combined Matrix
C = A\B;                                                           % Acc,Lambda
q_ddot = C(1:11);                                                  % Retrieving acc
dy = [q_dot;q_ddot];                                               % Return value
end