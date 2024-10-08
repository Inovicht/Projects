function phiq = Phi_q(params,q)
%PHI_Q
%    PHIQ = PHI_Q(L_2,L_3,THETA_2,THETA_3,X_2,Y_2)

%    This function was generated by the Symbolic Math Toolbox version 8.7.
%    08-May-2022 16:02:29

L_2 = params.L2;
L_3 = params.L3;
theta_2 = q(3);
theta_3 = q(6);
x_2 = q(1);
y_2 = q(2);


t2 = cos(theta_2);
t3 = cos(theta_3);
t4 = sin(theta_2);
t5 = sin(theta_3);
t6 = (L_2.*t2)./2.0;
t7 = (L_3.*t3)./2.0;
t8 = (L_2.*t4)./2.0;
t9 = (L_3.*t5)./2.0;
t10 = -t8;
t11 = -t9;
phiq = reshape([x_2.*2.0-L_2.*t2,1.0,0.0,0.0,0.0,y_2.*2.0-L_2.*t4,0.0,1.0,0.0,0.0,-L_2.*t4.*(t6-x_2)+L_2.*t2.*(t8-y_2),t10,t6,0.0,0.0,0.0,-1.0,0.0,1.0,0.0,0.0,0.0,-1.0,0.0,1.0,0.0,t11,t7,t11,t7],[5,6]);
