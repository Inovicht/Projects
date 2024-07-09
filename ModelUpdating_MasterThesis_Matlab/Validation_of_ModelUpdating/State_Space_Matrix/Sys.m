function sys = Sys(params,theta) 
%This function generates a structure that contain the system quadruplet:
%The state matrix A, the input matrix B, the output matrix C, and the
%transmission matrix D. The function can also do modal truncation to lower model
%order. REMEMBER TO KEEP EVERYTHING IN SI-UNITS

%input: 
% params:structure that contain the youngs modulus, crosssectional
%area, length, and area moment of inertia of each element

%theta: 5x1 vector that contain the boundary stiffnesses and tension of
%beam. [k1,k2,k3,k4,N]

[params.K_beam,params.Kg_beam,sys.M] =beam2D(params,theta); %Computing global stiffness matrix K_beam, global geometric stiffness Kg_beam, and global mass matrix M.
[sys.A,sys.K,sys.D] = Amat(params,sys.M,theta);             %Generation of state matrix
[sys.Cc,sys.Ca] = Cmat(sys.M,sys.K,sys.D,params.d,params.ddot,params.dddot);    %Generation of output matrix
[sys.B,sys.B2] = Bmat(sys.M,params.q);                                          %Generation of input matrix
if (isempty(params.dddot))                                                      %Generation of Transmission matrix
    sys.H = 0;
else
    sys.H = Hmat(sys.Ca,sys.M,sys.B2);
end
sys = Truncation(sys,params.n_trunk);