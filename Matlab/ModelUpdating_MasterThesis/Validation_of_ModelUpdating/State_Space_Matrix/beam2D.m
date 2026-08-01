%% Function to conduct element assembly
function [K,Kg,M]=beam2D(params,theta)
%The function generates the global mass, stiffness and geometric stiffness
%of the beam structure without implementation of boundary conditions.

%input:
%params structure that contains: 
    %n_nodes*1 vector that contains the youngs modulus E for each element
    %n_nodes*1 vector that contains the cross sectional area A for each element
    %n_nodes*1 vector that contains the area moment of inertia I for each element
    %n_nodes*2 matrix that contains the coordinates of each node
    %n_nodes*1 vector that contains the density rho for each element
    %row vector that contains all the degrees of freedom that are fixed. In
    %our case, it is the axial degrees of freedom. 
%theta: 5*1 vector that contain the boundary springs and tension of the
%beam [k1,k2,k3,k4,N]

%Output:
%K:Global stiffness matrix
%Kg: Global geometrix stiffness
%M: Global mass matrix


fixeddof = params.fixeddof;
coor = params.coor;
no = params.no;
if size(no,1)==1
    dofe=1:6;
else
    dofs=reshape(1:(max(no(:,2))*3),3,[])';
    dofno=num2cell(dofs,2);
    dofe=cell2mat(dofno(no(:,1:2)));
end
ndof=max(max(dofe));
Kglobal = zeros(ndof,ndof);
Kgglobal = zeros(ndof,ndof);
Mglobal = zeros(ndof,ndof);
for i=1:size(no,1)
    L = sqrt((coor(no(i,2),1)-coor(no(i,1),1))^2 + ...
             (coor(no(i,2),2)-coor(no(i,1),2))^2);
    I = no(i,3);
    E = no(i,4);
    A = no(i,5);
    rho = no(i,6);
    
    Ke = [ E*A/L        0           0               -E*A/L      0               0;
           0            12*E*I/L^3  6*E*I/L^2       0           -12*E*I/L^3     6*E*I/L^2;
           0            6*E*I/L^2   4*E*I/L         0           -6*E*I/L^2      2*E*I/L;
           -E*A/L       0           0               E*A/L       0               0;
           0            -12*E*I/L^3 -6*E*I/L^2      0           12*E*I/L^3      -6*E*I/L^2;
           0            6*E*I/L^2   2*E*I/L         0           -6*E*I/L^2      4*E*I/L ];
    
    Me = rho*A*L/420 * [ 140    0       0       70      0       0;
                         0      156     22*L    0       54      -13*L;
                         0      22*L    4*L^2   0       13*L    -3*L^2;
                         70     0       0       140     0       0;
                         0      54      13*L    0       156     -22*L;
                         0      -13*L   -3*L^2  0       -22*L   4*L^2 ];
    if theta(end) > 0
        Kg = ...
            [0,     0,          0,          0,      0,              0;
             0,  (6/5)*(1/L),   1/10,       0,      -(6/5)*(1/L),   1/10;
             0,    1/10,        (2/15)*L,   0,      -1/10,          -(1/30)*L;
             0,     0,          0,          0,      0,              0;
             0, -(6/5)*(1/L),   -1/10,      0,      (6/5)*(1/L),    -1/10;
             0,     1/10,       -(1/30)*L,  0,     -1/10,           (2/15)*L;
             ];
    else
        Kg = zeros(length(Ke),length(Ke));
    end
    
    Kglobal(dofe(i,:),dofe(i,:))=Ke+Kglobal(dofe(i,:),dofe(i,:));
    Kgglobal(dofe(i,:),dofe(i,:))= Kg+Kgglobal(dofe(i,:),dofe(i,:));
    Mglobal(dofe(i,:),dofe(i,:))=Me+Mglobal(dofe(i,:),dofe(i,:));
end
    totdof=1:ndof;
    freedof=totdof;
    freedof(fixeddof)=[];
    K=Kglobal(freedof,freedof);
    Kg=Kgglobal(freedof,freedof);
    M=Mglobal(freedof,freedof);
end