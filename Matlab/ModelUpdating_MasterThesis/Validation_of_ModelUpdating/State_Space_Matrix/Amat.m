function [A,K,D] = Amat(params,M,theta)
%This function computes the state matrix A, including with boundary
%confitions of the stiffness matrix

K_beam = params.K_beam;                %Global stiffness of the beam
Kg_beam = params.Kg_beam;              %Global geometric stiffness of the beam
zeta_nom = params.zeta_nom;            %model
K = Kmat(K_beam,Kg_beam,theta);        %Computing Global total stiffness matrix with boundary condiditions
D = Dmat(M,K,zeta_nom);                %Computing damping matrix
A = [zeros(length(M)),eye(length(M));  %Computing state matrix A
     -inv(M)*K, -inv(M)*D];
end