function [A,K,D] = Amat(params,M,theta)
K_beam = params.K_beam;
Kg_beam = params.Kg_beam;
zeta_nom = params.zeta_nom;
K = Kmat(K_beam,Kg_beam,theta);
D = Dmat(M,K,zeta_nom);
A = [zeros(length(M)),eye(length(M));
     -inv(M)*K, -inv(M)*D];
end