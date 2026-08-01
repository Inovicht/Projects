function D_mat = Dmat(M,K_nom,zeta_mod)
%function that computes the damping matrix. Modal damping is assumed

%inputs:
%Global mass matrix M
%Global total stiffness matrix with boundary conditions K_nom
%vector containing the modal damping of each DoF

[Phi,lambda] = eig(K_nom,M);                %Estimation of eigenvectors and eigenvalues
omega = diag(sqrt(lambda));
D_tilde_nom = diag(2*omega.*zeta_mod);      %Computing modal damping matrix
D_mat = (Phi.')\D_tilde_nom/(Phi);          %Transformation of damping matrix to physical coordinates
end

