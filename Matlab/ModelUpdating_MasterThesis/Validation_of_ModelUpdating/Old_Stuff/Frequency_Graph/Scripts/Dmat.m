function D_mat = Dmat(M,K_nom,zeta_mod)
D_nom = zeros(length(M));
[Phi,lambda] = eig(K_nom,M);
omega = diag(sqrt(lambda));
D_tilde_nom = diag(2*omega.*zeta_mod);
D_mat = (Phi.')\D_tilde_nom/(Phi);
end

