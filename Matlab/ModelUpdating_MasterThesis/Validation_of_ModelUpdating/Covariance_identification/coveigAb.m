function res_ = coveigAb(A,B,C,G,dA,dB,dC,nmax)

    In = speye(nmax);
    Pn = permutation_matrix(nmax,nmax);    

    Ab = A - B*G*C;

    dAb = dA - kron((G*C)',In)*dB - kron(In,B*G)*dC;

    [r_vec,e_val,l_vec] = eig(Ab);

    e_vald = diag(e_val);

    for i = 1:length(e_vald)

        psi = r_vec(:,i);
        chi = l_vec(:,i);

        J_l_Ab = 1/(chi'*psi) * kron(psi.',chi');

        J_psi_Ab = pinv(e_val(i)*In-Ab)*kron(psi.',( In - (psi*chi')/(chi'*psi)));
        J_chi_Ab = pinv(conj(e_val(i))*In-Ab')*( kron(chi.',In)*Pn - (Ab'*chi/(conj(e_val(i)))) * conj(J_l_Ab));            

        J_absl_Ab = 1/abs(e_vald(i))*[2*real(e_vald(i)) 2*imag(e_vald(i))]*...
            [real(J_l_Ab);imag(J_l_Ab)]*dAb;

        cov_abslambda_Ab(i) = J_absl_Ab*J_absl_Ab';
        
        J_psi_Ab = [real(J_psi_Ab*dAb);imag(J_psi_Ab*dAb)];
        J_chi_Ab = [real(J_chi_Ab*dAb);imag(J_chi_Ab*dAb)];

        cov_psi_Ab{i}= J_psi_Ab;
        cov_chi_Ab{i}= J_psi_Ab;

    end    
    
    res_.cov_abslambda_Ab = cov_abslambda_Ab; %eigenvalue covariance
    res_.cov_psi_Ab = cov_psi_Ab; %right eigenvector covariance
    res_.cov_chi_Ab = cov_chi_Ab; %left eigenvector covariance

    res_.e_vald = e_vald; %eigenvalues
    res_.r_vec = r_vec; %right eigenvectors
    res_.l_vec = l_vec; %left eigenvectors
    
end


function P = permutation_matrix(m,n)
    P = [];
    for ik = 1:m
        e_m = zeros(1,m);
        e_m(1,ik) = 1;
        e_m = sparse(e_m);
%             P = [P kron(speye(a),e_b)];
        P = [P;kron(speye(n),e_m)];
    end


end  
