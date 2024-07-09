function res_ = coveigAbcont(A,B,C,G,dA,dB,dC,nmax,dt)

    In = speye(nmax);
    Ir = speye(size(B,2));
    
    [V,DD] = eig(A);
        
    Ac = logm(A)/dt;
    Bc = Ac*inv(A-In)*B;
%     Bc = inv(A-In)*Ac*B;
    
%     Ab = Ac + Bc*G*C;
    Ab = Ac - Bc*G*C;
    
    J_lambda_all = zeros(nmax^2,nmax^2);
    J_phi_H_all = zeros(nmax^2,nmax^2);

    D = diag(DD);
    U = inv(V)';
    
    Dc = log(D)/dt;
    DDc = diag(Dc);
    
    for i = 1:nmax
        psi = V(:,i);
        chi = U(:,i);
        sval = D(i);

        J_lambda_H = 1/(sval*dt) * 1/(chi'*psi) * kron(psi.',chi');
        J_lambda_all((i-1)*nmax+i:(i-1)*nmax+i,:) = J_lambda_H;

        J_phi_H = pinv(sval*In-A)*kron(psi.',( In - (psi*chi')/(chi'*psi)));
        J_phi_H_all((i-1)*nmax+1:i*nmax,:) = J_phi_H;
    end 

    J_A_Ad = kron(inv(V).',V)*J_lambda_all + ( kron((DDc*inv(V)).',In) - kron(inv(V).',Ac))*J_phi_H_all;
%     
    J_B_A1 = kron((inv(A-In)*B).',In);
    J_B_A2 = kron(B.',Ac);
    J_B_A3 = -kron(inv(A-In).',inv(A-In));
    
    J_B_Ad = J_B_A1*J_A_Ad+J_B_A2*J_B_A3;
    
    J_B_Bd = kron(Ir,Ac*inv(A-In));
    
    J_Ab_ABC = [J_A_Ad-kron((G*C).',In)*J_B_Ad -kron((G*C).',In)*J_B_Bd -kron(In,Bc*G)];
%     J_Ab_ABC = [J_A_Ad+kron((G*C).',In)*J_B_Ad kron((G*C).',In)*J_B_Bd kron(In,Bc*G)];
    
    d_ABC = [dA;dB;dC];
    
%     % BEWARE OF THE COMPLEX CALCULUS
    
    cov_ABC = d_ABC*d_ABC';

    [r_vec,e_val] = eig(Ab);

    e_vald = diag(e_val);
    l_vec = inv(r_vec)';
    
    for i = 1:length(e_vald)

        psi = r_vec(:,i);
        chi = l_vec(:,i);

        J_l_Ab = 1/(chi'*psi) * kron(psi.',chi');

%         J_psi_Ab = pinv(e_val(i)*In-Ab)*kron(psi.',( In - (psi*chi')/(chi'*psi)));
%         J_chi_Ab = pinv(conj(e_val(i))*In-Ab')*( kron(chi.',In)*Pn - (Ab'*chi/(conj(e_val(i)))) * conj(J_l_Ab));    
        
        mm_J_l_Ab = [real(J_l_Ab*J_Ab_ABC);imag(J_l_Ab*J_Ab_ABC)];
        
        cov_lambda_Ab{i} = mm_J_l_Ab*cov_ABC*mm_J_l_Ab';

    end    
    
    res_.cov_lambda_Ab = cov_lambda_Ab; %eigenvalue covariance 2x2 [Re ReIm; ImRe Im]

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
