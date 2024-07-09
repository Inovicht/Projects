% I/O system identification 
% Ask your advisor for further explanation of this function
function res_ = system_id(Y,param)

    nmax = param.nmax;
    p = param.p;
    q = param.q;
    r = param.r;
    u = param.inputs;
    
    s_p1r = q*r;
    
    % identity matrices;
    Iqrmn = speye(q*r-nmax);
    Iqqrmn = speye(q*(q*r-nmax));
    Irpn = speye(r+nmax);
    Iu = speye(u);
    In = speye(nmax);  
    In2 = speye(nmax^2);  
    Pn = permutation_matrix(nmax,nmax);    
    
    % selection matrices;
    S1 = s1matrix(q,r,u,nmax);
    S4 = s4matrix(q,r,nmax);
    S6 = [sparse(r,(q-1)*r) ; speye( (q-1)*r)] * [speye( (q-1)*r) sparse((q-1)*r,r)];
    S7 = [sparse(q*r^2,q*r*nmax); kron( speye(nmax), S6 )];    
    
    S11 = [eye(s_p1r-r) zeros(s_p1r-r,r)];
    S22 = [zeros(s_p1r-r,r) eye(s_p1r-r)];      
    
    [covSSI,Htot,Ctot,hk]=covvecSSI(Y,param);

    % B and D computation; 
    %
    [U,S,V] = svd(hk,'econ');
    U1 = U(:,1:nmax);
    V1 = V(:,1:nmax);
    U2 = U(:,nmax+1:end);
    s1 = diag(S);
    s1 = s1(1:nmax);
    Ob  =  U(:,1:nmax);
    C  =  Ob(1:r,:);
    G_up  =  Ob(1:p*r,:);
    G_dn = Ob(r+1:s_p1r,:);
    A = pinv(G_up)*G_dn;
    %
    % Determine the matrices M and U2T
    R4 = Htot.YfUf;
    iR8 = inv(Htot.UfUf);
    %
    U2T = U2';
    M = U2T*((Htot.YfUf)*iR8); % U2^T R_4 inv(R_8)

    % Determine the set of equations
    Mv = zeros(q*(r*q-nmax),u);
    L = zeros(q*(r*q-nmax),q*r);
    for k=1:q
      Mv((k-1)*(r*q-nmax)+1:k*(r*q-nmax),:) = M(:,(k-1)*u+1:k*u);
      L((k-1)*(r*q-nmax)+1:k*(r*q-nmax),1:(q-k+1)*r) = U2T(:,(k-1)*r+1:r*q);
    end
    Os = [eye(r),zeros(r,nmax);zeros(r*(q-1),r),G_up];
    Ls = L*Os;
    pinvLs = pinv(Ls);
    %
    % Solve least squares
    sol = pinvLs*Mv;

    % Extract the system matrices
    D = sol(1:r,:);
    B = sol(r+1:r+nmax,:);
    %
    %% compute the jacobians for B and D
    
    % Left kernel transposed w.r.t. to subspace matrix
    iD_s = diag(1./s1);
    p1 = (U1*iD_s*V1');        
    J_H_UkerT = -kron(p1,U2T);
        
    % transpose the jacobian for later
    Pab = permutation_matrix(size(U2,2),size(U2,1)); % corrected 
    J_H_Uker = Pab*J_H_UkerT;
    
    % assembled jacobian for the first part of Delta [D;B] 
    %    
    p1 = kron((R4*iR8)',Iqrmn) * J_H_UkerT;
    R_1 = kron(iR8',U2T);
    R_2 = kron(iR8',-U2T*R4*iR8); 
    
    J_M = kron(Iu,pinvLs);    
    J_R_Mv = J_M*S1*p1*covSSI;
    
    R_part = J_M*S1*R_1*Ctot.YfUf + J_M*S1*R_2*Ctot.UfUf;
    
    %% Second part 
    
    J_L = kron(Mv',Irpn);    
    J_L_Ls = kron(Os',Iqqrmn);
    J_Os_Ls = kron(Irpn,L);
    
    Pef = permutation_matrix(size(L,2),size(L,1)); % corrected  
    
    [J_O_R] = deltaO_H_clean(p,r,nmax,V1,U2,s1,covSSI);
    
    J_A_O = kron(In,inv(G_up'*G_up))*(-(kron((A)',In))*(Pn+In2)*(kron(In,G_up'*S11)) + Pn*(kron(In,G_dn'*S11)) + kron(In,G_up'*S22) );
    J_A_R = J_A_O*J_O_R;

    J_C_O = kron(In,[eye(r) zeros(r,p*r)]);
    J_C_R = J_C_O*J_O_R;

    J_pinvLs_Ls = JJ_pinvLs_Ls(Ls,pinvLs);
    
    J_R_pLs = J_L*J_pinvLs_Ls*J_L_Ls*Pef*S4*J_H_Uker*covSSI + J_L*J_pinvLs_Ls*J_Os_Ls*S7*J_O_R;
    
    
    %% Sensitivity [D;B]
    
    J_BD_R = J_R_Mv + J_R_pLs + R_part;
    
    Sel_D = kron(speye(u),[speye(r) zeros(r,nmax)]);
    Sel_B = kron(speye(u),[zeros(nmax,r) speye(nmax)]);
    
    J_D_R = Sel_D*J_BD_R;
    J_B_R = Sel_B*J_BD_R;
    
    
    %%
    res_.dA = J_A_R;
    res_.dB = J_B_R;
    res_.dC = J_C_R;
    res_.dD = J_D_R;
    
    res_.A = A;
    res_.B = B;
    res_.C = C;
    res_.D = D;
    

end
    

function S1 = s1matrix(q,r,u,nmax)

    I_1 = speye(q);
    I_2 = speye(q*r-nmax);
    S1 = [];

    for i = 1:u
        e_i = zeros(1,u);
        e_i(1,i) = 1;
        e_i = sparse(e_i);
        t_1 = kron(kron(I_1,e_i),I_2);
        S1 = [S1;t_1];            
    end

end

function S4 = s4matrix(q,r,nmax)

    % make a selection matrix S4     
    k_vec = linspace(q,1,q);
    S4 = [];
    for i = 1:q
        k = k_vec(i);
        S_3k = [sparse((q-k)*r,k*r); speye(k*r)];
        S_4k = sparse(q*r,(q-k)*r);
        S_2 = ([S_3k S_4k])';
        S4 = [S4;kron(speye(q*r-nmax),S_2)];
    end 

end   

function J_pinvLs_Ls = JJ_pinvLs_Ls(Ls,pinvLs)

    I_1 = speye(size(Ls,1));
    p1 = kron(pinvLs',-pinvLs);
    p2 = (pinvLs*pinvLs');
    p3 = kron((I_1 - Ls*pinvLs)',p2);

    Pab = permutation_matrix(size(Ls,1),size(Ls,2)); %correct

    J_pinvLs_Ls = p1 + p3*Pab;

end

function [J_O_H] = deltaO_H_clean(p,r,ntrue,Vs,Un,s1,covSSI)
    
    d_Uv_H = zeros((p+1)*r*ntrue,(p+1)*r*(p+1)*r);
    Unn = Un*Un';

    for i = 1:ntrue

        d_uf3 = kron((Vs(:,i)/s1(i))',Unn);
        d_Uv_H((p+1)*r*(i-1)+1:(p+1)*r*i,:) = d_uf3;

    end    

    J_O_H = d_Uv_H*covSSI;

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


    
    
    
    
    
    