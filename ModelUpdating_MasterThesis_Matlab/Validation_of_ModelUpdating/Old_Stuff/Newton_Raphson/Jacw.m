function J = Jacw(params,sys,G,theta,delta)
    C = sys.Cc;
    B = sys.B;
    A = sys.A;
    n = length(theta);
    m = length(A);
    Ag_origin = A - B*G*C;
    [V,~,W] = eig(Ag_origin);
    W = W';
    J = zeros(m,n);
    for i = 1:n
        theta_p = theta;
        theta_p(i) = theta_p(i)+delta;
        sys_p = Sys(params,theta_p);
        Agp = sys_p.A + sys_p.B*G*sys_p.Cc;
        DeltaA = (Agp - Ag_origin)/delta;
        for j = 1:m
            J(j,i) = 1/(W(j,:)*V(:,j))*W(j,:)*DeltaA*V(:,j);
        end
    end
end
 