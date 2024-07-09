function sys = Truncation(sys,r)
B = sys.B;
C_c = sys.Cc;
Ac = sys.A;

%% We truncate nominal state space matrix
n_trunk = r;
[V, d] = eig(Ac);
[d,index1] = sort(diag(d));
V = V(:,index1);
P = inv(V);

% Trunkating the eigenspectrum and the eigenvalues
d=d(1:n_trunk);
V_trunk = V(:,1:n_trunk);
P_trunk=P(1:n_trunk,:);

%Trunkating
sys.A = diag(d);
sys.B = P_trunk*B;
sys.Cc = C_c*V_trunk;
sys.lambda = eig(sys.A);