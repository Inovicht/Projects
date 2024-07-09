function [C_c,C_a] = Cmat(M,K,D,d,d_dot,d_dotdot)
n = length(M);
m = length([d,d_dot,d_dotdot]);
k = 1;
C_d = zeros(m,n);
C_v = zeros(m,n);
C_a = zeros(m,n);

for i = 1:length(d)
    C_d(k,d(i)) = 1;
    k = k + 1;
end

for i = 1:length(d_dot)
    C_v(k,d_dot(i)) = 1;
    k = k + 1;
end

for i = 1:length(d_dotdot)
    C_a(k,d_dotdot(i)) = 1;
    k = k + 1;
end

C_c = [C_d-C_a/M*K, C_v - C_a/M*D];
end