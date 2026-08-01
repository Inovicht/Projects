function [B,B2] = Bmat(M,b)

n = length(M);
r = length(b);
B2 = zeros(n,r);
k = 1;

for i = 1:length(b)
    B2(b(i),k) = 1;
    k = k + 1;
end

B = [zeros(n,r);inv(M)*B2];
end
