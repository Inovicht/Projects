function K = Kmat(K_origin,Kg,theta)
    a = length(theta);
    K = K_origin;
    K(1,1) = K(1,1) + theta(1);
    K(2,2) = K(2,2) + theta(2);
    K(end-1,end-1) = K(end-1,end-1) + theta(3);
    K(end,end) = K(end,end) + theta(4);
    if length(theta) > 4
        K = K + theta(5)*Kg;
    end
end