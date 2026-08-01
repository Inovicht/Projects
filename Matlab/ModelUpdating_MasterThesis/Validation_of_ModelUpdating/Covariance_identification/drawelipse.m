function drawelipse(lambda,cov)

invcov = real(inv(cov));
r = real(eig(invcov));
x = linspace(real(lambda)-r(1)*sqrt(9.210),real(lambda)+r(1)*sqrt(9.210),1000);
y = linspace(imag(lambda)-r(2)*sqrt(9.210),imag(lambda)+r(2)*sqrt(9.210),1000);
[x,y] = meshgrid(x,y);
z = (x-real(lambda)).^2/(r(1)^2) + (y-imag(lambda)).^2/(r(2)^2);
contour(x,y,z,[0,9.210])
axis padded
% x = reshape(x.',1,[]);
% y = reshape(y.',1,[]);
% z = reshape(z.',1,[]);
% scatter3(x,y,z)