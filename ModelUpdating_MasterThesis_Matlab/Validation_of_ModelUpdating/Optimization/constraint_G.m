function [h,heq,dh,dheq] = constraint_G(params,theta)
sys_opd = Sys(params,theta);
p=0.99;
s=-2*log(1-p);
lambdatilt = params.lambdatilt;
h = zeros(length(lambdatilt),1);
lambda_temp = findl(lambdatilt,sys_opd,params);
dh = zeros(length(lambdatilt),length(theta));
dh = zeros(length(theta),length(lambdatilt));
for i = 1:length(lambdatilt)
    % Get the difference
    lambda_nom = [real(lambdatilt(i));imag(lambdatilt(i))];
    lambda_tilt =[real(lambda_temp(i));imag(lambda_temp(i))];
    diff = lambda_tilt-lambda_nom;
    % Get the constraint
    h(i) = real(diff'/params.cov{i}*diff) - s;
    % Get the gradient of the constraint
    Ag = sys_opd.A - sys_opd.B*params.G(i,:)*sys_opd.Cc;
    Ja = Jacw(params,sys_opd,params.G(i,:),theta,0.1);
    l = eig(Ag);
    indx = indxfind(l,lambdatilt(i));
    J = Ja(indx,:);
    Jr = real(J);
    Ji = imag(J);
    J = [Jr;Ji];
    for j = 1:length(theta)
        dh(j,i) = diff.'/real(params.cov{i})*J(:,j);
    end
end
h;
dh;
heq = [];
dheq = [];