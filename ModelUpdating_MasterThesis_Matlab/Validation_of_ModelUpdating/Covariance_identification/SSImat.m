function [hk,R11]=SSImat_Szymon(inMat,param)
%==========================================================================
%  Computation of the Hankel matrix for IO covariance-driven identification 
%==========================================================================
% N : numner of samples
% p : number of block rows
% q : number of block columns
% lag : offset of the correlations 
% hk : Hankel matrix (p+1,q) block of size (ny,npc)
% r0 : number of chosen sensors
% r : number of sensors
% inMat : data samples of size(N,r)

mthd=param.mthd;
p=param.p;
q=param.q;
pchannel=param.pchannel;

if isfield(param,'inputs')
    inp=param.inputs;
end

r0=length(pchannel);
r=size(inMat,2);
N=size(inMat,1); 


if 1==strcmp('IOcov_all',mthd) 
    %cov
    r=r-inp;
    pq = p + q;
    NN = N-pq;
    
    Y=inMat(:,1:r);
    U=inMat(:,r+1:end);
    
    Uf = zeros((p+1)*inp,NN);
    Up = zeros(q*inp,NN);
    Yf = zeros((p+1)*r,NN);
    Yp = zeros(q*r0,NN);
    
    for i=0:q-1
        Up(i*inp + (1:inp), :) = U(i+(1:NN),:)';
        Yp(i*r0 + (1:r0), :) = Y(i+(1:NN),pchannel)';
    end
    for i=0:p
        Uf(i*inp + (1:inp), :) = U(q+i+(1:NN),:)';
        Yf(i*r + (1:r), :) = Y(q+i+(1:NN),:)';
    end
%     hk=Yf*Yp'-Yf*Uf'/(Uf*Uf')*Uf*Yp';
    R11.YfYp=Yf*Yp'/NN;
    R11.YfUf=Yf*Uf'/NN;
    R11.UfUf=Uf*Uf'/NN;
    R11.UfYp=Uf*Yp'/NN;
    hk=R11.YfYp-R11.YfUf/R11.UfUf*R11.UfYp;
end


end



