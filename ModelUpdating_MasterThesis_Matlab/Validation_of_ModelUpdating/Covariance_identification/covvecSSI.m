function [covSSI,Htot,Ctot,hk]=covvecSSI(inMat,param)

% methd : Choice of methods 
% p, q : integers, number of blocks rows/columns for inMat
% inMat : data samples of size(s,r)
% r : number of sensors
% N : number of samples
% nb : number of blocks for Y+ and Y-
% r0  : number of chosen sensors
% pq = p + q 

methd=param.mthd;
p=param.p;
q=param.q;
pchannel=param.pchannel;
nb=param.nb;
Htot=0;

r0=size(pchannel,2);
r=size(inMat,2);

if isfield(param,'inputs')
    inp=param.inputs;
    r=r-inp;
end

N = size(inMat,1);
Nn = nb*floor((N-p-q)/nb);
Nb = Nn/nb;
pq=p+q;
covSSI=zeros((p+1)*r*q*r0,nb);
    
if 1==strcmp('IOcov_all',methd)
    [hk,Htot] = SSImat(inMat,param);
    YfYp=zeros((p+1)*r*q*r0,nb);
    YfUf=zeros((p+1)*r*(p+1)*inp,nb);
    UfUf=zeros((p+1)*inp*(p+1)*inp,nb);
    UfYp=zeros((p+1)*inp*q*r0,nb);
    
    for i = 1:nb
        [~,H] = SSImat(inMat(((i-1)*Nb+1):i*Nb+pq,:),param);
        YfYp(:,i)= H.YfYp(:)/sqrt(nb*(nb-1));
        YfUf(:,i)= H.YfUf(:)/sqrt(nb*(nb-1));
        UfUf(:,i)= H.UfUf(:)/sqrt(nb*(nb-1));
        UfYp(:,i)= H.UfYp(:)/sqrt(nb*(nb-1));
    end
    %substract mean
    M1=mean(YfYp,2);
    M2=mean(YfUf,2);
    M3=mean(UfUf,2);
    M4=mean(UfYp,2);
    for i = 1:nb
        YfYp(:,i) = YfYp(:,i) - M1;
        YfUf(:,i) = YfUf(:,i) - M2;
        UfUf(:,i) = UfUf(:,i) - M3;
        UfYp(:,i) = UfYp(:,i) - M4;
    end
    iUfUf=inv(Htot.UfUf);
    covSSI = YfYp ...
      - kron(Htot.UfYp'*iUfUf', speye((p+1)*r)) * YfUf ...
      + kron(Htot.UfYp'*iUfUf', Htot.YfUf*iUfUf) * UfUf ...
      - kron(speye(q*r0), Htot.YfUf*iUfUf) * UfYp;

    Ctot.YfYp = YfYp;
    Ctot.YfUf = YfUf;
    Ctot.UfUf = UfUf;
    Ctot.UfYp = UfYp;    
    
end



end