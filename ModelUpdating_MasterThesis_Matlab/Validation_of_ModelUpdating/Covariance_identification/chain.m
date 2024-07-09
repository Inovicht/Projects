% *************** CHAIN.m *********************
% Program to generate system matrices for a chain-like system.


function [M,C,K]=chain(m1,c1,k1)

dof=length(m1);
if dof~=1
    kd=zeros(dof,1);    % elements on diagonal
    cd=zeros(dof,1);
    kf=zeros(dof-1,1);  % elements on off-diagonal
    cf=zeros(dof-1,1);
    
    for j=1:dof-1
       kd(j)=k1(j)+k1(j+1);
       cd(j)=c1(j)+c1(j+1);
       kf(j)=-k1(j+1);
       cf(j)=-c1(j+1);
    end;
    kd(dof)=k1(dof);
    cd(dof)=c1(dof);

    % Create the system matrices
    K= diag(kf,1)+diag(kf,-1)+diag(kd);
    C=(diag(cf,1)+diag(cf,-1)+diag(cd));
    M=diag(m1);
else
    K=k1;
    M=m1;
    C=c1;
end;