function [cov_dAd] = Confidensinterval(res_id,G,fs)
%This function computes the closed-loop eigenvalues and their corresponding
%confidence matrices for one closed-loop realization thus A_g = A-B*G*C

%Input: 
%res_id: Structure that contain the system quadruplet in discrete time
%domain and jacobians
%Gain matrix (Note this is one Gain-realization, not the augmented gain matrix!)
%Sample frequency fs

%Outputs:
%A structure cov_dAd that contain a vector that contain the closed-loop
%eigenvalues
% and a cell that contain the confidence matrix of each CL eigevalues

param.mthd='IOcov_all'; % I/O covariance-driven identification
param.nmax=length(res_id.A);

dt = 1/fs;
Ad_ID=res_id.A;
Bd_ID = res_id.B;
Cd_ID = res_id.C;
Dd_ID = res_id.D;
dAd_ID=res_id.dA;
dBd_ID = res_id.dB;
dCd_ID = res_id.dC;
dDd_ID = res_id.dD;
A_ID=logm(Ad_ID)/dt;
B_ID=A_ID*(Ad_ID-eye(size(Ad_ID)))^-1*Bd_ID;
C_ID=Cd_ID;
D_ID=Dd_ID;
% Eigenvalue uncertainty
cov_dAd=coveigAbcont(Ad_ID,Bd_ID,Cd_ID,G,dAd_ID,dBd_ID,dCd_ID,param.nmax,dt);   
vec=cov_dAd.r_vec;
poles=cov_dAd.e_vald;
cov_poles=cov_dAd.cov_lambda_Ab; 

%% Small function(s)
function KL=gain_design(A,B,C,poles)
    r=size(B,2);
    m=size(C,1);
    N=size(A,1);
    P1=[];
    P2=[];
    for j=1:r
        Z=A.'-poles(j)*eye(N);
        ZZ=[Z -C.']; % Wide matrix, so null command is OK
        V=null(ZZ);
        F=V*randn(m,1);       
        FF=F(:,1);  
        part1=FF(1:N,:);  
        part2=FF(N+1:N+m,:);
        P1=[P1 part1];
        P2=[P2 part2];
    end
    CP1=B.'*P1;
    KL=pinv(CP1.')*P2.';
end
end







    