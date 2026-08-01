function [A_ID,B_ID,C_ID,D_ID] = n4sid(uk,yk,ns,params)
% Ask your advisor for further explanation of n4sid

% We have inserted params in this function
   % params.sysidNhankel

%% Dimension check of captured measurements
[nu,Nu]=size(uk);
if Nu<nu
    uk=uk';
    [nu,Nu]=size(uk);
end
[ny,Ny]=size(yk);
if Ny < ny
    yk=yk';
    [ny,Ny]=size(yk);
end
if Nu~=Ny
    disp('Discrepancy in input and output sample lengths')
    return
end

% ii=ns*2; % User-defined number of block rows
ii = params.sysidNhankel;
s=Nu;
m=size(uk,1);
l=size(yk,1);
j=s-2*ii+1; % Number of columns in block Hankel matrix

%% Block Hankel matrix
uk1=uk(:,ii+1:s)/sqrt(j);
yk1=yk(:,ii+1:s)/sqrt(j);
uk=uk/sqrt(j);
yk=yk/sqrt(j);


Up = zeros(ii*nu,j);
Yp = zeros(ii*ny,j);
Uf = zeros(ii*nu,j);
Yf = zeros(ii*ny,j);
for k=1:ii 
    Up((k-1)*nu+1:k*nu,:)=uk(:,k:k+j-1);
    Yp((k-1)*ny+1:k*ny,:)=yk(:,k:k+j-1);
    Uf((k-1)*nu+1:k*nu,:)=uk1(:,k:k+j-1);
    Yf((k-1)*ny+1:k*ny,:)=yk1(:,k:k+j-1);
end

%% LQ Decomposition
U = [Up;Uf];
Y = [Yp;Yf];
L = triu(qr([U;Y]'))';
L = L(1:2*ii*(m+l),1:2*ii*(m+l));

%% Oblique projection
Lf = L((2*m+l)*ii+1:2*(m+l)*ii,:);              
Lp = [L(1:m*ii,:);L(2*m*ii+1:(2*m+l)*ii,:)];    
Lu=L(m*ii+1:2*m*ii,:);                         
Lfp = Lf*(eye(2*(m+l)*ii)-Lu'*pinv(Lu*Lu')*Lu);
Lpp = Lp*(eye(2*(m+l)*ii)-Lu'*pinv(Lu*Lu')*Lu);

Oi = Lfp*pinv(Lpp)*Lp;                          

%% SVD of Oblique projection and selection of model order
[U1,S,V1] = svd(Oi);
% ss=diag(S);
% figure;
% bar((1:l*ii), ss./max(ss), 'k');
% set(gca,'YScale','log')
% title('Singular Values');
% xlabel('Order');
% ylabel('Singular value (normalized)');
% zoom('on');
% YScale = 'log';
%% Retrieve the system matrices
U2=U1(:,1:ns);
S1=S(1:ns,1:ns);

gammai=U2*sqrt(S1);
gammaip1=gammai(1:l*(ii-1),:);

Lf=L((2*m+l)*ii+l+1:2*(m+l)*ii,:);              
Lp=[L(1:m*(ii+1),:);L(2*m*ii+1:(2*m+l)*ii+l,:)]; 
Lu=L(m*ii+m+1:2*m*ii,:);                        
Lfp=Lf*(eye(2*(m+l)*ii)-Lu'*pinv(Lu*Lu')*Lu);   
Lpp=Lp*(eye(2*(m+l)*ii)-Lu'*pinv(Lu*Lu')*Lu);    

Oip=Lfp*pinv(Lpp)*Lp;                            

%% Determine the state sequeces
Xi=pinv(gammai)*Oi;
Xip=pinv(gammaip1)*Oip;

%% Solve linear equation for system matrices
Left=[Xip ; L((2*m+l)*ii+1:(2*m+l)*ii+l,:)];
Right=[Xi ; L(m*ii+1:m*(ii+1),:)]; 
sol=Left/Right;

%% Retrieve system matrices from solution
A_ID=sol(1:ns,1:ns);
B_ID=sol(1:ns,ns+1:ns+m);
C_ID=sol(ns+1:ns+l,1:ns);
D_ID=sol(ns+1:ns+l,ns+1:ns+m);
    

