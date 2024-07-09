function drawelipse(MS_data,params,lambdapl,lambda_new)
% This function draws the uncertainty ellipse for the covariance. Furthermore plotes the estimated
% eigenvalues, pertubated eigenvalues ,and the exact eigenvalues

%% Defining the quantile for ellipse
p=0.99;
s=-2*log(1-p);
%% Extracting Monte Carlo estimates
MS_length = size(MS_data,2);
MS_real = MS_data(:,1:MS_length/2);
MS_imag = MS_data(:,MS_length/2+1:end);

% Points from the exact and pertubated eigenvalues
lambda = params.lambdatilt; % Pertubated eigenvalues
lambda_estimate = lambda_new; % Exact eigenvalues
lambda_exact = lambdapl;    % Exact eigenvalues
ewq = 1; % Variables for defining name in title
figure(1)
% Iterating through all the placed eigenvalues
for i = 1:length(params.cov)
%% Making Data
% Parameters for uncertainty ellipse
[V,D]=eig(abs(params.cov{i})*s);
ttt=linspace(0,2*pi);
a=(V*sqrt(D))*[cos(ttt(:))'; sin(ttt(:))'];

% Restructuring the eigenvalues for (real,imag) plot
mu=[real(lambda(i)); imag(lambda(i))];                            % Placement of pertubated eigenvalues
mu_estimate=[real(lambda_estimate(i)); imag(lambda_estimate(i))]; % Placement of estimated eigenvalue
mu_exact=[real(lambda_exact(i)); imag(lambda_exact(i))];          % Placement of exact eigenvalue

% Creating subplots for each placed eigenvalue
subplot(3,3,i)
hold on
grid on

% Plotting the eigenvalues from Monte Carlo, exact, pertubated and
% estimated
scatter(MS_real(:,i),MS_imag(:,i),'k');                             % Placing all the Monte Carlo estimates
plot(mu(1),mu(2),'dr',linewidth=5, markersize=5);                   % Pertubated eigenvalue
plot(mu_estimate(1),mu_estimate(2),'db',linewidth=5, markersize=5); % Estimated eigenvalue
plot(mu_exact(1),mu_exact(2),'b',linewidth=5, markersize=5);        % Exact eigenvalues
% Plotting the uncertainty ellipse
plot(a(1, :) + mu(1), a(2, :) + mu(2),'k--');                       % Plotting Ellipse

%% Defining subplot title and labeling
set(get(gcf,'CurrentAxes'),'FontName','Times New Roman','FontSize',12) % Settings for the plot
uf = mod(i,1+params.n_CL);
if uf == 1 % Open-Loop
title('\textbf{}' + string(ewq) + '\textbf{ Eigenvalue - Open-loop}','FontName','Times New Roman','FontSize',12,'Interpreter','latex')
elseif uf == 2
title('\textbf{}' + string(ewq) + '\textbf{ Eigenvalue - $1$st closed-loop}','FontName','Times New Roman','FontSize',12,'Interpreter','latex')
else
title("    " + string(ewq) + '\textbf{ Eigenvalue - $2$nd closed-loop}','FontName','Times New Roman','FontSize',12,'Interpreter','latex')
ewq = ewq+1;
end
xlabel('$\Re(\hat{\lambda})$','FontName','Times New Roman','FontSize',12,'Interpreter','latex')
ylabel('$\Im(\hat{\lambda})$','FontName','Times New Roman','FontSize',12,'Interpreter','latex')
hold off
end
legend('Delta method-based $99~\%$ confidence ellipse','Sys ID eigenvalue','model updated eigenvalue','MS without model increase','MS with increase','FontName','Times New Roman','FontSize',10,'Interpreter','latex','location','best')