function stop = StopCriteria(theta,optimValues,state,params) 
%This function serves as a stopping criteria for the optimization problem.
%If all inequality constraints are satisfied, then the optimzation problem
%halts.

%inputs:
%structure params (Check params.m for content and description)
%vector that contains the values of the boundary springs and axial tension
%[k1,k2,k3,k4,N]
%optimValues: We dont use it, but is nesessary to have as input to run
%fmincon
%state: We dont use it, but is nesessary to have as input to run
%fmincon

stop = false;
h = constraint(params,theta);

if ~any(h>0) 
stop = true; 
disp('Stopping,');
end