function stop = outfun(theta,optimValues,state,params,cov) 
stop = false;

sys1 = Sys(params,theta);
h = constraintnew(params,theta,cov);
if ~any(h>0) 
stop = true; 
disp('Stopping,');
end