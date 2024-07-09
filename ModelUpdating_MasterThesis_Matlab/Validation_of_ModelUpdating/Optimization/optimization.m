function theta_new = optimization(params,theta)
%This function computes new values of the boundary springs and axial
%tension, such that the discrepancy between the target eigenvalues and the
%preditions of the model is minimized. 

%inputs:
%structure params (Check params.m for content and description)
%vector that contains the values of the boundary springs and axial tension
%[k1,k2,k3,k4,N]

%Outputs: 
%theta_new: model updated boundary springs and bolt tension

%% Setting up const functions, non linear constraints, and stopping criteria
f = @(x) costfunc(params,x);                    %Defining the cost function
nonlcon = @(x) constraint(params,x);            %Defining the inequality constraints of the optimization problem
%nonlcon = [];                                  
stop = @(x,optim,state) StopCriteria(x,optim,state,params); %A stop condition that if all model predictions are withing the constraints, then the optimization stops

%% Generating start guess:
radius = 0.2*theta;
select_Neighbor = @(x)([x(1)+radius(1)*2*(rand(1)-0.5),x(2)+radius(2)*2*(rand(1)-0.5),x(3)+radius(1)*2*(rand(1)-0.5),x(4)+radius(2)*2*(rand(1)-0.5),x(5)+radius(3)*2*(rand(1)-0.5)]');
x0 = select_Neighbor(theta);

%% Setting up the options of the optimization
options=optimoptions('fmincon','Algorithm','sqp','OutputFcn',stop,'Display','iter-detailed','SpecifyObjectiveGradient',false,'SpecifyConstraintGradient',false);
options.MaxFunEvals = 5e3;
options.StepTolerance = 1e-10;
options.OptimalityTolerance = 1e-10;
options.ConstraintTolerance = 1e-10;
options.UseParallel = true;
%options.EnableFeasibilityMode = true;
%options.SubproblemAlgorithm = 'cg';

theta_new = fmincon(f,x0,params.Ab,params.b,params.Aeq,params.beq,params.lb,params.ub,nonlcon,options);         %We use the fmincon function the run the optimization problem. Specifically, we use sequential quadratic programming.
end