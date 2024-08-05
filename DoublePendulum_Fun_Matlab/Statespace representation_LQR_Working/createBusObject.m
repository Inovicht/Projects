
createBusObject();

function createBusObject()
    % Create a temporary instance of the params structure
    tempParams.m = [1, 1];
    tempParams.L = [1; 1];
    tempParams.g = 9.81;
    tempParams.q = [1, 2];
    tempParams.epsilon = 1e-5;
    tempParams.Q0 = diag([10, 10, 1, 1]);
    tempParams.R0 = 1;

    % Create a bus object
    paramsBusInfo = Simulink.Bus.createObject(tempParams);
    
    % Assign the name to the bus object
    paramsBus = eval(paramsBusInfo.busName);
    paramsBus.description = 'Bus object for params structure';

    % Save the bus object to the base workspace
    assignin('base', 'paramsBus', paramsBus);
end