function [Hmat] = Hmat(Ca,M,B2)
%Function that computes the Transmission matrix H

%Inputs:
% Global mass matrix M
% input distribution matrix B2 in second order formulation
% Output distribution matrix for acceleration outputs Ca

%Outputs:
%Transmission matrix H

    Hmat = Ca/M*B2;
end
