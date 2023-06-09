% Compute dynamic viscosity of air using Sutherland's law.
%
% mu0 and T0 changed for the square duct problem. mu0 is based on the given
% Re_D value.
%
% A point of error may be that the meaning of S1 is unknown, so it is left
% unchanged.
% 
% INPUTS
% T   : Temperature
% mu0 : reference viscosity, taken at inlet based on Re_D
% 
% OUTPUTS
% mu : dynamic viscosity


function [ mu ] = sutherland( T, mu0 )

% Sutherland's law
mu = mu0*(T/const.T0).^(3/2) .* ( (const.T0+const.S1) ./ (T+const.S1) );

end


%% EOF
