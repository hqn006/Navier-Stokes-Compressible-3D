% Properties of air at standard (sea-level) conditions (default)
%
% - p0, T0 changed to IC of square duct problem
% - rho0 updated accordingly
% - mu0 updated based on Re_D


classdef const
    properties (Constant)
    p0 = 416 * 10^3;  % pressure, [Pa]
    T0 = 300;  % Temperature,     [K]
    
    R = 287;   % ideal gas constant,                 [J/kg/K]
    cp = 1005; % specific heat at constant pressure, [J/kg/K]
    cv = 718;  % specific heat at constant volume,   [J/kg/K]
    gamma = 1.4; % ratio of specific heats,          [1]

    rho0 = const.p0/const.R/const.T0; % density, [kg/m^3]

    S1 = 110.4; % Sutherland Temperature, [K]
    Pr = 0.71; % Prandtl number,          [1]
    end
end


%% EOF
