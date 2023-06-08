% Compute numerical schlieren image


function [ S ] = schlieren( rho, dx,dy,dz )

% Constant fudge factors
Beta = 0.8;
kappa = 10;


% Compute magnitude of the density gradient
drhodx = ddx_central(rho,dx);
drhody = ddy_central(rho,dy);
drhodz = ddz_central(rho,dz);

grad_rho = sqrt( drhodx.^2 + drhody.^2 +drhodz.^2 ); % magnitude


S = Beta * exp( -kappa * grad_rho/max(grad_rho(:)) );

end


%% EOF
