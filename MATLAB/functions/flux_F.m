% Construct F vector (which has outer derivative in y)
%
% INPUTS
%  F : flux vector at previous step
%  Primitive variables : see prim2cons()
%  dx,dy,dz : intervals in space
%
% PARAMETERS
%  whichStep : 'Predictor' or 'Corrector', determines whether to use
%              forward or backward FDs, respectively
%  Adiabatic : (bool) zero heat flux at walls
%
% OUTPUT
%  F : flux vector at next step, used in dF/dy
%
%
% DEPENDENCIES
%  * FD approximation functions in the same folder


function [ F ] = flux_F( F, rho,u,v,w, T,p,Et, mu,k, dx,dy,dz, ...
                         whichStep, Adiabatic )

if     strcmp(whichStep,'Predictor')
    ddy = @ddy_bwd;
elseif strcmp(whichStep,'Corrector')
    ddy = @ddy_fwd;
end

% Precompute gradients
% Backward FDs for nested derivatives
dudy = ddy(u,dy);
dvdy = ddy(v,dy);
dwdy = ddy(w,dy);
dTdy = ddy(T,dy);

% Adiabatic walls, both boundaries of y,z
if (Adiabatic)
    dTdy(:,1,:) = 0;
    dTdy(:,end,:) = 0;
    dTdy(:,:,1) = 0;
    dTdy(:,:,end) = 0;
end

% Central FDs for derivatives without an outer 2nd derivative
dudx = ddx_central(u,dx);
dvdx = ddx_central(v,dx);
dvdz = ddz_central(v,dz);
dwdz = ddz_central(w,dz);

div_u = dudx + dvdy + dwdz; % divergence of velocity vector


% Compute mass, stress, and heat fluxes
rhov = rho.*v;
tau_yy = 2*mu.*( dvdy - div_u/3 );
tau_xy =   mu.*( dudy + dvdx );
tau_yz =   mu.*( dvdz + dwdy );
qdot_y = -k .* dTdy;


% Construct F vector
F(1,:,:,:) = rhov;
F(2,:,:,:) = rhov.*u     - tau_xy;
F(3,:,:,:) = rhov.*v + p - tau_yy;
F(4,:,:,:) = rhov.*w     - tau_yz;
F(5,:,:,:) = (Et+p).*v - u.*tau_xy - v.*tau_yy - w.*tau_yz + qdot_y;

end


%% EOF
