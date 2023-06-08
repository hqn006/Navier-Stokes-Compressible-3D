% Construct E vector (which has outer derivative in x)
%
% INPUTS
%  E : flux vector at previous step
%  Primitive variables : see prim2cons()
%  dx,dy,dz : intervals in space
%
% PARAMETERS
%  whichStep : 'Predictor' or 'Corrector', determines whether to use
%              forward or backward FDs, respectively
%  Adiabatic : (bool) zero heat flux at walls
%
% OUTPUT
%  E : flux vector at next step, used in dE/dx
%
%
% DEPENDENCIES
%  * FD approximation functions in the same folder


function [ E ] = flux_E( E, rho,u,v,w, T,p,Et, mu,k, dx,dy,dz, ...
                         whichStep, Adiabatic )

if     strcmp(whichStep,'Predictor')
    ddx = @ddx_bwd;
elseif strcmp(whichStep,'Corrector')
    ddx = @ddx_fwd;
end

% Precompute gradients
% Backward FDs for nested derivatives
dudx = ddx(u,dx);
dvdx = ddx(v,dx);
dwdx = ddx(w,dx);
dTdx = ddx(T,dx);

% Adiabatic walls, both boundaries of y,z
if (Adiabatic)
    dTdx(:,1,:) = 0;
    dTdx(:,end,:) = 0;
    dTdx(:,:,1) = 0;
    dTdx(:,:,end) = 0;
end

% Central FDs for derivatives without an outer 2nd derivative
dudy = ddy_central(u,dy);
dvdy = ddy_central(v,dy);
dudz = ddz_central(u,dz);
dwdz = ddz_central(w,dz);

div_u = dudx + dvdy + dwdz; % divergence of velocity vector


% Compute mass, stress, and heat fluxes
rhou = rho.*u;
tau_xx = 2*mu.*( dudx - div_u/3 );
tau_xy =   mu.*( dudy + dvdx );
tau_xz =   mu.*( dudz + dwdx );
qdot_x = -k .* dTdx;


% Construct E vector
E(1,:,:,:) = rhou;
E(2,:,:,:) = rhou.*u + p - tau_xx;
E(3,:,:,:) = rhou.*v     - tau_xy;
E(4,:,:,:) = rhou.*w     - tau_xz;
E(5,:,:,:) = (Et+p).*u - u.*tau_xx - v.*tau_xy - w.*tau_xz + qdot_x;

end


%% EOF
