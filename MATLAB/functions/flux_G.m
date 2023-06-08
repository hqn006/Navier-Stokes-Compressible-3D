% Construct G vector (which has outer derivative in z)
%
% INPUTS
%  G : flux vector at previous step
%  Primitive variables : see prim2cons()
%  dx,dy,dz : intervals in space
%
% PARAMETERS
%  whichStep : 'Predictor' or 'Corrector', determines whether to use
%              forward or backward FDs, respectively
%  Adiabatic : (bool) zero heat flux at walls
%
% OUTPUT
%  G : flux vector at next step, used in dG/dz
%
%
% DEPENDENCIES
%  * FD approximation functions in the same folder


function [ G ] = flux_G( G, rho,u,v,w, T,p,Et, mu,k, dx,dy,dz, ...
                         whichStep, Adiabatic )

if     strcmp(whichStep,'Predictor')
    ddz = @ddz_bwd;
elseif strcmp(whichStep,'Corrector')
    ddz = @ddz_fwd;
end

% Precompute gradients
% Backward FDs for nested derivatives
dudz = ddz(u,dz);
dvdz = ddz(v,dz);
dwdz = ddz(w,dz);
dTdz = ddz(T,dz);

% Adiabatic walls, both boundaries of y,z
if (Adiabatic)
    dTdz(:,1,:) = 0;
    dTdz(:,end,:) = 0;
    dTdz(:,:,1) = 0;
    dTdz(:,:,end) = 0;
end

% Central FDs for derivatives without an outer 2nd derivative
dudx = ddx_central(u,dx);
dwdx = ddx_central(w,dx);
dvdy = ddy_central(v,dy);
dwdy = ddy_central(w,dy);

div_u = dudx + dvdy + dwdz; % divergence of velocity vector


% Compute mass, stress, and heat fluxes
rhow = rho.*w;
tau_zz = 2*mu.*( dwdz - div_u/3 );
tau_xz =   mu.*( dudz + dwdx );
tau_yz =   mu.*( dvdz + dwdy );
qdot_z = -k .* dTdz;


% Construct F vector
G(1,:,:,:) = rhow;
G(2,:,:,:) = rhow.*u     - tau_xz;
G(3,:,:,:) = rhow.*v     - tau_yz;
G(4,:,:,:) = rhow.*w + p - tau_zz;
G(5,:,:,:) = (Et+p).*w - u.*tau_xz - v.*tau_yz - w.*tau_zz + qdot_z;

end


%% EOF
