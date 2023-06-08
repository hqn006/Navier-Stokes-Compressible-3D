% Enforce Boundary Conditions on primitives
%
% No-slip, adiabatic at 4 walls in y,z
% Inlet  : Far-field conditions
% Outlet : Extrapolate


function [ rho,u,v,w, T,p,e,Et ] = BC( rho,u,v,w, T,p,e,Et, u_Inf, ...
                                       Adiabatic )


%% Update BCs that are explicitly defined
% @WALLS (both boundaries in y & z):
if (Adiabatic)
    % Zero heat flux at walls
    % From the 2nd-Order 1st-Difference Approximation set equal to 0
    % Equation: f'(y1) = 1/(2*h) * ( -3*f(y1) + 4*f(y2) - f(y3) ) = 0
    T(:,1,:) = (4*T(:,2,:) - T(:,3,:)) / 3;
    T(:,end,:) = (4*T(:,end-1,:) - T(:,end-2,:)) / 3;
    T(:,:,1) = (4*T(:,:,2) - T(:,:,3)) / 3;
    T(:,:,end) = (4*T(:,:,end-1) - T(:,:,end-2)) / 3;
else
    % Constant Temperature
    T(:,1,:) = const.T0;   % bottom
    T(:,end,:) = const.T0; % top
    T(:,:,1) = const.T0;   % front
    T(:,:,end) = const.T0; % back
end
u(:,1,:) = 0; v(:,1,:) = 0; w(:,1,:) = 0;
u(:,end,:) = 0; v(:,end,:) = 0; w(:,end,:) = 0;
u(:,:,1) = 0; v(:,:,1) = 0; w(:,:,1) = 0;
u(:,:,end) = 0; v(:,:,end) = 0; w(:,:,end) = 0;
% Extrapolation from 2 neighbor points
p(:,1,:) = 2*p(:,2,:) - p(:,3,:);
p(:,end,:) = 2*p(:,end-1,:) - p(:,end-2,:);
p(:,:,1) = 2*p(:,:,2) - p(:,:,3);
p(:,:,end) = 2*p(:,:,end-1) - p(:,:,end-2);


% @INLET:
u(1,:,:) = u_Inf; v(1,:,:) = 0; w(1,:,:) = 0;
p(1,:,:) = const.p0; T(1,:,:) = const.T0;

% @OUTLET:
% Extrapolation from 2 neighbor points
u(end,:,:) = 2*u(end-1,:,:) - u(end-2,:,:);
v(end,:,:) = 2*v(end-1,:,:) - v(end-2,:,:);
w(end,:,:) = 2*w(end-1,:,:) - w(end-2,:,:);
p(end,:,:) = 2*p(end-1,:,:) - p(end-2,:,:);
T(end,:,:) = 2*T(end-1,:,:) - T(end-2,:,:);


%% Update boundaries for dependent variables
% @WALLS
rho(:,1,:) = p(:,1,:)./T(:,1,:) /const.R;       % bottom
rho(:,end,:) = p(:,end,:)./T(:,end,:) /const.R; % top
rho(:,:,1) = p(:,:,1)./T(:,:,1) /const.R;       % front
rho(:,:,end) = p(:,:,end)./T(:,:,end) /const.R; % back

rho(1,:,:) = p(1,:,:)./T(1,:,:) /const.R;       % @INLET
rho(end,:,:) = p(end,:,:)./T(end,:,:) /const.R; % @OUTLET


% @WALLS
e(:,1,:) = const.cv.*T(:,1,:); % bottom
e(:,end,:) = const.cv.*T(:,end,:); % top
e(:,:,1) = const.cv.*T(:,:,1); % front
e(:,:,end) = const.cv.*T(:,:,end); % back

e(1,:,:) = const.cv.*T(1,:,:);     % @INLET
e(end,:,:) = const.cv.*T(end,:,:); % @OUTLET


% @WALLS
Et(:,1,:) = rho(:,1,:).*( e(:,1,:) + 0.5*(u(:,1,:).^2 + v(:,1,:).^2) );
Et(:,end,:) = rho(:,end,:).*( e(:,end,:) + 0.5*(u(:,end,:).^2 + v(:,end,:).^2) );
Et(:,:,1) = rho(:,:,1).*( e(:,:,1) + 0.5*(u(:,:,1).^2 + v(:,:,1).^2) );
Et(:,:,end) = rho(:,:,end).*( e(:,:,end) + 0.5*(u(:,:,end).^2 + v(:,:,end).^2) );

Et(1,:,:) = rho(1,:,:).*...
            ( e(1,:,:) + 0.5*(u(1,:,:).^2 + v(1,:,:).^2) ); % @INLET
Et(end,:,:) = rho(end,:,:).*...
              ( e(end,:,:) + 0.5*(u(end,:,:).^2 + v(end,:,:).^2) ); % @OUTLET


end


%% EOF
