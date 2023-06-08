% MAE 290C
% Project : MacCormack for compressible Navier-Stokes
%
% - Conservative vector (flux) form of 3D compressible Navier-Stokes
% - Assume calorically perfect gas
%
%
% Validation case:
% Yang, R., Wang, Z., Zhao, Y., Wang, Q., amp; Feng, W. (2020). Numerical 
% investigation on spatial development of the secondary flow in a supersonic 
% turbulent square duct. Aerospace Science and Technology, 100, 105832. 
% https://doi.org/10.1016/j.ast.2020.105832
% 
%
% Huy Nguyen, Blake Carter
% Created: 5 June 2023


clearvars
close all
clc

addpath('classes')
addpath('functions')

format long
set(0, 'DefaultTextInterpreter','Latex')




%% Input

M = 3.9; % Inflow Mach number = u_Inf/a, [1]


% Number of grid points
Nx = 30;
Ny = 20;
Nz = 20;

% Domain (x,y,z : L x H x W)
% 1 in. x 1 in. (1 in. = 25.4 mm)
H = 25.4 * 10^-3; % Height, [m]
W = 25.4 * 10^-3; % Width,  [m]
L = 6 * H; % Length, [m]


% Time span
% dt = 2.35 * 10^-11; % Time step, [s]
dt = 1 * 10^-7;
final_time = 1 * 10^-3;

% Number of iterations
max_iter = floor(final_time / dt)


% Convergence variable
converge_name = 'u';

% Update every _ iterations
% update_rate = ceil(max_iter/10); % variable field plots
update_rate = 100;
% update_conv = update_rate/5;     % convergence plot
update_conv = 1;
print_rate = 1;
video_rate = update_rate;
% The convergence plot will be updated to a smooth plot at the end


%% Enable/Disable

% Numerical schlieren ('streaks')
useSchlieren = true; % Replace density with numerical schlieren

% Adiabatic wall
Adiabatic = true; % Zero heat flux through wall BC




%% Setup

% Grid spacing
dx = L / (Nx-1);
dy = H / (Ny-1);
dz = H / (Ny-1);
x = 0 : dx : L;
y = 0 : dy : H;
z = 0 : dz : W;
[xx,yy,zz] = ndgrid(x,y,z);


% Properties @Infinity
a_Inf = sqrt(const.gamma * const.R * const.T0); % speed of sound
u_Inf = M*a_Inf;                                % from input Mach number

% Initial Condition @(x,y, t=0)
rho = ones(Nx,Ny,Nz) * const.rho0;
% p   = ones(Nx,Ny,Nz) * const.p0; % Unused, calculated later
T   = ones(Nx,Ny,Nz) * const.T0;

u = ones(Nx,Ny,Nz) * u_Inf;
v = zeros(Nx,Ny,Nz);
w = zeros(Nx,Ny,Nz);


% Preallocate vector form of Navier-Stokes
U = prim2cons( rho, u,v,w, T, const.cv );
E    = zeros(size(U));
F    = zeros(size(U));
G    = zeros(size(U));
Ubar = zeros(size(U));
Ebar = zeros(size(U));
Fbar = zeros(size(U));
Gbar = zeros(size(U));

% Preallocate derivatives of flux vectors
dEdx = zeros(size(U));
dFdy = zeros(size(U));
dGdz = zeros(size(U));
dEbardx = zeros(size(U));
dFbardy = zeros(size(U));
dGbardz = zeros(size(U));


% Set up vector to check convergence over time
t = 0 : dt : (max_iter-1)*dt;
converge = zeros(max_iter,1);

% Example string: "max(u,[],'all')"
converge_str = strcat("max(", converge_name, ",[],'all')");
title_str = ['Convergence of max ', converge_name];

converge(1) = eval( converge_str );


% Compute primitive variables
[ rho,u,v,w,T,p,e,Et ] = cons2prim( U, const.R,const.cv );

% Enforce Boundary Conditions on primitives
[ rho,u,v,w, T,p,~,Et ] = BC( rho,u,v,w, T,p,e,Et, u_Inf, Adiabatic );


% Write video
video = VideoWriter('nse-3d.mp4', 'MPEG-4');
open(video);




%% Solve
% MacCormack Method for compressible Navier-Stokes

% Compute numerical schlieren image
if (useSchlieren)
    S = schlieren(rho, dx,dy,dz);
else
    S = [];
end

% Plot Initial Condition
fig = figure;
plot_fields(rho,u,v, T,p,e, xx,yy, useSchlieren,S)


% Loop over time
timekeep_calc = zeros(max_iter,1);
timekeep_plot = zeros(max_iter,1);
time_all = tic;
iteration = 0; % number of completed iterations
while iteration < max_iter
    time = tic;


    if ~( isreal(rho) )
        error('Imaginary density found.')
    end

    
    %% PREDICTOR STEP
    % Forward FDs for outer derivatives
    % Backward FDs for inner nested derivatives

    % Update physical parameters
    mu = sutherland(T);         % dynamic viscosity from Sutherland's law
    k = const.cp/const.Pr * mu; % thermal conductivity


    % Construct flux vectors
    E = flux_E(E, rho,u,v,w, T,p,Et, mu,k, dx,dy,dz, 'Predictor', Adiabatic);
    F = flux_F(F, rho,u,v,w, T,p,Et, mu,k, dx,dy,dz, 'Predictor', Adiabatic);
    G = flux_G(G, rho,u,v,w, T,p,Et, mu,k, dx,dy,dz, 'Predictor', Adiabatic);


    % Precompute flux vector derivatives
    dEdx(1,:,:,:) = ddx_fwd(squeeze(E(1,:,:,:)),dx);
    dEdx(2,:,:,:) = ddx_fwd(squeeze(E(2,:,:,:)),dx);
    dEdx(3,:,:,:) = ddx_fwd(squeeze(E(3,:,:,:)),dx);
    dEdx(4,:,:,:) = ddx_fwd(squeeze(E(4,:,:,:)),dx);
    dEdx(5,:,:,:) = ddx_fwd(squeeze(E(5,:,:,:)),dx);

    dFdy(1,:,:,:) = ddy_fwd(squeeze(F(1,:,:,:)),dy);
    dFdy(2,:,:,:) = ddy_fwd(squeeze(F(2,:,:,:)),dy);
    dFdy(3,:,:,:) = ddy_fwd(squeeze(F(3,:,:,:)),dy);
    dFdy(4,:,:,:) = ddy_fwd(squeeze(F(4,:,:,:)),dy);
    dFdy(5,:,:,:) = ddy_fwd(squeeze(F(5,:,:,:)),dy);

    dGdz(1,:,:,:) = ddz_fwd(squeeze(G(1,:,:,:)),dz);
    dGdz(2,:,:,:) = ddz_fwd(squeeze(G(2,:,:,:)),dz);
    dGdz(3,:,:,:) = ddz_fwd(squeeze(G(3,:,:,:)),dz);
    dGdz(4,:,:,:) = ddz_fwd(squeeze(G(4,:,:,:)),dz);
    dGdz(5,:,:,:) = ddz_fwd(squeeze(G(5,:,:,:)),dz);

    % Compute Ubar, convert to primitives
    Ubar = U - dt*(dEdx + dFdy + dGdz);
    [ rho,u,v,w, T,p,e,Et ] = cons2prim( Ubar, const.R,const.cv );


    % Enforce Boundary Conditions on primitives
    [ rho,u,v,w, T,p,~,Et ] = BC( rho,u,v,w, T,p,e,Et, u_Inf, Adiabatic );
    % Update Ubar
    Ubar = prim2cons( rho, u,v,w, T, const.cv );




    %% CORRECTOR STEP
    % Backward FDs for outer derivatives
    % Forward FDs for inner nested derivatives

    % Update physical parameters, primitives updated from PREDICTOR
    mu = sutherland(T);         % dynamic viscosity from Sutherland's law
    k = const.cp/const.Pr * mu; % thermal conductivity


    % Construct flux vectors
    Ebar = flux_E(Ebar, rho,u,v,w, T,p,Et, mu,k, dx,dy,dz, 'Corrector', Adiabatic);
    Fbar = flux_F(Fbar, rho,u,v,w, T,p,Et, mu,k, dx,dy,dz, 'Corrector', Adiabatic);
    Gbar = flux_G(Gbar, rho,u,v,w, T,p,Et, mu,k, dx,dy,dz, 'Corrector', Adiabatic);


    % Precompute flux vector derivatives
    dEbardx(1,:,:,:) = ddx_bwd(squeeze(Ebar(1,:,:,:)),dx);
    dEbardx(2,:,:,:) = ddx_bwd(squeeze(Ebar(2,:,:,:)),dx);
    dEbardx(3,:,:,:) = ddx_bwd(squeeze(Ebar(3,:,:,:)),dx);
    dEbardx(4,:,:,:) = ddx_bwd(squeeze(Ebar(4,:,:,:)),dx);
    dEbardx(5,:,:,:) = ddx_bwd(squeeze(Ebar(5,:,:,:)),dx);

    dFbardy(1,:,:,:) = ddy_bwd(squeeze(Fbar(1,:,:,:)),dy);
    dFbardy(2,:,:,:) = ddy_bwd(squeeze(Fbar(2,:,:,:)),dy);
    dFbardy(3,:,:,:) = ddy_bwd(squeeze(Fbar(3,:,:,:)),dy);
    dFbardy(4,:,:,:) = ddy_bwd(squeeze(Fbar(4,:,:,:)),dy);
    dFbardy(5,:,:,:) = ddy_bwd(squeeze(Fbar(5,:,:,:)),dy);

    dGbardz(1,:,:,:) = ddz_bwd(squeeze(Gbar(1,:,:,:)),dz);
    dGbardz(2,:,:,:) = ddz_bwd(squeeze(Gbar(2,:,:,:)),dz);
    dGbardz(3,:,:,:) = ddz_bwd(squeeze(Gbar(3,:,:,:)),dz);
    dGbardz(4,:,:,:) = ddz_bwd(squeeze(Gbar(4,:,:,:)),dz);
    dGbardz(5,:,:,:) = ddz_bwd(squeeze(Gbar(5,:,:,:)),dz);

    % Compute U at time (n+1)
    U = 0.5*( (U + Ubar) - dt*(dEbardx + dFbardy + dGbardz) );
    [ rho,u,v,w, T,p,e,Et ] = cons2prim( U, const.R,const.cv );


    % Enforce Boundary Conditions on primitives
    [ rho,u,v,w, T,p,e,Et ] = BC( rho,u,v,w, T,p,e,Et, u_Inf, Adiabatic );
    % Update U
    U = prim2cons( rho, u,v,w,  T, const.cv );


    time = toc(time);
    timekeep_calc(iteration+1) = time;




    %% Animate (n+1)
    time = tic;

    % Output iteration count
    iteration = iteration + 1;
    if mod(iteration, print_rate) == 0
        clc; fprintf('Current iteration: %d/%d\n\n', iteration,max_iter)
    end


    % Compute numerical schlieren image
    if (useSchlieren)
        S = schlieren(rho, dx,dy,dz);
    else
        S = [];
    end

    % Update primitive variable fields
    if mod(iteration, update_rate) == 0
        plot_fields(rho,u,v, T,p,e, xx,yy, useSchlieren,S);
    end


    % Store value of convergence variable
    converge(iteration) = eval( converge_str );

    % Update convergence plot (rough)
    if mod(iteration, update_conv) == 0
        plot_convergence(        t(1 : update_conv : iteration), ...
                          converge(1 : update_conv : iteration), ...
                          title_str, converge_str);
    end


    % Write video
    if mod(iteration, video_rate) == 0
        frame = getframe(fig);
        writeVideo(video, frame);
    end


    time = toc(time);
    timekeep_plot(iteration) = time;
    
end


% Update convergence plot with smooth plot
plot_convergence(t,converge, title_str,converge_str);


% Display average time and total time elapsed
time_avg_calc = mean(timekeep_calc)
time_avg_plot = mean(timekeep_plot)
time_all = toc(time_all)


close(video);




%% Functions

function [] = plot_fields( rho,u,v, T,p,e, xx,yy, useSchlieren,S )
% Pseudocolor plot of primitive variables in x,y domain

% Obtain zplane slices
sz = size(rho);
zplane = floor(sz(3)/2);
xx = squeeze(xx(:,:,zplane));
yy = squeeze(yy(:,:,zplane));

subplot(3,3,1);
if ~(useSchlieren)
    pcolor(xx,yy,rho(:,:,zplane));
    title('Density');
    cb_str = '$\rho$ $[\frac{kg}{m^3}]$';
    plot_settings(cb_str, jet,'auto')
else
    pcolor(xx,yy,S(:,:,zplane));
    title('Numerical schlieren image');
    cb_str = '$S(x,y)$ $[1]$';
    plot_settings(cb_str, gray,[0 1])
end

subplot(3,3,2);
pcolor(xx,yy,u(:,:,zplane));
title('Velocity x-component');
cb_str = '$u$ $[\frac{m}{s}]$';
plot_settings(cb_str, jet,'auto')

subplot(3,3,3);
pcolor(xx,yy,v(:,:,zplane));
title('Velocity y-component');
cb_str = '$v$ $[\frac{m}{s}]$';
plot_settings(cb_str, jet,'auto')

subplot(3,3,4);
pcolor(xx,yy,e(:,:,zplane));
title('Internal energy');
cb_str = '$e$ $[J]$';
plot_settings(cb_str, jet,'auto')

subplot(3,3,5);
pcolor(xx,yy,p(:,:,zplane));
title('Pressure');
cb_str = '$p$ $[\frac{N}{m^2}]$';
plot_settings(cb_str, jet,'auto')

subplot(3,3,6);
pcolor(xx,yy,T(:,:,zplane));
title('Temperature');
cb_str = '$T$ $[K]$';
plot_settings(cb_str, hot,'auto')


% % Full size pressure plot
% subplot(3,3,1:6);
% pcolor(xx,yy,p(:,:,zplane));
% title('Pressure');
% cb_str = '$p$ $[\frac{N}{m^2}]$';
% plot_settings(cb_str, jet,'auto')


drawnow

end


function [] = plot_settings( cb_str, cmap,limits )
% Common plot settings

ax = gca;

xlabel('x'); ylabel('y');
cb = colorbar;
cb.Label.Interpreter = 'Latex';
ylabel(cb, cb_str, 'FontSize',14);
colormap(cb,cmap);

colormap(ax,cmap);
clim(ax,limits);

axis equal tight; shading interp

end




function [] = plot_convergence( t,y, title_str,y_str )
% Plot of convergence variable over time

subplot(3,3,7:9);
plot(t,y);
title(title_str); xlabel('t [s]'); ylabel(y_str);
axis tight

% drawnow

end




%% EOF
