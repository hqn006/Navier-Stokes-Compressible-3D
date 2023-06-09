% MAE 290C
% Project : MacCormack for compressible Navier-Stokes
%
% - Conservative vector (flux) form of 3D compressible Navier-Stokes
% - Assume calorically perfect gas
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
set(0, 'DefaultAxesFontSize', 14);




%% Input

M = 3.9; % Inflow Mach number = u_Inf/a, [1]


% Number of grid points
Nx = 220;
Ny = 200;
Nz = 200;

% Domain (x,y,z : L x H x W)
% 1 in. x 1 in. (1 in. = 25.4 mm)
H = 1 * units.in2mm*10^-3; % Height, [m]
W = 1 * units.in2mm*10^-3; % Width,  [m]
L = 6 * H; % Length, [m]


% Time span
dt = 4 * 10^-7;
final_time = 1.2 * 10^-4;

% Number of iterations
max_iter = floor(final_time / dt) + 1


% Convergence variable
converge_name = 'u';

% Update every _ iterations
update_rate = 30; % variable field plots
update_conv = 1; % convergence plot
print_rate = 1;
video_rate = update_rate;


% Figure settings
fig_size = [1920 1080];


%% Enable/Disable

% Adiabatic wall
Adiabatic = true; % Zero heat flux through wall BC

% Numerical schlieren ('streaks')
useSchlieren = true; % Additional plot & video of numerical schlieren




%% Setup

% Grid spacing
dx = L / (Nx-1);
dy = H / (Ny-1);
dz = W / (Nz-1);
x = 0 : dx : L;
y = 0 : dy : H;
z = 0 : dz : W;

% Plotting grid
[xx,yy,zz] = ndgrid(x,y,z);
xx = xx ./ units.in2mm.*10^3;
yy = yy ./ units.in2mm.*10^3;
zz = zz ./ units.in2mm.*10^3;


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
Ubar = zeros(size(U));
E = zeros(size(U));
F = zeros(size(U));
G = zeros(size(U));

% Preallocate derivatives of flux vectors
dEdx = zeros(size(U));
dFdy = zeros(size(U));
dGdz = zeros(size(U));


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
video1 = VideoWriter('Navier-Stokes_3D_QoI.mp4', 'MPEG-4');
video2 = VideoWriter('Navier-Stokes_3D_schlieren.mp4', 'MPEG-4');

open(video1); fprintf("Opened Navier-Stokes_3D_QoI.mp4\n")
open(video2); fprintf("Opened Navier-Stokes_3D_schlieren.mp4")




%% Solve
% MacCormack Method for compressible Navier-Stokes

% Plot Initial Condition
fig(1) = figure('Position', [0 0 fig_size]);
plot_fields(rho,u,v, T,p,e, xx,yy)


% Compute and plot numerical schlieren
if (useSchlieren)
    S = schlieren(rho, dx,dy,dz);

    fig(2) = figure('Position', [0 100 fig_size]);
    plot_schlieren(S, xx,yy)
end


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


    % Construct flux vectors: E,F,G
    E = flux_E(E, rho,u,v,w, T,p,Et, mu,k, dx,dy,dz, 'Predictor', Adiabatic);
    F = flux_F(F, rho,u,v,w, T,p,Et, mu,k, dx,dy,dz, 'Predictor', Adiabatic);
    G = flux_G(G, rho,u,v,w, T,p,Et, mu,k, dx,dy,dz, 'Predictor', Adiabatic);


    % Precompute flux vector derivatives: dE/dx, dF/dy, dG/dz
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


    % Construct flux vectors: Ebar,Fbar,Gbar
    E = flux_E(E, rho,u,v,w, T,p,Et, mu,k, dx,dy,dz, 'Corrector', Adiabatic);
    F = flux_F(F, rho,u,v,w, T,p,Et, mu,k, dx,dy,dz, 'Corrector', Adiabatic);
    G = flux_G(G, rho,u,v,w, T,p,Et, mu,k, dx,dy,dz, 'Corrector', Adiabatic);


    % Precompute flux vector derivatives: dEbar/dx, dFbar/dy, dGbar/dz
    dEdx(1,:,:,:) = ddx_bwd(squeeze(E(1,:,:,:)),dx);
    dEdx(2,:,:,:) = ddx_bwd(squeeze(E(2,:,:,:)),dx);
    dEdx(3,:,:,:) = ddx_bwd(squeeze(E(3,:,:,:)),dx);
    dEdx(4,:,:,:) = ddx_bwd(squeeze(E(4,:,:,:)),dx);
    dEdx(5,:,:,:) = ddx_bwd(squeeze(E(5,:,:,:)),dx);

    dFdy(1,:,:,:) = ddy_bwd(squeeze(F(1,:,:,:)),dy);
    dFdy(2,:,:,:) = ddy_bwd(squeeze(F(2,:,:,:)),dy);
    dFdy(3,:,:,:) = ddy_bwd(squeeze(F(3,:,:,:)),dy);
    dFdy(4,:,:,:) = ddy_bwd(squeeze(F(4,:,:,:)),dy);
    dFdy(5,:,:,:) = ddy_bwd(squeeze(F(5,:,:,:)),dy);

    dGdz(1,:,:,:) = ddz_bwd(squeeze(G(1,:,:,:)),dz);
    dGdz(2,:,:,:) = ddz_bwd(squeeze(G(2,:,:,:)),dz);
    dGdz(3,:,:,:) = ddz_bwd(squeeze(G(3,:,:,:)),dz);
    dGdz(4,:,:,:) = ddz_bwd(squeeze(G(4,:,:,:)),dz);
    dGdz(5,:,:,:) = ddz_bwd(squeeze(G(5,:,:,:)),dz);

    % Compute U at time (n+1)
    U = 0.5*( (U + Ubar) - dt*(dEdx + dFdy + dGdz) );
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
        clc
        fprintf('Current iteration: %d/%d\n', iteration,max_iter)
        fprintf('Current time: %.4e [s]\n\n', t(iteration))
    end


    % Store value of convergence variable
    converge(iteration) = eval( converge_str );

    % Update all field plots
    if mod(iteration, update_rate) == 0

        % Plot primitive variables
        figure(fig(1));
        plot_fields(rho,u,v, T,p,e, xx,yy);
        plot_convergence(        t(1 : update_conv : iteration), ...
                          converge(1 : update_conv : iteration), ...
                          title_str, converge_str);
        drawnow


        % Compute and plot numerical schlieren
        if (useSchlieren)
            S = schlieren(rho, dx,dy,dz);
        
            figure(fig(2));
            plot_schlieren(S, xx,yy)
            plot_convergence(        t(1 : update_conv : iteration), ...
                              converge(1 : update_conv : iteration), ...
                              title_str, converge_str);
            drawnow
        end

    end


    % Write video
    if mod(iteration, video_rate) == 0
        frame = getframe(fig(1));
        writeVideo(video1, frame);

        frame = getframe(fig(2));
        writeVideo(video2, frame);
    end


    time = toc(time);
    timekeep_plot(iteration) = time;
    
end


%% Post-Run

% Update convergence plot with smooth plot
plot_convergence(t,converge, title_str,converge_str);


% Display average time and total time elapsed
time_avg_calc = seconds( mean(timekeep_calc) );
time_avg_plot = seconds( mean(timekeep_plot) );
time_all = seconds( toc(time_all) );
time_all.Format = 'mm:ss';

fprintf("Average calculation time: " + char(time_avg_calc) + "\n")
fprintf("Average plotting time:    " + char(time_avg_plot) + "\n")
fprintf("\nTotal time: " + char(time_all) + " min\n\n\n")


% Save videos
close(video1); fprintf("Saved 'Navier-Stokes_3D_QoI.mp4'.\n")
close(video2); fprintf("Saved 'Navier-Stokes_3D_schlieren.mp4'.\n\n")


% Save variables
save prim.mat  rho u v w T p e Et  xx yy zz  Nx Ny Nz  dx dy dz
fprintf("Saved variables to 'prim.mat'.\n\n")




%% Functions

function [] = plot_fields( rho,u,v, T,p,e, xx,yy)
% Pseudocolor plot of primitive variables in x,y domain; slice of z domain

% Obtain zplane slices
sz = size(rho);
zplane = floor(sz(3)/2);
xx = squeeze(xx(:,:,zplane));
yy = squeeze(yy(:,:,zplane));


% Density
subplot(3,3,1);
pcolor(xx,yy,rho(:,:,zplane));
title('Density');
cb_str = '$\rho$ $[\frac{kg}{m^3}]$';
plot_settings(cb_str, jet,'auto')

% Velocity x-component (u)
subplot(3,3,2);
pcolor(xx,yy,u(:,:,zplane));
title('Velocity x-component');
cb_str = '$u$ $[\frac{m}{s}]$';
plot_settings(cb_str, jet,'auto')

% Velocity y-component (v)
subplot(3,3,3);
pcolor(xx,yy,v(:,:,zplane));
title('Velocity y-component');
cb_str = '$v$ $[\frac{m}{s}]$';
plot_settings(cb_str, jet,'auto')

% Internal energy
subplot(3,3,4);
pcolor(xx,yy,e(:,:,zplane));
title('Internal energy');
cb_str = '$e$ $[J]$';
plot_settings(cb_str, jet,'auto')

% Pressure
subplot(3,3,5);
pcolor(xx,yy,p(:,:,zplane));
title('Pressure');
cb_str = '$p$ $[\frac{N}{m^2}]$';
plot_settings(cb_str, jet,'auto')

% Temperature
subplot(3,3,6);
pcolor(xx,yy,T(:,:,zplane));
title('Temperature');
cb_str = '$T$ $[K]$';
plot_settings(cb_str, hot,'auto')

end


function [] = plot_schlieren( S, xx,yy )
% Pseudocolor plot of schlieren in x,y domain; slice of z domain

% Obtain zplane slices
sz = size(S);
zplane = floor(sz(3)/2);
xx = squeeze(xx(:,:,zplane));
yy = squeeze(yy(:,:,zplane));


% Plot schlieren
subplot(3,3,1:6);
pcolor(xx,yy,S(:,:,zplane));
title('Numerical schlieren image');
cb_str = '$S(x,y)$ $[1]$';
plot_settings(cb_str, gray,[0 1])

end


function [] = plot_settings( cb_str, cmap,limits )
% Common plot settings

ax = gca;

xlabel('x [in.]'); ylabel('y [in.]');
cb = colorbar;
cb.Label.Interpreter = 'Latex';
ylabel(cb, cb_str, 'FontSize',16);
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

end




%% EOF
