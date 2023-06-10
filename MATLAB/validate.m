% MAE 290C
% Project : MacCormack for compressible Navier-Stokes
%
% Extract variables at specific locations for validation.
% Requires 'prim.mat' to be generated from 'main.m'.
%
%
% Huy Nguyen, Blake Carter
% Created: 8 June 2023


clearvars
close all
clc

addpath('classes')
addpath('functions')


load('prim.mat')


format long
set(0, 'DefaultTextInterpreter','Latex')
set(0, 'DefaultAxesFontSize', 14);




%% Input

search = [0.5 1.5 2.5 3.5 5.3 20];


% Figure settings
contour_size = [600 600];
bisector_size = [500 700];


%% Setup

x = squeeze(xx(:,1,1));
y = squeeze(yy(1,:,1));
z = squeeze(zz(1,1,:));

% For corner bisector
y_diag = diag(squeeze(yy(1,:,:)));
z_diag = diag(squeeze(zz(1,:,:)));

delta = dx ./ units.in2m; % range to find search(i)




%% Plot

% Axial velocity contours, transverse velocity quiver
streamwise = zeros(length(search),1);
current = 0;
index = 1;
for i = 1:length(search)

    while current < search(i)
        current = current + delta;
        index = index + 1;
    end
    streamwise(i) = index;
    if index > Nx
        break
    end


    % Compute
    u_cl = u(streamwise(i),Ny/2,Nz/2); % centerline velocity

    axial_vel = squeeze( u(streamwise(i),:,:) ) ./ u_cl; % u contour
    pressure = squeeze( p(streamwise(i),:,:) ); % pressure contour

    y_vel = squeeze( v(streamwise(i),:,:) );
    z_vel = squeeze( w(streamwise(i),:,:) );


    % Plot Velocity
    fig_contour(i) = figure('Position', [0 0 contour_size]);
    contourf(z,y, axial_vel); hold on
    colorbar
    quiver(z,y, z_vel,y_vel)

    title("Velocity contours at x/D=" + x(streamwise(i)))
    xlabel('z')
    ylabel('y')
    axis square


    % Plot Pressure
    fig_p_contour(i) = figure('Position', ...
                               [contour_size(1) 0 contour_size]);
    contourf(z,y, pressure);
    colorbar

    title("Pressure contours at x/D=" + x(streamwise(i)))
    xlabel('z')
    ylabel('y')
    axis square

end

%% Plot

% Axial velocity at wall bisector and corner bisector
j = length(streamwise);
for i = 1:2

    % Compute
    u_cl = u(streamwise(j),Ny/2,Nz/2); % centerline velocity

    axial_vel_through_y = ...
             squeeze( u(streamwise(j),1:floor(Ny/2),floor(Nz/2)) ) ./ u_cl;
    axial_vel_corner = ...
     diag(squeeze( u(streamwise(j),1:floor(Ny/2),1:floor(Nz/2)) )) ./ u_cl;


    % Plot wall bisector
    fig_bisector(j) = figure('Position', [0 contour_size(2) bisector_size]);
    plot(axial_vel_through_y, y(1:floor(Ny/2)).*2)

    title("Wall bisector velocity at x/D=" + x(streamwise(j)))
    xlabel('$U/U_{cl}$')
    ylabel('$y/a$')


    % Plot corner bisector
    fig_corner(j) = figure('Position', [bisector_size(1) contour_size(2) ...
                                        bisector_size]);
    plot(axial_vel_corner, y_diag(1:floor(Ny/2)).*2)

    title("Corner bisector velocity at x/D=" + x(streamwise(j)))
    xlabel('$U/U_{cl}$')
    ylabel('$y/a$')


    % Next iteration
    j = j - 1;

end


%%




%% EOF
