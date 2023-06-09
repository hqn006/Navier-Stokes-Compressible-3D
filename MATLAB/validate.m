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

search = [0.5 1.5 2.5 3.5];


%% Setup

x = squeeze(xx(:,1,1));
y = squeeze(yy(1,:,1));
z = squeeze(zz(1,1,:));

u_cl = u(Nx/2,Ny/2,Nz/2);


H = 1 * units.in2mm*10^-3; % Height, [m]
W = 1 * units.in2mm*10^-3; % Width,  [m]
L = 6 * H; % Length, [m]
dx = L / (Nx-1);
delta = dx ./ units.in2mm.*10^3;




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


    % Compute
    x_D(i) = x(streamwise(i));

    axial_vel = squeeze( u(streamwise(i),:,:) ) ./ u_cl;
    y_vel = squeeze( v(streamwise(i),:,:) );
    z_vel = squeeze( w(streamwise(i),:,:) );


    % Plot
    fig_contour(i) = figure;

    contourf(z,y, axial_vel); hold on
    colorbar

    quiver(z,y, z_vel,y_vel)


    title("Velocity contours at x/D=" + x_D(i))

end

%% Plot

% Axial velocity at wall bisector
for j = length(streamwise)

    axial_vel = squeeze( u(streamwise(j),:,floor(Nz/2)) ) ./ u_cl;


    % Plot
    fig_bisector(j) = figure;

    plot(axial_vel,y)

end


%%




%% EOF
