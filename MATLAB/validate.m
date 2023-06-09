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

streamwise = floor(Nx.* (0.125:0.2:1));


%% Setup

x = squeeze(xx(:,1,1));
y = squeeze(yy(1,:,1));
z = squeeze(zz(1,1,:));

u0 = u(1,1,1);




%% Plot

% Axial velocity contours, transverse velocity quiver
for i = 1:length(streamwise)

    % Compute
    x_D = x(streamwise);

    axial_vel = squeeze( u(streamwise(i),:,:) );
    y_vel = squeeze( v(streamwise(i),:,:) );
    z_vel = squeeze( w(streamwise(i),:,:) );


    % Plot
    fig_contour(i) = figure;

    contourf(y,z, axial_vel); hold on
    colorbar

    quiver(y,z, y_vel,z_vel)


    title("Velocity contours at x/D=" + x_D(i))

end

%% Plot

% Axial velocity at wall bisector
for j = length(streamwise)

    axial_vel = squeeze( u(streamwise(j),:,floor(Nz/2)) ) ./ u0;


    % Plot
    fig_bisector(j) = figure;

    plot(axial_vel,y)

end


%%




%% EOF
