%**************************************************************************
% Created    : 30.03.2020
% Author     : Andre Breuer
%**************************************************************************
echo off; clear; clc; close all;

%------------------%
% Input Parameters %
%------------------%
nX = 128;
nY = 128;
gamma = 1.4;
% Error in implementation of van Leer Flux... 
fluxType = 'HLL';
testCase = 1;
T = 1;

sim = runSimulation(nX, nY, gamma, T, fluxType, testCase);

[rho, vX, vY, p, E] = sim.getResults();

function sim = runSimulation(nX, nY, gamma, T, fluxType, testCase)
    switch testCase
        case 1
            dx = 2*pi/nX;
            x = -pi:dx:pi;
            x = x(1:end-1);

            dy = 2*pi/nY;
            y = -pi:dy:pi;
            y = y(1:end-1);

            [X, Y] = meshgrid(x,y);

            % Initial conditions
            rho0 = 1 + 0.2 * sin(X) .* cos(Y);
            vX0 = ones(size(X));
            vY0 = ones(size(X));
            p0 = ones(size(X));

            sim = solver(gamma, X, Y, rho0, vX0, vY0, p0, T, fluxType);
        case 2
            x = linspace(0, 0.25, nX);
            y = linspace(0, 1, nY);
            [X, Y] = meshgrid(x, y);
            
            c = 1;
            %Initial conditions
            rho0 = 2*ones(size(X)) .* (Y < 0.5) + ones(size(X)) .* (Y >= 0.5);
            vX0 = zeros(size(X));
            vY0 = -0.025 * c* cos(8*pi*X);
            p0 = (2*Y+1) .* (Y < 0.5) + (Y + 1.5) .* (Y >= 0.5);

            sim = solver(gamma, X, Y, rho0, vX0, vY0, p0, T, fluxType, 'reflectiveFull');
            %sim = solver(gamma, X, Y, rho0, vX0, vY0, p0, T, fluxType, 'RT');
    end
    sim.plotInitialConditions();
    sim.solve();
    sim.animate('vX');
end 

% return
% 
% % Exact solution in density
% nT = 20;
% Z = 1 + 0.5 * sin(2*pi*X) .* cos(2*pi*Y);
% [~, h] = contourf(X, Y, Z);
% %colormap(autumn);
% colorbar;
% for t = linspace(0, T, nT)
% Z = 1 + 0.5 * sin(2*pi*(X-t)) .* cos(2*pi*(Y-t));
% set(h, 'ZData', Z);
% drawnow;
% pause(0.1);
% end
% 
