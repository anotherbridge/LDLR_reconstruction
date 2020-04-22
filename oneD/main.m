%**************************************************************************
% Created    : 30.03.2020
% Author     : Andre Breuer
%**************************************************************************
echo off; clear; clc; close all;

%------------------%
% Input Parameters %
%------------------%
n = 1000;
gamma = 1.4;
testCase = 5;
numericalFlux = 'HLL';
reconstructionMethod = 'LDLR';
integrationMethod = 'RK3';

sim = runSimulation(n, gamma, testCase, numericalFlux, reconstructionMethod, integrationMethod);

function sim = runSimulation(n, gamma, testCase, numericalFlux, reconstructionMethod, integrationMethod, animationTime)
    if nargin == 7
        animTime = animationTime; 
    else
        animTime = 10; % in seconds
    end
    
    switch testCase
        case 1 % Sod problem
            T = 0.2;
            x = linspace(0, 1, n);
            rho0 = ones(1, n);
            v0 = zeros(1, n);
            p0 = ones(1, n);
            rho0(x >= 0.5) = 0.125;
            p0(x >= 0.5) = .1;
            sim = solver(gamma, x, rho0, v0, p0, T, numericalFlux, ...
                         'transmissive', reconstructionMethod, ...
                         integrationMethod);
            sim.plotInitialConditions();
        case 2 
            T = 1;
            % Correction so that the endpoint is only included once for PBC
            dx = 2*pi/n;
            x = -pi:dx:pi;
            x = x(1:end-1);
            rho0 = 1 + 0.2*sin(x);
            v0 = ones(size(x));
            p0 = ones(size(x));
            sim = solver(gamma, x, rho0, v0, p0, T, numericalFlux, ...
                         'periodic', reconstructionMethod, ...
                         integrationMethod);
            sim.plotInitialConditions();
            subplot(2,2,1);
        case 3 % Lax problem
            T = 1;
            x = linspace(0, 1, n);
            rho0 = 0.445 * ones(size(x));
            v0 = 0.698 * ones(size(x));
            p0 = 3.528 * ones(size(x));
            rho0(x >= 0.5) = 0.5;
            v0(x >= 0.5) = 0;
            p0(x >= 0.5) = 0.571;
            sim = solver(gamma, x, rho0, v0, p0, T, numericalFlux, ...
                         'reflectiveRight', reconstructionMethod, ...
                         integrationMethod);
            sim.plotInitialConditions();
        case 4 % Lax problem
            T = 1;
            x = linspace(0, 1, n);
            rho0 = 0.445 * ones(size(x));
            v0 = 0.698 * ones(size(x));
            p0 = 3.528 * ones(size(x));
            rho0(x >= 0.5) = 0.5;
            v0(x >= 0.5) = 0;
            p0(x >= 0.5) = 0.571;
            sim = solver(gamma, x, rho0, v0, p0, T, numericalFlux, ...
                         'reflectiveFull', reconstructionMethod, ...
                         integrationMethod);
            sim.plotInitialConditions();
        case 5 % Shu-Osher shock-acoustic problem
            T = 1.8;
            x = linspace(-5, 5, n);
            rho0 = 3.857143 * ones(size(x));
            v0 = 2.629369 * ones(size(x));
            p0 = 10.333333 * ones(size(x));
            rho0(x >= -4) = 1 + 0.2*sin(5*x(x >= -4));
            v0(x >= -4) = 0;
            p0(x >= -4) = 1;
            sim = solver(gamma, x, rho0, v0, p0, T, numericalFlux, ...
                         'transmissive', reconstructionMethod, ...
                         integrationMethod);
            sim.plotInitialConditions();
        otherwise
            error('Given test case is not valid!')
    end
    
    sim.solve();
    % sim.animate('rho');
    % sim.animate('v');
    % sim.animate('p');
    % sim.animate('E');
    sim.animate('all', animTime);
end
