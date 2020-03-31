%**************************************************************************
% Created    : 30.03.2020
% Author     : Andre Breuer
%**************************************************************************
echo off; clear; clc; close all;

%------------------%
% Input Parameters %
%------------------%
n = 500;
gamma = 1.4;
testCase = 1;
numericalFlux = 'vanLeer';

runSimulation(n, gamma, testCase, numericalFlux);

function runSimulation(n, gamma, testCase, numericalFlux, animationTime)
    if nargin == 6
        animTime = animationTime; 
    else
        animTime = 10; % in seconds
    end
    
    switch testCase
        case 1
            T = 0.2;
            x = linspace(0, 1, n);
            rho0 = ones(1, n);
            v0 = zeros(1, n);
            p0 = ones(1, n);
            rho0(x >= 0.5) = 0.125;
            p0(x >= 0.5) = .1;
            sim = solver(gamma, x, rho0, v0, p0, T, numericalFlux);
            sim.plotInitialConditions();
        case 2
            T = 1;
            x = linspace(-pi, pi, n);
            rho0 = 1 + 0.2*sin(x);
            v0 = ones(size(x));
            p0 = ones(size(x));
            sim = solver(gamma, x, rho0, v0, p0, T, numericalFlux, 'periodic');
            sim.plotInitialConditions();
            subplot(2,2,1);
            ylim([0.8,1.2]);
            subplot(2,2,2);
            ylim([0,2]);
            subplot(2,2,3);
            ylim([0,2]);
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
