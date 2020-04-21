%**************************************************************************
% Created    : 30.03.2020
% Author     : Andre Breuer
%**************************************************************************
echo off; clear; clc; close all;

n = 10;
dx = 2*pi/n;
x = -pi:dx:pi;
x = x(1:end-1);
u = 1 + 0.2*sin(x);

% x = linspace(0,1, n);
% u = ones(1, n);
% u(x >= 0.5) = 0.125;

r = reconstructor(dx);
[uL, uR] = r.reconstructValuesLinear(u, 'periodic', 'vanLeer');
%[uL, uR] = r.reconstructValuesLDLR(u, 'periodic');

xx = linspace(-1.1*pi, pi, 100);
uu = 1 + 0.2*sin(xx);

% m = 200;
% xx = linspace(0,1,m);
% uu = ones(1, m);
% uu(xx >= 0.5) = 0.125;

figure(1)
xP = x - 0.5 * dx;
stairs(xP, u, 'b-');
grid on
hold on
plot(xP, circshift(uL,1,2), 'ro')
plot(xP, uR, 'go')
plot(xx, uu, 'b--')
plot(xP, circshift(uL,1,2), 'r--')
plot(xP, uR, 'g--')
xlim([min(xP)-0.05, max(xx)])
legend({'$u_{\Delta x}$', '$u^{\rm L}$', '$u^{\rm R}$', '$u_{\rm{exact}}$'},'Interpreter','latex', 'FontSize', 14, 'Location', 'best')
xlabel('$x$', 'Interpreter', 'latex');
ylabel('$u$', 'Interpreter', 'latex');
title('Visualization of LDLR');