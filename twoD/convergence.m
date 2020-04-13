%**************************************************************************
% Created    : 31.03.2020
% Author     : Andre Breuer
%**************************************************************************
echo off; clear; clc; close all;

% Convergence is going to be checked at time T=1

lblSize = 12;
nVec = [8, 16, 32, 64, 128];%, 256];
gamma = 1.4;
T = 1;
numericalFlux = 'HLL';

dxVec  = [];
err    = [];

disp('Error calculation...')
for n = nVec
    fprintf('\nn = %d\n', n);
    % Correction so that the endpoint is only included once for PBC    
    dx = 2*pi/n;
    x = -pi:dx:pi;
    x = x(1:end-1);
    y = -pi:dx:pi;
    y = y(1:end-1);
    dxVec = [dxVec, dx];
    
    [X, Y] = meshgrid(x,y);

    % Compute numerical solution
    rho0 = 1 + 0.5 * sin(X) .* cos(Y);
    vX0 = ones(size(X));
    vY0 = ones(size(X));
    p0 = ones(size(X));

    sim = solver(gamma, X, Y, rho0, vX0, vY0, p0, T, numericalFlux);
    sim.solve();
    [rho, vX, vY, p, ~, ~, ~, ~, tEnd] = sim.getResults();
    
    % Compute analytic solution in density (v = v0, p = p0)
    rhoEx = 1 + 0.5 * sin((X-tEnd)) .* cos((Y-tEnd));
    
    uErr = [sum(sum(abs(rho(:,:,end) - rhoEx)))*dx^2;
             sum(sum(abs(vX(:,:,end) - vX0)))*dx^2;
             sum(sum(abs(vY(:,:,end) - vY0)))*dx^2;
             sum(sum(abs(p(:,:,end) - p0)))*dx^2];
    % Calculate L1 error
    err = [err, uErr];
end

% Plot the error
for i = 1:4
    logx = log(dxVec);
    logy = log(err(i,:));
    p = polyfit(logx, logy, 1);
    x_fitlin = linspace(log(min(dxVec)),log(max(dxVec)),2);
    y_fitlin = p(1)*x_fitlin + p(2);
    x_fitlog = exp(x_fitlin);
    y_fitlog = exp(y_fitlin);

    figure(i);
    loglog(dxVec, err(i,:), 'rx');
    hold on
    xlabel('$\Delta x$', 'Interpreter', 'latex', 'FontSize', lblSize);
    switch i
        case 1 
            str = '$|| \rho - \rho_{\rm exact} ||_{L^1 ([-\pi,\pi]xT)}$';
        case 2
            str = '$|| v_{x \rm exact} - v ||_{L^1 ([-\pi,\pi]xT)}$';
        case 3
            str = '$|| v_{y \rm exact} - v ||_{L^1 ([-\pi,\pi]xT)}$';
        case 4
            str = '$|| p_{\rm exact} - p ||_{L^1 ([-\pi,\pi]xT)}$';
    end
    ylabel(str, 'Interpreter', 'latex', 'FontSize', lblSize);
    plot(x_fitlog, y_fitlog, 'b-');
    xlim([min(dxVec),max(dxVec)]);
    title('Order of convergence: $\Delta x^{' + string(abs(p(1))) + '}$', 'Interpreter','latex', 'FontSize', lblSize+4);
    grid on
end

% figure(4)
% plot(x, rho(end,:), 'r-')
% hold on
% grid on
% xlim([min(x), max(x)])
% xlabel('$x$', 'Interpreter', 'latex', 'Fontsize', lblSize);
% ylabel('$\rho$', 'Interpreter', 'latex', 'Fontsize', lblSize);
% plot(x, 1 + 0.2 * sin(x - tEnd), 'b-')
% legend({'$\rho$', '$\rho_{\rm{exact}}$'},'Interpreter','latex', 'FontSize', 14, 'Location', 'best')