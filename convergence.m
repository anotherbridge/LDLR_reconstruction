%**************************************************************************
% Created    : 31.03.2020
% Author     : Andre Breuer
%**************************************************************************
echo off; clear; clc; close all;

% Convergence is going to be checked at time T=1

lblSize = 12;
nVec = [8, 16, 32, 64, 128, 256, 512, 1024];
gamma = 1.4;
T = 1;
numericalFlux = 'vanLeer';

dxVec  = [];
err    = [];

disp('Error calculation...')
for n = nVec
    fprintf('\nn = %d\n', n);
    % Correction so that the endpoint is only included once for PBC
    dx = 2*pi/n;
    x = -pi:dx:pi;
    x = x(1:end-1);
    dxVec = [dxVec, dx];
        
    % Compute numerical solution
    rho0 = 1 + 0.2*sin(x);
    v0 = ones(size(x));
    p0 = ones(size(x));
    sim = solver(gamma, x, rho0, v0, p0, T, numericalFlux, 'periodic'); 
    sim.solve();
    [rho, v, p, ~, ~, ~, tEnd] = sim.getResults();
    
    % Compute analytic solution in density (v = v0, p = p0)
    rhoEx = 1 + 0.2 * sin(x - tEnd);
    
    uErr = [abs(rho(end,:) - rhoEx);
             abs(v(end,:) - v0);
             abs(p(end,:) - p0)];
    % Calculate L1 error
    err = [err, sum(uErr, 2)*dx];
end

% Plot the error
for i = 1:3
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
            str = '$|| p_{\rm exact} - p ||_{L^1 ([-\pi,\pi]xT)}$';
        case 3
            str = '$|| v_{\rm exact} - v ||_{L^1 ([-\pi,\pi]xT)}$';
    end
    ylabel(str, 'Interpreter', 'latex', 'FontSize', lblSize);
    plot(x_fitlog, y_fitlog, 'b-');
    xlim([min(dxVec),max(dxVec)]);
    title('Order of convergence: $\Delta x^{' + string(abs(p(1))) + '}$', 'Interpreter','latex', 'FontSize', lblSize+4);
    grid on
end

figure(4)
plot(x, rho(end,:), 'r-')
hold on
grid on
xlim([min(x), max(x)])
xlabel('$x$', 'Interpreter', 'latex', 'Fontsize', lblSize);
ylabel('$\rho$', 'Interpreter', 'latex', 'Fontsize', lblSize);
plot(x, 1 + 0.2 * sin(x - tEnd), 'b-')
legend({'$\rho$', '$\rho_{\rm{exact}}$'},'Interpreter','latex', 'FontSize', 14, 'Location', 'best')