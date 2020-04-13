%**************************************************************************
% Created    : 30.03.2020
% Author     : Andre Breuer
%**************************************************************************

classdef solver < handle
    properties (Access = private)
        dx;         % step width in x
        dy;         % step width in y
        dt;         % step width in time
        X;          % x vector
        Y;          % x vector
        t;          % time vector
        rho;        % density matrix (space & time)
        vX;         % velocity matrix in x (space & time)
        vY;         % velocity matrix in y (space & time)
        p;          % presuure matrix (space & time)
        E;          % energy matrix (space & time)
        U;          % vector of variables in Euler eqns.
        dUdt;       % time dervative of U
        c;          % speed of sound vector
        MX;         % Mach number matrix (space & time)
        MY;         % Mach number matrix (space & time)
        gamma;      % isentropic coefficient
        nX;         % number of steps in x
        nY;         % number of steps in y
        nT;         % number of time steps
        T;          % end time for simulation
        figP;       % figure handle for pressure
        figRho;     % figure handle for density
        figVX;      % figure handle for velocity
        figVY;      % figure handle for velocity
        figE;       % figure handle for energy
        BC;         % boundary conditions
        r;         % reconstruction handle
        fluxHandle; % flux handle
        fluxType;   % string which chooses the flux function to use
    end

    properties (Constant)
        CFL = 0.8;  % max- CFL-number allowed 
                    % (for stability CFL = 1 is required)
        q = 1.4;    % exponent for controlling resolution and robustness
    end
    
    methods (Access = public)
        function obj = solver(gamma, X, Y, rho0, vX0, vY0, p0, T, fluxType, BC)
            % Constructor
            if nargin < 8
                error('Not enough input arguments!');
            elseif nargin < 9
                fluxType = 'vanLeer';
                BC = 'periodic';
            elseif nargin < 10
                BC = 'periodic';
            end
            
            obj.fluxType = fluxType;
            obj.BC = BC;
            if not(strcmp(obj.BC, 'transmissive')    | ...
                   strcmp(obj.BC, 'periodic')        | ...
                   strcmp(obj.BC, 'reflectiveX')     | ...
                   strcmp(obj.BC, 'reflectiveFull')  | ...
                   strcmp(obj.BC, 'RT'))
                error('Invalid boundary condition given!');
            end
            
            [obj.nY, obj.nX] = size(X);
            obj.gamma = gamma;
            obj.X = X;
            obj.Y = Y;
            obj.dx = X(1,2) - X(1,1);
            obj.dy = Y(2,1) - Y(1,1);
            obj.rho = rho0;
            obj.vX = vX0;
            obj.vY = vY0;
            obj.p = p0;
            obj.E = obj.p / (obj.gamma - 1) + 0.5 * obj.rho .* (obj.vX.^2 + obj.vY.^2);
            obj.c = sqrt(obj.gamma * obj.p ./ obj.rho);
            obj.MX = obj.vX ./ obj.c;
            obj.MY = obj.vY ./ obj.c;
            obj.U = cat(3,obj.rho, obj.rho.*obj.vX, obj.rho.*obj.vY, obj.E);
            [aL, aR] = obj.calculateEigenvalue();
            aMax = max(max(max(abs(aL))), max(max(abs(aR))));
            obj.dt = min(obj.dx, obj.dt) * obj.CFL / aMax;
            %obj.t = (0:obj.dt:T);
            %obj.nT = length(obj.t);
            obj.t = 0;
            obj.nT = 0;
            obj.T = T;
            obj.r = reconstructor(obj.dx, obj.dy, obj.gamma, obj.q);
            obj.fluxHandle = numericalFluxesEuler2D(obj.dx, obj.dy, obj.dt, obj.gamma);
        end
        
        function plotInitialConditions(obj)
            obj.assignFigures();
            obj.plot(obj.figRho, obj.rho, '$\rho$');
            obj.plot(obj.figVX, obj.vX, '$v$');
            obj.plot(obj.figVY, obj.vY, '$v$');
            obj.plot(obj.figP, obj.p, '$p$');
            obj.plot(obj.figE, obj.E, '$E$');
        end
        
        function solve(obj)
            disp('Solving the PDE...')
            f = waitbar(0,'Solving the PDE...');
            tic;
            %for i = 1:obj.nT
            while obj.t(end) < obj.T
                disp(obj.t(end))
                obj.performUpdateStep();
                obj.assignResults();
                waitbar(obj.t(end)/obj.T,f,'Solving the PDE...');
            end
            toc;
            close(f);
            disp('Succesfully computed the solution with ' + string(obj.nT) + ' time steps.')
        end
        
        function [rho, vX, vY, p, E, c, MX, MY, tEnd] = getResults(obj)
            rho = obj.rho;
            vX = obj.vX;
            vY = obj.vY;
            p = obj.p;
            E = obj.E;
            c = obj.c;
            MX = obj.MX;
            MY = obj.MY;
            tEnd = obj.t(end);
        end
        
        function animate(obj, property, playTime)
            if nargin < 3
                playTime = 10;
            end
            switch property
                case 'rho'
                    obj.animatePlot(obj.figRho, obj.rho, playTime);
                case 'vX'
                    obj.animatePlot(obj.figVX, obj.vX, playTime);
                case 'vY'
                    obj.animatePlot(obj.figVY, obj.vY, playTime);
                case 'p'
                    obj.animatePlot(obj.figP, obj.p, playTime);
                case 'E'
                    obj.animatePlot(obj.figE, obj.E, playTime);
                case 'all'
                    obj.animateAll(playTime);
                otherwise
                    error('Invalid property for plotting!')
            end
        end
    end
    
    methods (Access = private)        
        function assignFigures(obj)
            obj.figRho = figure(1);
            obj.setPosition(obj.figRho);
            obj.figVX = figure(2);
            obj.setPosition(obj.figVX);
            obj.figVY = figure(3);
            obj.setPosition(obj.figVY);
            obj.figP = figure(4);
            obj.setPosition(obj.figP);
            obj.figE = figure(5);
            obj.setPosition(obj.figE);
        end
        
        function setPosition(obj, f)
            portion = 0.6;
            screenSizePortion = portion * get(0, 'ScreenSize');
            set(f, 'units','pixels','position', screenSizePortion);
            movegui(f, 'center');
        end
        
        function plot(obj, figureHandle, Z, zLabel)
            labelSize = 12;
            set(0, 'CurrentFigure', figureHandle);
            surf(obj.X, obj.Y, Z);
            %contourf(obj.X, obj.Y, Z);
            xlabel('$x$', 'Interpreter', 'latex', 'Fontsize', labelSize);
            ylabel('$y$', 'Interpreter', 'latex', 'Fontsize', labelSize);
            zlabel(zLabel, 'Interpreter', 'latex', 'Fontsize', labelSize);
            grid minor;
            %hold on;
            xlim([min(min(obj.X)), max(max(obj.X))]);
            ylim([min(min(obj.Y)), max(max(obj.Y))]);
            colorbar;
        end
        
        function animatePlot(obj, figureHandle, Z, playTime)
            labelSize = 12;
            figure(figureHandle);
            dec = 10; % Factor for rounding up or down to decimal places
            set(0, 'CurrentFigure', figureHandle);
            plt = surf(obj.X, obj.Y, Z(:,:,1));
            shading interp;
            zlim([floor(dec*min(min(min(Z))))/dec, ceil(dec*max(max(max(Z))))/dec]);
            colorbar;
            caxis([floor(dec*min(min(min(Z))))/dec, ceil(dec*max(max(max(Z))))/dec]);
%             [~, plt] = contour(obj.X, obj.Y, Z(:,:,1));
            xlabel('$x$', 'Interpreter', 'latex', 'Fontsize', labelSize);
            ylabel('$y$', 'Interpreter', 'latex', 'Fontsize', labelSize);
            
            nTime = length(obj.t);
            
            pauseTime = playTime / nTime;
    
            for i = (2:nTime)
                set(plt, 'ZData', Z(:,:,i))
                drawnow;
                pause(pauseTime);
            end
        end
        
        function animateAll(obj, playTime)
            figure(1);
            dec = 10; % Factor for rounding up or down to decimal places
            set(1, 'CurrentAxes', obj.figRho);
            pltRho = plot(obj.x, obj.rho(1,:), 'Color', 'blue');
            ylim([floor(dec*min(min(obj.rho)))/dec, ceil(dec*max(max(obj.rho)))/dec]);
            set(1, 'CurrentAxes', obj.figV);
            pltV = plot(obj.x, obj.v(1,:), 'Color', 'blue');
            ylim([floor(dec*min(min(obj.v)))/dec, ceil(dec*max(max(obj.v)))/dec]);
            set(1, 'CurrentAxes', obj.figP);
            pltP = plot(obj.x, obj.p(1,:), 'Color', 'blue');
            ylim([floor(dec*min(min(obj.p)))/dec, ceil(dec*max(max(obj.p)))/dec]);
            set(1, 'CurrentAxes', obj.figE);
            pltE = plot(obj.x, obj.E(1,:), 'Color', 'blue');
            ylim([floor(dec*min(min(obj.E)))/dec, ceil(dec*max(max(obj.E)))/dec]);
            
            nTime = length(obj.t);
            
            pauseTime = playTime / nTime;
    
            for i = (2:nTime)
                set(pltRho, 'YData', obj.rho(i,:))
                set(pltV, 'YData', obj.v(i,:))
                set(pltP, 'YData', obj.p(i,:))
                set(pltE, 'YData', obj.E(i,:))
                drawnow;
                pause(pauseTime);
            end
            
            str = strcat('\bf Solution at $t = ', num2str(obj.t(end),2), '$');
            annotation('textbox', [0 0.9 1 0.1], ...
                       'String', str, ...
                       'Interpreter', 'latex', ...
                       'EdgeColor', 'none', ...
                       'FontSize', 14, ...
                       'HorizontalAlignment', 'center');
        end
        
        function performUpdateStep(obj)
            obj.setTimeStep();
            obj.dUdt = obj.calculateRightHandSide(obj.U);
            if obj.CFLconditionSatisfied()
                obj.performTimeUpdate();
            else
                warning('CFL condition is not satisfied!')
            end
        end
        
        function L = calculateRightHandSide(obj, U)
            [UL, UR] = obj.r.reconstructValuesLDLR(U, 'x', obj.BC);
            UL = obj.setGhostCells(UL);
            UR = obj.setGhostCells(UR);
            %FL = F_{i-0.5}; FR = F_{i+0.5}
            FL = obj.fluxHandle.calculateNumericalFlux(UL(1:end-2,2:end-1,:), UR(2:end-1,2:end-1,:), 'F', obj.fluxType);
            FR = obj.fluxHandle.calculateNumericalFlux(UL(2:end-1,2:end-1,:), UR(3:end,2:end-1,:), 'F', obj.fluxType);
            [UL, UR] = obj.r.reconstructValuesLDLR(U, 'y', obj.BC);
            UL = obj.setGhostCells(UL);
            UR = obj.setGhostCells(UR);
            %GL = G_{i-0.5}; GR = G_{i+0.5}
            GL = obj.fluxHandle.calculateNumericalFlux(UL(2:end-1,1:end-2,:), UR(2:end-1,2:end-1,:), 'G', obj.fluxType);
            GR = obj.fluxHandle.calculateNumericalFlux(UL(2:end-1,2:end-1,:), UR(2:end-1,3:end,:), 'G', obj.fluxType);
            L = - (FR - FL) / obj.dx - (GR - GL) / obj.dy;
        end
        
        function UGhost = setGhostCells(obj, U)
            switch obj.BC
                case 'periodic'
                    % Ghost cells in x-direction
                    UGhost = cat(1, U(end,:,:), U, U(1,:,:));
                    % Ghost cells in y-direction
                    UGhost = cat(2, UGhost(:,end,:), UGhost, UGhost(:,1,:));
                case 'transmissive'
                    % Ghost cells in x-direction
                    UGhost = cat(1, U(1,:,:), U, U(end,:,:));
                    % Ghost cells in y-direction
                    UGhost = cat(2, UGhost(:,1,:), UGhost, UGhost(:,end,:));
                case 'reflectiveX'
                    UGhost = cat(1, ...
                             cat(3, U(1,:,1), -U(1,:,2), -U(1,:,3), U(1,:,4)), ...
                             U, ...
                             cat(3, U(end,:,1), -U(end,:,2), -U(end,:,3), U(end,:,4)));
                    % Ghost cells in y-direction
                    UGhost = cat(2, UGhost(:,end,:), UGhost, UGhost(:,1,:));
                case 'reflectiveFull'
                    UGhost = cat(1, ...
                             cat(3, U(1,:,1), -U(1,:,2), -U(1,:,3), U(1,:,4)), ...
                             U, ...
                             cat(3, U(end,:,1), -U(end,:,2), -U(end,:,3), U(end,:,4)));
                    % Ghost cells in y-direction
                    UGhost = cat(2, ...
                             cat(3, UGhost(:,1,1), -UGhost(:,1,2), -UGhost(:,1,3), UGhost(:,1,4)), ...
                             UGhost, ...
                             cat(3, UGhost(:,end,1), -UGhost(:,end,2), -UGhost(:,end,3), UGhost(:,end,4)));
                case 'RT'
                    UGhost = cat(1, ...
                             cat(3, U(1,:,1), -U(1,:,2), -U(1,:,3), U(1,:,4)), ...
                             U, ...
                             cat(3, U(end,:,1), -U(end,:,2), -U(end,:,3), U(end,:,4)));
                    UGhost = cat(2, ...
                             cat(3, 2*ones(size(UGhost(:,1,1))), ...
                                    zeros(size(UGhost(:,1,2))), ...
                                    zeros(size(UGhost(:,1,3))), ...
                                    ones(size(UGhost(:,1,4))) / (obj.gamma - 1)), ...
                             UGhost, ...
                             cat(3, ones(size(UGhost(:,end,1))), ...
                                    zeros(size(UGhost(:,end,2))), ...
                                    zeros(size(UGhost(:,end,3))), ...
                                    2.5*ones(size(UGhost(:,end,4))) / (obj.gamma - 1)));                    
            end
        end

        function setTimeStep(obj)
             [aL, aR] = obj.calculateEigenvalue();
             aMax = max(max(max(abs(aL))), max(max(abs(aR))));
             obj.dt = min(obj.dx, obj.dy) * obj.CFL / aMax;
             obj.fluxHandle.setTimeStep(obj.dt);
             if (obj.t(end) + obj.dt) < obj.T
                obj.t = [obj.t, obj.t(end)+obj.dt];
             else
                obj.dt = obj.T - obj.t(end);
                obj.t = [obj.t, obj.T];
             end
             obj.nT = obj.nT + 1;
        end
       
        function b = CFLconditionSatisfied(obj)
             [aL, aR] = obj.calculateEigenvalue();
             aMax = max(max(abs(aL)), max(abs(aR)));
             if (aMax * obj.dt / obj.dx) <= 1
                b = true;
                return;
             end
             b = false;
        end
        
        function [aL, aR] = calculateEigenvalue(obj)
            v = obj.U(:,:,2) ./ obj.U(:,:,1);
            c = sqrt(obj.gamma * ((obj.gamma - 1) * (obj.U(:,:,4) - 0.5 * ...
                (obj.U(:,:,2).^2 + obj.U(:,:,3).^2) ./ ... 
                obj.U(:,:,1))) ./ obj.U(:,:,1));
            aL = v - c;
            aR = v + c;
        end
        
        function performTimeUpdate(obj)
            % 3rd order SSP-RK method
            U1 = obj.U + obj.dt * obj.dUdt;
            U2 = (3 * obj.U + U1 + ...
                 obj.dt * obj.calculateRightHandSide(U1)) / 4;
            obj.U = (obj.U  + 2 * U2 + ...
                     2 * obj.dt * obj.calculateRightHandSide(U2)) / 3;
        end
        
        function assignResults(obj)
            obj.rho = cat(3, obj.rho, obj.U(:,:,1));
            obj.vX = cat(3, obj.vX,  obj.U(:,:,2) ./ obj.U(:,:,1)); 
            obj.vY = cat(3, obj.vY,  obj.U(:,:,3) ./ obj.U(:,:,1));
            obj.p = cat(3, obj.p, (obj.gamma - 1) * (obj.U(:,:,4) - 0.5 * ...
                    (obj.U(:,:,2).^2 + obj.U(:,:,3).^2) ./ obj.U(:,:,1)));
            obj.E = cat(3, obj.E, obj.U(:,:,4));
            %obj.c = sqrt(obj.gamma * obj.p ./ obj.rho);
            %obj.MX = cat(3, obj.MX, obj.U(:,:,2) ./ (obj.U(:,:,1) .* obj.c));
            %obj.MY = cat(3, obj.MY, obj.U(:,:,3) ./ (obj.U(:,:,1) .* obj.c));
        end
    end
end