%**************************************************************************
% Created    : 30.03.2020
% Author     : Andre Breuer
%**************************************************************************

classdef solver < handle
    properties (Access = private)
        dx;         % step width in space
        dt;         % step width in time
        x;          % x vector
        t;          % time vector
        rho;        % density matrix (space & time)
        v;          % velocity matrix (space & time)
        p;          % presuure matrix (space & time)
        E;          % energy matrix (space & time)
        U;          % vector of variables in Euler eqns.
        dUdt;       % time dervative of U
        c;          % speed of sound vector
        M;          % Mach number matrix (space & time) 
        gamma;      % isentropic coefficient
        n;          % number of steps
        nT;         % number of time steps
        T;          % end time for simulation
        testCase;   % number of test case
        figP;       % figure handle for pressure
        figRho;     % figure handle for density
        figV;       % figure handle for velocity
        figE;       % figure handle for energy
        BC;         % boundary conditions
        r;          % reconstruction handle
        fluxHandle; % flux handle
        fluxType;   % string which chooses the flux function to use
    end

    properties (Constant)
        CFL = 0.6;  % max- CFL-number allowed 
                    % (for stability CFL = 1 is required)
        q = 1.4;    % exponent for controlling resolution and robustness
    end
    
    methods (Access = public)
        function obj = solver(gamma, x, rho0, v0, p0, T, fluxType, BC)
            % Constructor
            if nargin < 6
                error('Not enough input arguments!');
            elseif nargin < 7
                fluxType = 'HLL';
            elseif nargin < 8
                BC = 'none';
            end
            
            obj.fluxType = fluxType;
            obj.BC = BC;
            if not(strcmp(obj.BC,'none') | strcmp(obj.BC, 'periodic'))
                error('Invalid boundary condition given!');
            end
            
            obj.n = length(x);
            obj.gamma = gamma;
            obj.x = x;
            obj.dx = obj.x(2) - obj.x(1);
            obj.rho = rho0;
            obj.v = v0;
            obj.p = p0;
            obj.E = obj.p / (obj.gamma - 1) + 0.5 * obj.rho .* obj.v.^2;
            obj.c = sqrt(obj.gamma * obj.p ./ obj.rho);
            obj.M = obj.v ./ obj.c;
            obj.U = [obj.rho; obj.rho .* obj.v; obj.E];
            [aL, aR] = obj.calculateEigenvalue();
            aMax = max(max(abs(aL)), max(abs(aR)));
            obj.dt = obj.dx * obj.CFL / aMax;
            %obj.t = (0:obj.dt:T);
            %obj.nT = length(obj.t);
            obj.t = 0;
            obj.nT = 0;
            obj.T = T;
            obj.r = reconstructor(obj.dx, obj.q);
            obj.fluxHandle = numericalFluxesEuler(obj.dx, obj.dt, obj.gamma);
        end
        
        function plotInitialConditions(obj)
            obj.assignFigures();
            obj.plot(obj.figRho, obj.rho, '$\rho$');
            obj.plot(obj.figV, obj.v, '$v$');
            obj.plot(obj.figP, obj.p, '$p$');
            obj.plot(obj.figE, obj.E, '$E$');
        end
        
        function solve(obj)
            disp('Solving the System...')
            %f = waitbar(0,'Solving the System...');
            tic;
            while obj.t(end) < obj.T
                obj.performUpdateStep();
                obj.assignResults();
                %waitbar(i/obj.nT,f,'Solving the System...');
            end
            toc;
            %close(f);
            disp('Succesfully computed the solution with ' + string(obj.nT) + ' time steps.')
        end
        
        function animate(obj, property, playTime)
            if nargin < 3
                playTime = 10;
            end
            switch property
                case 'rho'
                    obj.animatePlot(obj.figRho, obj.rho, playTime);
                case 'v'
                    obj.animatePlot(obj.figV, obj.v, playTime);
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
            f = figure(1);
            portion = 0.6;
            screenSizePortion = portion * get(0, 'ScreenSize');
            set(f, 'units','pixels','position', screenSizePortion);
            movegui(f, 'center');
            obj.figRho = subplot(2,2,1);
            obj.figV = subplot(2,2,2);
            obj.figP = subplot(2,2,3);
            obj.figE = subplot(2,2,4);
        end
        
        function plot(obj, figureHandle, y, yLabel)
            labelSize = 12;
            set(1, 'CurrentAxes', figureHandle);
            plot(obj.x, y, 'Color', 'red');
            xlabel('$x$', 'Interpreter', 'latex', 'Fontsize', labelSize);
            ylabel(yLabel, 'Interpreter', 'latex', 'Fontsize', labelSize);
            grid minor;
            hold on;
            xlim([min(obj.x), max(obj.x)]);
        end
        
        function animatePlot(obj, figureHandle, y, playTime)
            figure(1);
            set(1, 'CurrentAxes', figureHandle);
            plt = plot(obj.x, y(1,:), 'Color', 'blue'); 
            
            nTime = length(obj.t);
            
            pauseTime = playTime / nTime;
    
            for i = (2:nTime)
                set(plt, 'YData', y(i,:))
                drawnow;
                pause(pauseTime);
            end
        end
        
        function animateAll(obj, playTime)
            figure(1);
            set(1, 'CurrentAxes', obj.figRho);
            pltRho = plot(obj.x, obj.rho(1,:), 'Color', 'blue');
            set(1, 'CurrentAxes', obj.figV);
            pltV = plot(obj.x, obj.v(1,:), 'Color', 'blue');
            set(1, 'CurrentAxes', obj.figP);
            pltP = plot(obj.x, obj.p(1,:), 'Color', 'blue');
            set(1, 'CurrentAxes', obj.figE);
            pltE = plot(obj.x, obj.E(1,:), 'Color', 'blue');
            
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
                error('CFL condition is not satisfied!')
            end
        end
        
        function L = calculateRightHandSide(obj, U)
            if strcmp(obj.BC, 'periodic')
                [UL, UR] = obj.r.reconstructValuesPeriodicLDLR(U);
                FL = obj.fluxHandle.calculateNumericalFlux(circshift(UL,1,2), UR, obj.fluxType);
                FR = obj.fluxHandle.calculateNumericalFlux(UL, circshift(UR,-1,2), obj.fluxType);
            else
                [UL, UR] = obj.r.reconstructValuesLDLR(U);
                [~, URshiftR] = obj.r.reconstructValuesLDLR([U(:,2:end), U(:,end)]);
                [ULshiftL, ~] = obj.r.reconstructValuesLDLR([U(:,1), U(:,1:end-1)]);
                FL = obj.fluxHandle.calculateNumericalFlux(ULshiftL, UR, obj.fluxType);
                FR = obj.fluxHandle.calculateNumericalFlux(UL, URshiftR, obj.fluxType);
            end
            L = - (FR - FL) / obj.dx;
        end
        
        function setTimeStep(obj)
             [aL, aR] = obj.calculateEigenvalue();
             aMax = max(max(abs(aL)), max(abs(aR)));
             obj.dt = obj.dx * obj.CFL / aMax;
             obj.fluxHandle.setTimeStep(obj.dt);
             obj.t = [obj.t, obj.t(end)+obj.dt];
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
            v = obj.U(2,:) ./ obj.U(1,:);
            c = sqrt(obj.gamma * ((obj.gamma - 1) * (obj.U(3,:) - 0.5 * obj.U(2,:).^2 ./ obj.U(1,:))) ./ obj.U(1,:));
            aL = v - c;
            aR = v + c;
        end
        
        function performTimeUpdate(obj)
            % 3rd order SSP-RK method
            U1 = obj.U + obj.dt * obj.dUdt;
            U2 = 0.75 * obj.U + 0.25 * U1 + ...
                 0.25 * obj.dt * obj.calculateRightHandSide(U1);
            obj.U = (obj.U  + 2 * U2 + ...
                     2 * obj.dt * obj.calculateRightHandSide(U2)) / 3;
        end
        
        function assignResults(obj)
            obj.rho = [obj.rho; obj.U(1,:)];
            obj.v = [obj.v; obj.U(2,:) ./ obj.U(1,:)]; 
            obj.p = [obj.p; (obj.gamma - 1) * (obj.U(3,:) - 0.5 * obj.U(2,:).^2 ./ obj.U(1,:))];
            obj.E = [obj.E; obj.U(3,:)];
            obj.c = sqrt(obj.gamma * obj.p ./ obj.rho);
            obj.M = [obj.M; obj.U(2,:) ./ (obj.U(1,:) .* obj.c)];
        end
    end
end