%**************************************************************************
% Created    : 31.03.2020
% Author     : Andre Breuer
%**************************************************************************

classdef numericalFluxesEuler2D < handle
     properties (Access = private)
        dx;         % step width in space
        dy;         % step width in space
        dt;         % step width in time
        gamma;      % Exponent for controlling stability and accuracy
    end
    
    
    methods (Access = public)
        function obj = numericalFluxesEuler2D(dx, dy, dt, gamma)
            obj.dx = dx;
            obj.dy = dy;
            obj.dt = dt;
            obj.gamma = gamma;
        end
        
        function setTimeStep(obj, dt)
            obj.dt = dt;
        end
        
        function F = calculateNumericalFlux(obj, UL, UR, fluxFunction, fluxType)
            switch fluxType
                case 'HLL'
                    F = obj.calculateHLLFlux(UL, UR, fluxFunction);
                case 'LF'
                    F = obj.calculateLFFlux(UL, UR, fluxFunction);
                case 'vanLeer'
                    F = obj.calculateVanLeerFlux(UL, UR, fluxFunction);
                otherwise
                    error('Invalid flux function choosen...');
            end
        end
    end
    
    methods (Access = private)
        function F = calculateHLLFlux(obj, UL, UR, fluxFunction)
            [aL, aR] = obj.calculateEigenvalues(UL, UR);
            switch fluxFunction
                case 'F'
                    fL = obj.fluxFunctionF(UL);
                    fR = obj.fluxFunctionF(UR);
                case 'G'
                    fL = obj.fluxFunctionG(UL);
                    fR = obj.fluxFunctionG(UR);
            end
             
            aL = repmat(aL, [1 1 4]);
            aR = repmat(aR, [1 1 4]);
            F = (aR .* fL - aL .* fR + aL .* aR .* (UR - UL)) ./ (aR - aL);
            F = F .* (aL < 0 & aR > 0) + fL .* (aL >= 0) + fR .* (aR <= 0);
        end
        
        function F = calculateVanLeerFlux(obj, UL, UR, fluxFunction)
            % All eigenvalues are negative for  M  < -1
            % All eigenvalues are positive for  M  >  1
            % Eigenvalues have mixed signs for |M| <= 1
            switch fluxFunction
                case 'F'
                    F = obj.vanLeerFluxPlusF(UL) + obj.vanLeerFluxMinusF(UR);
                case 'G'
                    F = obj.vanLeerFluxPlusG(UL) + obj.vanLeerFluxMinusG(UR);
            end
        end
        
        function F = calculateLFFlux(obj, UL, UR, fluxFunction)
            switch fluxFunction
                case 'F'
                    fL = obj.fluxFunctionF(UL);
                    fR = obj.fluxFunctionF(UR);
                case 'G'
                    fL = obj.fluxFunctionG(UL);
                    fR = obj.fluxFunctionG(UR);
            end
            aMax = max(obj.calculateMaxEigenvalue(UL), obj.calculateMaxEigenvalue(UR));
            
            switch fluxFunction
                case 'F'
                    F = 0.5 * (fL + fR - obj.dt / obj.dx * abs(aMax) .* (UR - UL));
                case 'G'
                    F = 0.5 * (fL + fR - obj.dt / obj.dy * abs(aMax) .* (UR - UL));
            end
        end
        
        function aMax = calculateMaxEigenvalue(obj, U)
            v = U(:,:,2) ./ U(:,:,1);
            c = sqrt(obj.gamma * ((obj.gamma - 1) * (U(:,:,4) - 0.5 * (U(:,:,2).^2 + U(:,:,3).^2) ./ U(:,:,1))) ./ U(:,:,1));
            aMax = abs(v) + c;
        end
        
        function [aL, aR] = calculateEigenvalues(obj, lhs, rhs, fluxFunction)            
            [aLLhs, aRLhs] = obj.calculateEigenvalue(lhs);
            [aLRhs, aRRhs] = obj.calculateEigenvalue(rhs);
            
            aL = min(aLLhs, aLRhs);
            aR = max(aRLhs, aRRhs);
        end
        
        function [aL, aR] = calculateEigenvalue(obj, U)
            v = U(:,:,2) ./ U(:,:,1);
            c = sqrt(obj.gamma * ((obj.gamma - 1) * (U(:,:,4) - 0.5 * (U(:,:,2).^2 + U(:,:,3).^2) ./ U(:,:,1))) ./ U(:,:,1));
            aL = v - c;
            aR = v + c;
        end
        
        function F = fluxFunctionF(obj, U)
            p = (obj.gamma - 1) * (U(:,:,4) - 0.5 * ...
                (U(:,:,2).^2 + U(:,:,3).^2) ./ U(:,:,1));
            F = cat(3, U(:,:,2), ...
                       U(:,:,2).^2 ./ U(:,:,1) + p, ...
                       U(:,:,2) .* U(:,:,3) ./ U(:,:,1), ...
                       U(:,:,2) .* (U(:,:,4) + p) ./ U(:,:,1));
        end
        
        function F = fluxFunctionG(obj, U)
            p = (obj.gamma - 1) * (U(:,:,4) - 0.5 * ...
                (U(:,:,2).^2 + U(:,:,3).^2) ./ U(:,:,1));
            F = cat(3, U(:,:,3), ...
                       U(:,:,2) .* U(:,:,3) ./ U(:,:,1), ...
                       U(:,:,3).^2 ./ U(:,:,1) + p, ...
                       U(:,:,3) .* (U(:,:,4) + p) ./ U(:,:,1));
        end
        
        function b = CFLconditionSatisfied(obj, aMax)
             if (aMax * obj.dt / obj.dx) <= 1
                b = true;
                return;
             end
             b = false;
        end
        
        function FP = vanLeerFluxPlusF(obj, U)
            rho = U(:,:,1);
            vX = U(:,:,2) ./ U(:,:,1);
            vY = U(:,:,3) ./ U(:,:,1);
            p = (obj.gamma - 1) * (U(:,:,4) - 0.5 * (vX.^2 + vY.^2));
            c = sqrt(obj.gamma * p ./ U(:,:,1));
            MX = vX ./ c;
            
            fMass = 0.25 * rho .* c .* (1 + MX).^2;
            fMomentumX = fMass .* ((obj.gamma - 1) .* vX + 2*c) / obj.gamma;
            %fMomentumX = c .* fMass .* ((obj.gamma-1)*MX + 2) / obj.gamma;
            fMomentumY = fMass .* vY;
            fEnergy = fMass .* (((obj.gamma - 1)*vX + 2*c).^2 /(2*(obj.gamma^2 -1)) ...
                      + 0.5 * (vX.^2 + vY.^2));
%             fEnergy = 0.5 * fMass .* (c.^2 .* ((obj.gamma-1)*MX + 2).^2 ./ ...
%                       (obj.gamma^2 -1) + (U(:,:,2).^2 + U(:,:,3).^2) ./ (U(:,:,1).^2));
            FP = cat(3, fMass, fMomentumX, fMomentumY, fEnergy);
            Mx4 = repmat(MX, [1 1 4]);
            FReal = obj.fluxFunctionF(U);
            % Since the shape should be conserved it's a rather unreadable
            % way to compute the flux. It uses logical 3D-matrices to
            % filter out the desired values
            FP = (FP .* (Mx4 < 1) + FReal .* (Mx4 >= 1)) .* not(Mx4 <= -1);
        end
        
        function FM = vanLeerFluxMinusF(obj, U)
            rho = U(:,:,1);
            vX = U(:,:,2) ./ U(:,:,1);
            vY = U(:,:,3) ./ U(:,:,1);
            p = (obj.gamma - 1) * (U(:,:,4) - 0.5 * (vX.^2 + vY.^2));
            c = sqrt(obj.gamma * p ./ U(:,:,1));
            MX = vX ./ c;
            
            fMass = - 0.25 * rho .* c .* (1 - MX).^2;
            fMomentumX = fMass .* ((obj.gamma - 1) .* vX - 2*c) / obj.gamma;
            %fMomentumX = c .* fMass .* ((obj.gamma-1)*MX - 2) / obj.gamma;
            fMomentumY = fMass .* vY;
            fEnergy = fMass .* (((obj.gamma - 1)*vX - 2*c).^2 /(2*(obj.gamma^2 -1)) ...
                      + 0.5 * (vX.^2 + vY.^2));
%             fEnergy = 0.5 * fMass .* (c.^2 .* ((obj.gamma-1)*MX - 2).^2 ./ ...
%                       (obj.gamma^2 -1) + (U(:,:,2).^2 + U(:,:,3).^2) ./ (U(:,:,1).^2));
            FM = cat(3, fMass, fMomentumX, fMomentumY, fEnergy);
            Mx4 = repmat(MX, [1 1 4]);
            FReal = obj.fluxFunctionF(U);
            % Since the shape should be conserved it's a rather unreadable
            % way to compute the flux. It uses logical 3D-matrices to
            % filter out the desired values
            FM = (FM .* (Mx4 > -1) + FReal .* (Mx4 <= -1)) .* not(Mx4 >= 1);
        end
        
        function GP = vanLeerFluxPlusG(obj, U)
            rho = U(:,:,1);
            vX = U(:,:,2) ./ U(:,:,1);
            vY = U(:,:,3) ./ U(:,:,1);
            p = (obj.gamma - 1) * (U(:,:,4) - 0.5 * (vX.^2 + vY.^2));
            c = sqrt(obj.gamma * p ./ U(:,:,1));
            MY = vY ./ c;
            
            fMass = 0.25 * rho .* c .* (1 + MY).^2;
            fMomentumX = fMass .* vX;
            fMomentumY = fMass .* ((obj.gamma - 1) .* vY + 2*c) / obj.gamma;
%             fMomentumY = c .* fMass .* ((obj.gamma-1)*MY + 2) / obj.gamma;
            fEnergy = fMass .* (((obj.gamma - 1)*vY + 2*c).^2 /(2*(obj.gamma^2 -1)) ...
                      + 0.5 * (vX.^2 + vY.^2));
%             fEnergy = 0.5 * fMass .* (c.^2 .* ((obj.gamma-1)*MY + 2).^2 ./ ...
%                       (obj.gamma^2 -1) + (U(:,:,2).^2 + U(:,:,3).^2) ./ (U(:,:,1).^2));
            GP = cat(3, fMass, fMomentumX, fMomentumY, fEnergy);
            My4 = repmat(MY, [1 1 4]);
            FReal = obj.fluxFunctionG(U);
            % Since the shape should be conserved it's a rather unreadable
            % way to compute the flux. It uses logical 3D-matrices to
            % filter out the desired values
            GP = (GP .* (My4 < 1) + FReal .* (My4 >= 1)) .* not(My4 <= -1);
        end
        
        function GM = vanLeerFluxMinusG(obj, U)
            rho = U(:,:,1);
            vX = U(:,:,2) ./ U(:,:,1);
            vY = U(:,:,3) ./ U(:,:,1);
            p = (obj.gamma - 1) * (U(:,:,4) - 0.5 * (vX.^2 + vY.^2));
            c = sqrt(obj.gamma * p ./ U(:,:,1));
            MY = vY ./ c;
            
            fMass = - 0.25 * rho .* c .* (1 - MY).^2;
            fMomentumX = fMass .* vX;
            fMomentumY = fMass .* ((obj.gamma - 1) .* vY - 2*c) / obj.gamma;
            fEnergy = fMass .* (((obj.gamma - 1)*vY - 2*c).^2 /(2*(obj.gamma^2 -1)) ...
                      + 0.5 * (vX.^2 + vY.^2));
%             fMomentumY = c .* fMass .* ((obj.gamma-1)*MY - 2) / obj.gamma;
%             fEnergy = 0.5 * fMass .* (c.^2 .* ((obj.gamma-1)*MY - 2).^2 ./ ...
%                       (obj.gamma^2 -1) + (U(:,:,2).^2 + U(:,:,3).^2) ./ (U(:,:,1).^2));
            GM = cat(3, fMass, fMomentumX, fMomentumY, fEnergy);
            My4 = repmat(MY, [1 1 4]);
            FReal = obj.fluxFunctionG(U);
            % Since the shape should be conserved it's a rather unreadable
            % way to compute the flux. It uses logical 3D-matrices to
            % filter out the desired values
            GM = (GM .* (My4 > -1) + FReal .* (My4 <= -1)) .* not(My4 >= 1);
        end
        
        function b = CFLconditionVanLeerSatisfied(obj, U, c, M)
            vMaxAbs = abs(max(U(2,:) ./ U(1,:)));
            CFLmax = (2 * obj.gamma + abs(max(M)) * (3 - obj.gamma)) ./ ...
                      (obj.gamma + 3);
            if (obj.dt / obj.dx * (vMaxAbs + max(c))) <= CFLmax
               b = true;
               return;
            end
            b = false; 
        end 
    end
end