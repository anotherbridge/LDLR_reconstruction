%**************************************************************************
% Created    : 31.03.2020
% Author     : Andre Breuer
%**************************************************************************

classdef numericalFluxesEuler < handle
     properties (Access = private)
        dx;         % step width in space
        dt;         % step width in time
        gamma;      % Exponent for controlling stability and accuracy
    end
    
    
    methods (Access = public)
        function obj = numericalFluxesEuler(dx, dt, gamma)
            obj.dx = dx;
            obj.dt = dt;
            obj.gamma = gamma;
        end
        
        function setTimeStep(obj, dt)
            obj.dt = dt;
        end
        
        function F = calculateNumericalFlux(obj, UL, UR, fluxType)
            switch fluxType
                case 'HLL'
                    F = obj.calculateHLLFlux(UL, UR);
                case 'LF'
                    F = obj.calculateLFFlux(UL, UR);
                case 'vanLeer'
                    F = obj.calculateVanLeerFlux(UL, UR);
                otherwise
                    error('Invalid flux function choosen...');
            end
        end
    end
    
    methods (Access = private)
        function F = calculateHLLFlux(obj, UL, UR)
            [aL, aR] = obj.calculateEigenvalues(UL, UR);
            fL = obj.fluxFunction(UL);
            fR = obj.fluxFunction(UR);
             
            aL = repmat(aL, 3, 1);
            aR = repmat(aR, 3, 1);
            F = (aR .* fL - aL .* fR + aL .* aR .* (UR - UL)) ./ (aR - aL);             
            F(aL >= 0) = fL(aL >= 0);
            F(aR <= 0) = fR(aR <= 0);
        end
        
        function F = calculateVanLeerFlux(obj, UL, UR)
            % All eigenvalues are negative for  M  < -1
            % All eigenvalues are positive for  M  >  1
            % Eigenvalues have mixed signs for |M| <= 1

            F = obj.vanLeerFluxPlus(UL) + obj.vanLeerFluxMinus(UR);
        end
        
        function F = calculateLFFlux(obj, UL, UR)
            fL = obj.fluxFunction(UL);
            fR = obj.fluxFunction(UR);
            aMax = max(obj.calculateMaxEigenvalue(UL), obj.calculateMaxEigenvalue(UR));

            F = 0.5 * (fL + fR - obj.dt / obj.dx * abs(aMax) .* (UR - UL));
        end
        
        function aMax = calculateMaxEigenvalue(obj, U)
            v = U(2,:) ./ U(1,:);
            c = sqrt(obj.gamma * ((obj.gamma - 1) * (U(3,:) - 0.5 * U(2,:).^2 ./ U(1,:))) ./ U(1,:));
            aMax = abs(v) + c;
        end
        
        function [aL, aR] = calculateEigenvalues(obj, lhs, rhs)            
            [aLLhs, aRLhs] = obj.calculateEigenvalue(lhs);
            [aLRhs, aRRhs] = obj.calculateEigenvalue(rhs);
            
            aL = min(aLLhs, aLRhs);
            aR = max(aRLhs, aRRhs);
        end
        
        function [aL, aR] = calculateEigenvalue(obj, U)
            v = U(2,:) ./ U(1,:);
            c = sqrt(obj.gamma * ((obj.gamma - 1) * (U(3,:) - 0.5 * U(2,:).^2 ./ U(1,:))) ./ U(1,:));
            aL = v - c;
            aR = v + c;
        end
        
        function F = fluxFunction(obj, U)
            F = [U(2,:);
                 0.5 * (3 - obj.gamma) * U(2,:).^2 ./ U(1,:) + (obj.gamma - 1) * U(3,:);
                 U(2,:) ./ U(1,:) .* (obj.gamma * U(3,:) - 0.5 * (obj.gamma - 1) * U(2,:).^2 ./ U(1,:))];
        end
        
        function b = CFLconditionSatisfied(obj, aMax)
             if (aMax * obj.dt / obj.dx) <= 1
                b = true;
                return;
             end
             b = false;
        end
        
        function FP = vanLeerFluxPlus(obj, U)
            rho = U(1,:);
            c = sqrt(obj.gamma * ((obj.gamma - 1) * (U(3,:) - 0.5 * U(2,:).^2 ./ U(1,:))) ./ U(1,:)); 
            M = U(2,:) ./ (U(1,:) .* c);
            
            fMass = 0.25 * rho .* c .* (1 + M).^2;
            fMomentum = c .* fMass .* ((obj.gamma-1)*M + 2) / obj.gamma;
            fEnergy = 0.5 * obj.gamma^2 * fMomentum.^2 ./ ...
                      ((obj.gamma^2 - 1) * fMass);
            FP = [fMass; fMomentum; fEnergy];
            FP(:,M >= 1) = obj.fluxFunction(U(:,M >= 1));
            FP(:,M <= -1) = 0;
        end
        
        function FM = vanLeerFluxMinus(obj, U)
            rho = U(1,:);
            c = sqrt(obj.gamma * ((obj.gamma - 1) * (U(3,:) - 0.5 * U(2,:).^2 ./ U(1,:))) ./ U(1,:)); 
            M = U(2,:) ./ (U(1,:) .* c);
            
            fMass = - 0.25 * rho .* c .* (1 - M).^2;
            fMomentum = c .* fMass .* ((obj.gamma-1)*M - 2) / obj.gamma;
            fEnergy = 0.5 * obj.gamma^2 * fMomentum.^2 ./ ...
                      ((obj.gamma^2 - 1) * fMass);
            FM = [fMass; fMomentum; fEnergy];
            FM(:,M <= -1) = obj.fluxFunction(U(:,M <= -1));
            FM(:,M >= 1) = 0;
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