%**************************************************************************
% Created    : 30.03.2020
% Author     : Andre Breuer
%**************************************************************************

classdef reconstructor < handle
    properties (Access = private)
        dx;         % step width in space
        q;          % Exponent for controlling stability and accuracy
    end

    methods (Access = public)
        function obj = reconstructor(dx, q)
            if nargin < 2
                q = 1.4;
            end
            
            obj.dx = dx;
            obj.q = q;
        end
        
        function setExponentQ(obj, q)
            obj.q = q; 
        end
        
        function [UL, UR] = reconstructValuesLinear(obj, U, limiter)
            switch limiter
                case 'minMod'
                    Phi = obj.minModLimiter(U - circshift(U,1,2), ...
                          circshift(U,-1,2) - U);
                case 'vanLeer'
                    Phi = obj.vanLeerLimiter(U - circshift(U,1,2), ...
                          circshift(U,-1,2) - U);
            end
            sigma = 0.5 *(circshift(U,-1,2) - circshift(U,1,2)) .* Phi / obj.dx;
            UL = U + 0.5 * sigma * obj.dx;
            UR = U - 0.5 * sigma * obj.dx;
        end
        
        function [UL, UR] = reconstructValuesLDLR(obj, U, BC)
            switch BC
                case 'periodic'
                    [UL, UR] = obj.reconstructValuesPeriodicLDLR(U);
                case 'transmissive'
                    [UL, UR] = obj.reconstructValuesTransmissiveLDLR(U);
                case 'reflectiveRight'
                    [UL, UR] = obj.reconstructValuesReflectiveRightLDLR(U);
                case 'reflectiveFull'
                    [UL, UR] = obj.reconstructValuesReflectiveFullLDLR(U);
            end
        end
    end
        
    methods (Access = private)
        function Phi = minModLimiter(obj, UL, UR)
            Theta = UL ./ UR;
            Phi = max(0, min(1, Theta));
        end
        
        function Phi = vanLeerLimiter(obj, UL, UR)
            Theta = UL ./ UR;
            Phi = (Theta + abs(Theta)) ./ (1 + Theta);
        end
        
        function [UL, UR] = reconstructValuesTransmissiveLDLR(obj, U)
            % Lateral derivatives
            d1 = (U - [U(:,1), U(:,1:end-1)]) / obj.dx;
            d2 = ([U(:,2:end), U(:,end)] - U) / obj.dx;
            
            [a, ~, c, d] = obj.calculateReconstructionCoefficients(d1, d2);
            [etaAP, etaBP, etaAM, etaBM] = obj.etaFunction(a);
            
            UL = U + obj.dx * (c .* etaAP + d .* etaBP);
            UR = U + obj.dx * (c .* etaAM + d .* etaBM);
        end
        
        function [UL, UR] = reconstructValuesReflectiveRightLDLR(obj, U)
            % Lateral derivatives
            d1 = (U - [U(:,1), U(:,1:end-1)]) / obj.dx;
            d2 = ([U(:,2:end), [U(1,end); -U(2,end); U(3,end)]] - U) / obj.dx;
            
            [a, ~, c, d] = obj.calculateReconstructionCoefficients(d1, d2);
            [etaAP, etaBP, etaAM, etaBM] = obj.etaFunction(a);
            
            UL = U + obj.dx * (c .* etaAP + d .* etaBP);
            UR = U + obj.dx * (c .* etaAM + d .* etaBM);
        end
        
        function [UL, UR] = reconstructValuesReflectiveFullLDLR(obj, U)
            % Lateral derivatives
            d1 = (U - [[U(1,1); -U(2,1); U(3,1)], U(:,1:end-1)]) / obj.dx;
            d2 = ([U(:,2:end), [U(1,end); -U(2,end); U(3,end)]] - U) / obj.dx;
            
            [a, ~, c, d] = obj.calculateReconstructionCoefficients(d1, d2);
            [etaAP, etaBP, etaAM, etaBM] = obj.etaFunction(a);
            
            UL = U + obj.dx * (c .* etaAP + d .* etaBP);
            UR = U + obj.dx * (c .* etaAM + d .* etaBM);
        end
        
        function [UL, UR] = reconstructValuesPeriodicLDLR(obj, U)
            % Lateral derivatives
            UM1 = circshift(U,1,2);
            %UM1(:,1) = UM1(:,end);
            UP1 = circshift(U,-1,2);
            %UP1(:,end) = UP1(:,1);
            d1 = (U - UM1) / obj.dx;
            d2 = (UP1 - U) / obj.dx;
            
            [a, ~, c, d] = obj.calculateReconstructionCoefficients(d1, d2);
            [etaAP, etaBP, etaAM, etaBM] = obj.etaFunction(a);

            UL = U + obj.dx * (c .* etaAP + d .* etaBP);
            UR = U + obj.dx * (c .* etaAM + d .* etaBM);
        end
        
        function [a, b, c, d] = calculateReconstructionCoefficients(obj, d1, d2)
             cTol = 0.1 * obj.dx^obj.q;
             a = (1 - cTol) * (1 + cTol - ...
                 (2 * abs(d1).^obj.q .* abs(d2).^obj.q + cTol) ./ ...
                 (abs(d1).^(2*obj.q) + abs(d2).^(2*obj.q) + cTol));
             b = a ./ (a - 1);
             c = (a - 1) .* (d2 .* (1 - b) - d1) ./ (b - a);
             d = d1 - c;
        end
        
        function [etaAP, etaBP, etaAM, etaBM] = etaFunction(obj, a)
            etaAP = obj.etaFunctionPlus(a);
            etaAM = obj.etaFunctionMinus(a);
            etaBP = (a - 1) .* etaAM;
            etaBM = (a - 1) .* etaAP;
        end
        
        function eta = etaFunctionPlus(obj, t)
            eta = - (log(1 - t) + t) ./ (t.^2); 
            eta(t == 0) = 0.5;
        end
        
        function eta = etaFunctionMinus(obj, t)
            eta = ((t - 1).*log(1 - t) - t) ./ (t.^2);
            eta(t == 0) = - 0.5;
        end
    end
end