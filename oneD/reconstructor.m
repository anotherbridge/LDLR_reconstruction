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
            UGhost = obj.setGhostCells(U, BC);
            [UL, UR] = obj.calculateValuesLDLR(UGhost);
        end
    end
        
    methods (Access = private)
        function UGhost = setGhostCells(obj, U, BC)
            switch BC
                case 'periodic'
                    UGhost = [U(:,end), U, U(:,1)];
                case 'transmissive'
                    UGhost = [U(:,1), U, U(:,end)];
                case 'reflectiveRight'
                    UGhost = [U(:,1), U, [U(1,end); -U(2,end); U(3,end)]];
                case 'reflectiveFull'
                    UGhost = [[U(1,1); -U(2,1); U(3,1)], U, [U(1,end); -U(2,end); U(3,end)]];
            end
        end
        
        function Phi = minModLimiter(obj, UL, UR)
            Theta = UL ./ UR;
            Phi = max(0, min(1, Theta));
        end
        
        function Phi = vanLeerLimiter(obj, UL, UR)
            Theta = UL ./ UR;
            Phi = (Theta + abs(Theta)) ./ (1 + Theta);
        end
        
        function [UL, UR] = calculateValuesLDLR(obj, UGhost)
            % Lateral derivatives
            d1 = (UGhost(:,2:end-1) - UGhost(:,1:end-2)) / obj.dx;
            d2 = (UGhost(:,3:end) - UGhost(:,2:end-1)) / obj.dx;
            
            [a, ~, c, d] = obj.calculateReconstructionCoefficients(d1, d2);
            [etaAP, etaBP, etaAM, etaBM] = obj.etaFunction(a);
            
            UL = UGhost(:,2:end-1) + obj.dx * (c .* etaAP + d .* etaBP);
            UR = UGhost(:,2:end-1) + obj.dx * (c .* etaAM + d .* etaBM); 
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