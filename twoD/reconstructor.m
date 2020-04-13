%**************************************************************************
% Created    : 30.03.2020
% Author     : Andre Breuer
%**************************************************************************

classdef reconstructor < handle
    properties (Access = private)
        dx;         % step width in space
        dy;         % step width in space
        gamma;
        q;          % Exponent for controlling stability and accuracy
    end

    methods (Access = public)
        function obj = reconstructor(dx, dy, gamma, q)
            if nargin < 4
                q = 1.4;
            end
            
            obj.dx = dx;
            obj.dy = dy;
            obj.gamma = gamma;
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
        
        function [UL, UR] = reconstructValuesLDLR(obj, U, direction, BC)
            UGhost = obj.setGhostCells(U, direction, BC);
            [UL, UR] = obj.calculateValuesLDLR(UGhost, direction);
        end
    end
        
    methods (Access = private)
        function UGhost = setGhostCells(obj, U, direction, BC)
            switch direction
                case 'x'
                    UGhost = obj.setGhostCellsX(U, BC);
                case 'y'
                    UGhost = obj.setGhostCellsY(U, BC);
            end
        end
        
        function UGhost = setGhostCellsX(obj, U, BC)
            switch BC
                case 'periodic'
                    % Ghost cells in x-direction
                    UGhost = cat(1, U(end,:,:), U, U(1,:,:));
                case 'transmissive'
                    % Ghost cells in x-direction
                    UGhost = cat(1, U(1,:,:), U, U(end,:,:));
                case 'reflectiveX'
                    % Ghost cells in x-direction
                    UGhost = cat(1, ...
                             cat(3, U(1,:,1), -U(1,:,2), -U(1,:,3), U(1,:,4)), ...
                             U, ...
                             cat(3, U(end,:,1), -U(end,:,2), -U(end,:,3), U(end,:,4)));
                case 'reflectiveFull'
                    UGhost = cat(1, ...
                             cat(3, U(1,:,1), -U(1,:,2), -U(1,:,3), U(1,:,4)), ...
                             U, ...
                             cat(3, U(end,:,1), -U(end,:,2), -U(end,:,3), U(end,:,4)));
                case 'RT'
                    UGhost = cat(1, ...
                             cat(3, U(1,:,1), -U(1,:,2), -U(1,:,3), U(1,:,4)), ...
                             U, ...
                             cat(3, U(end,:,1), -U(end,:,2), -U(end,:,3), U(end,:,4)));
            end
        end
        
        function UGhost = setGhostCellsY(obj, U, BC)
            switch BC
                case 'periodic'
                    % Ghost cells in y-direction
                    UGhost = cat(2, U(:,end,:), U, U(:,1,:));
                case 'transmissive'
                    % Ghost cells in y-direction
                    UGhost = cat(2, U(:,1,:), U, U(:,end,:));
                case 'reflectiveX'
                    % Ghost cells in y-direction
                    UGhost = cat(2, U(:,end,:), U, U(:,1,:));
                case 'reflectiveFull'
                    % Ghost cells in y-direction
                    UGhost = cat(2, ...
                             cat(3, U(:,1,1), -U(:,1,2), -U(:,1,3), U(:,1,4)), ...
                             U, ...
                             cat(3, U(:,end,1), -U(:,end,2), -U(:,end,3), U(:,end,4)));
                case 'RT'
                    UGhost = cat(2, ...
                             cat(3, 2*ones(size(U(:,1,1))), ...
                                    zeros(size(U(:,1,2))), ...
                                    zeros(size(U(:,1,3))), ...
                                    ones(size(U(:,1,4))) / (obj.gamma - 1)), ...
                             U, ...
                             cat(3, ones(size(U(:,end,1))), ...
                                    zeros(size(U(:,end,2))), ...
                                    zeros(size(U(:,end,3))), ...
                                    2.5*ones(size(U(:,end,4))) / (obj.gamma - 1)));  
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
        
        function [UL, UR] = calculateValuesLDLR(obj, UGhost, direction)
            % Lateral derivatives
            [d1, d2] = obj.calculateLateralDerivatives(UGhost, direction);
            
            [a, ~, c, d] = obj.calculateReconstructionCoefficients(d1, d2);
            [etaAP, etaBP, etaAM, etaBM] = obj.etaFunction(a);

            switch direction
                case 'x'
                    UL = UGhost(2:end-1,:,:) + obj.dx * (c .* etaAP + d .* etaBP);
                    UR = UGhost(2:end-1,:,:) + obj.dx * (c .* etaAM + d .* etaBM); 
                case 'y'
                    UL = UGhost(:,2:end-1,:) + obj.dy * (c .* etaAP + d .* etaBP);
                    UR = UGhost(:,2:end-1,:) + obj.dy * (c .* etaAM + d .* etaBM); 
            end
        end
        
        function [d1, d2] = calculateLateralDerivatives(obj, U, direction)
            switch direction
                case 'x'
                    d1 = (U(2:end-1,:,:) - U(1:end-2,:,:)) / obj.dx;
                    d2 = (U(3:end,:,:) - U(2:end-1,:,:)) / obj.dx;
                case 'y'
                    d1 = (U(:,2:end-1,:) - U(:,1:end-2,:)) / obj.dy;
                    d2 = (U(:,3:end,:) - U(:,2:end-1,:)) / obj.dy;
            end
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