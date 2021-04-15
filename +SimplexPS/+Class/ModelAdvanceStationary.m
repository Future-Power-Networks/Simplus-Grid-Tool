% This class defines advanced properties and methods for state space models
% in stationary reference frame used in script and simulink.

% Author(s): Yitong Li, Yunjie Gu

%% Notes
%
% The models in this class achieves the transformation between dq and abc
% in functions "updateImpl" and "outputImpl".
%
% This class makes the model fit in both simulation (simulink) and
% theoratical analysis (script).
%
% Subclass of this class contains the specific models of different devices.

%% Class
classdef ModelAdvanceStationary < SimplexPS.Class.ModelAdvance
    
methods(Access = protected)
    
    % ### Update discreate states
    function updateImpl(obj, u)
        
        % Convert electrical u_abc to u_dq
    	theta = obj.x(end);
        u_dq = SimplexPS.abc2dq(u(1:3),theta);
        u_ = [u_dq;u(4:end)];
        
        switch obj.DiscreMethod
            
            % ### Case 1: Forward Euler 
          	case 1
                delta_x = obj.Ts * obj.StateSpaceEqu(obj, obj.x, u_, 1);
                obj.x = delta_x + obj.x;
                
            % ### Case 2 : Hybrid Euler-Trapezoidal (Yunjie's Method)
            % ### Case 2' : General virtual dissipation (Yitong's Method)
            case 2               
                if obj.LinearizationTimes == 2    
                    obj.SetDynamicSS(obj,obj.x,obj.u);
                end
                delta_x = obj.Wk * obj.StateSpaceEqu(obj,obj.x,u_,1);               
                obj.x = obj.x + delta_x;             
        end
        
        obj.Timer = obj.Timer + obj.Ts;
    end
        
    % ### Calculate output y
	function y = outputImpl(obj,u)
        
         % Convert electrical u_abc to u_dq
    	theta = obj.x(end);
        u_dq = SimplexPS.abc2dq(u(1:3),theta);
        u_ = [u_dq;u(4:end)];
        
        switch obj.DiscreMethod
            
            % ### Case 1: Forward Euler
          	case 1
                if obj.DirectFeedthrough
                    y_ = obj.StateSpaceEqu(obj,obj.x,u_,2);
                else
                    if obj.VirtualResistor
                        y_ = obj.StateSpaceEqu(obj,obj.x,obj.uk,2) - obj.Gk*obj.uk + obj.Fk*(u_-obj.uk);
                    else
                        y_ = obj.StateSpaceEqu(obj,obj.x,obj.uk,2) + obj.Fk*(u_-obj.uk);
                    end
                end
                
            % ### Case 2 : Hybrid Euler-Trapezoidal (Yunjie's Method)
            % ### Case 2' : General virtual dissipation (Yitong's Method)
            case 2
                
                if obj.DirectFeedthrough
                    y_ = obj.StateSpaceEqu(obj,obj.x,u_,2) + 1/2*obj.Ck*obj.Wk*obj.StateSpaceEqu(obj,obj.x,u_,1);
                else
                    if obj.VirtualResistor
                        y_ = obj.StateSpaceEqu(obj,obj.x,obj.uk,2) + 1/2*obj.Ck*obj.Wk*obj.StateSpaceEqu(obj,obj.x,obj.uk,1) - obj.Gk*obj.uk + obj.Fk*(u_-obj.uk);
                    else
                        y_ = obj.StateSpaceEqu(obj,obj.x,obj.uk,2) + 1/2*obj.Ck*obj.Wk*obj.StateSpaceEqu(obj,obj.x,obj.uk,1) + obj.Fk*(u_-obj.uk);
                    end
                end
        end
        
        obj.uk = u_;                 % store the current u=u[k+1/2]
        obj.Start = 1;
        
        % Convert electrical y_dq to y_abc
        y_abc = SimplexPS.dq2abc(y_(1:2),theta);
        y = [y_abc;y_(3:end)];
                            
    end
    
 	function [size] = getOutputSizeImpl(obj)
        [~,~,Output] = obj.SignalList(obj);
        size = length(Output)+1;
    end
    
end
    
end     % End class definition