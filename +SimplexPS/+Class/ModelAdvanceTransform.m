% This class defines advanced properties and methods for state space models
% in stationary reference frame used in script and simulink.

% Author(s): Yitong Li, Yunjie Gu

%% Notes
%
% The models in this class achieves the transformation between dq and abc
% in functions "updateImpl" and "outputImpl", i.e., reloading these two
% functions in its parent class ModelAdvance.
%
% This class makes the model fit in both simulation (simulink) and
% theoratical analysis (script).
%
% Subclass of this class contains the specific models of different devices.

%% Class
classdef ModelAdvanceTransform < SimplexPS.Class.ModelAdvance
    
methods(Access = protected)
    
    % ### Update discreate states
    function updateImpl(obj, u)
        
        % % Convert electrical u_abc to u_dq
    	theta = obj.x(end);
        % u_dq = SimplexPS.abc2dq(u(1:3),theta);
        % u_ = [u_dq;u(4:end)];

        % % Update in dq frame
        % updateImpl@SimplexPS.Class.ModelAdvance(obj,u_);
        updateImpl@SimplexPS.Class.ModelAdvance(obj,u);
        
    end
        
    % ### Calculate output y
    function y = outputImpl(obj,u)
    
        % Convert electrical u_abc to u_dq
        theta = obj.x(end);
        % u_dq = SimplexPS.abc2dq(u(1:3),theta);
        u_dq = u(1:2);
        %u_ = [u_dq;u(4:end)];

        % Output in dq frame
        y_ = outputImpl@SimplexPS.Class.ModelAdvance(obj,u);
        
        % Convert electrical y_dq to y_abc
        y_abc = SimplexPS.dq2abc(y_(1:2),theta);
        y = [u_dq;y_abc;y_];
                            
    end
    
end
    
end     % End class definition