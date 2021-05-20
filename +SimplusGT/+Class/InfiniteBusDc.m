% This class defines the model of a dc infinite bus.

% Notes: An infinite bus is short-circuit in small-signal analysis.

% Author(s): Yitong Li

classdef InfiniteBusDc < SimplusGT.Class.ModelAdvance
    
    methods(Static)
        
        % Set the strings of state, input, output
        function SetString(obj)
            obj.StateString  = {};    	% x
         	obj.InputString  = {'i'};  	% u
        	obj.OutputString = {'v'};  	% y
        end
        
        % Calculate the equilibrium
        function [x_e,u_e,xi] = Equilibrium(obj)
         	% Get the power PowerFlow values
            P 	= obj.PowerFlow(1);
            Q	= obj.PowerFlow(2);
            V	= obj.PowerFlow(3);
            xi	= obj.PowerFlow(4);
            
            % Calculate
            v = V;
            i = P/V;
            
            % Set equilibrium
            u_e = i;
            x_e = [];
        end
        
        % State space model
        function [Output] = StateSpaceEqu(obj,x,u,CallFlag)     
            if CallFlag == 1
                % State equations: dx/dt = f(x,u)
              	f_xu = [];
                Output = f_xu;
            elseif CallFlag == 2
                % Output equations: y = g(x,u)
                v = obj.PowerFlow(3);
                g_xu = v;
                Output = g_xu;              
            end
        end
        
    end
end