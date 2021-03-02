% This class defines the model of a dc infinite bus.

% Notes: An infinite bus is short-circuit in small-signal analysis.

% Author(s): Yitong Li

classdef InfiniteBusDc < SimplexPS.Class.ModelAdvance
    
    methods(Static)
        
        % Set the strings of state, input, output
        function SetString(obj)
            obj.StateString  = {};                  % x
         	obj.InputString  = {'i'};     	% u
        	obj.OutputString = {'v','w'};  	% y
        end
        
        % Calculate the equilibrium
        function [x_e,u_e,xi] = Equilibrium(obj)
         	% Get the power PowerFlow values
            P 	= obj.PowerFlow(1);
            Q	= obj.PowerFlow(2);
            V	= obj.PowerFlow(3);
            xi	= obj.PowerFlow(4);
            w   = obj.PowerFlow(5);
            
            % Calculate
            v = V;
            i = P/V;
            
            u_e = [i];
            x_e = [];
            xi  = [xi];
        end
        
        % State space model
        function [Output] = StateSpaceEqu(obj,x,u,CallFlag)     
            w= 0;
            if CallFlag == 1
              	f_xu = [];
                Output = f_xu;
            elseif CallFlag == 2
                % Output equations: y = g(x,u)
                V	= obj.PowerFlow(3);
                v = V;
                g_xu = [v; w];
                Output = g_xu;              
            end
        end
        
    end
end