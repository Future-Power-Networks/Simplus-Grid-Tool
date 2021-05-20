% This class defines the model of an AC infinite bus.

% Notes: An infinite bus is short-circuit in small-signal analysis.

% Author(s): Yitong Li

classdef InfiniteBusAc < SimplusGT.Class.ModelAdvance
    
    methods(Static)
        
        % Set the strings of state, input, output
        function SetString(obj)
            obj.StateString  = {};                  % x
         	obj.InputString  = {'i_d','i_q'};     	% u
        	obj.OutputString = {'v_d','v_q','w'};  	% y
        end
        
        % Calculate the equilibrium
        function [x_e,u_e,xi] = Equilibrium(obj)
         	% Get the power PowerFlow values
            P 	= obj.PowerFlow(1);
            Q	= obj.PowerFlow(2);
            V	= obj.PowerFlow(3);
            xi	= obj.PowerFlow(4);
            
            % Calculate
            i_d = P/V;
            i_q = -Q/V;
            
            u_e = [i_d; i_q];
            x_e = [];
            xi  = [xi];
        end
        
        % State space model
        function [Output] = StateSpaceEqu(obj,x,u,CallFlag)     
            w	= obj.PowerFlow(5);
            if CallFlag == 1
                % State equation: dx/dt = f(x,u)
              	f_xu = [];
                Output = f_xu;
            elseif CallFlag == 2
                % Output equation: y = g(x,u)
                V	= obj.PowerFlow(3);
                v_d = V;
                v_q = 0;
                g_xu = [v_d; v_q; w];
                Output = g_xu;              
            end
        end
        
    end
end