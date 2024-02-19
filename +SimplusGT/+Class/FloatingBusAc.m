% This class defines the model of an AC floating bus, i.e., open circuit
% bus or no apparatus connected to this bus.

% Author(s): Yitong Li

classdef FloatingBusAc < SimplusGT.Class.ModelAdvance
    
    methods(Static)
        
        % Set the strings of input, output, state
        function SetString(obj)
         	obj.InputString  = {'v_d','v_q'};  	% u
        	obj.OutputString = {'i_d','i_q'};  	% y
        	obj.StateString  = {};           	% x
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
            v_d = V;
            v_q = 0;
            
            u_e = [v_d; v_q];
            x_e = [];
        end
        
        % State space model
        function [Output] = StateSpaceEqu(obj,x,u,CallFlag)        
            if CallFlag == 1
              	f_xu = [];
                Output = f_xu;
            elseif CallFlag == 2
                % Output equations: y = g(x,u)
                i_d = 0;
                i_q = 0;
                g_xu = [i_d; i_q];
                Output = g_xu;              
            end
        end
        
    end
end