% This class defines the model of a single phase inductor for test.

% Author(s): Yitong Li

classdef Class_FloatingBus < Class_Model_Advance
    
    methods(Static)
        
        % Set the strings of input, output, state
        function SetString(obj)
         	obj.InputString  = {'v_d','v_q'};  	% u
        	obj.OutputString = {'i_d','i_q'};  	% y
        	obj.StateString  = {};           	% x
        end
        
        % Calculate the equilibrium
        function Equilibrium(obj)
         	% Get the power PowerFlow values
            P 	= obj.PowerFlow(1);
            Q	= obj.PowerFlow(2);
            V	= obj.PowerFlow(3);
            xi	= obj.PowerFlow(4);
            w   = obj.PowerFlow(5);
            
            % Calculate
            v_d = V;
            v_q = 0;
            
            obj.u_e = [v_d; v_q];
            obj.x_e = [];
            obj.xi = [xi];
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