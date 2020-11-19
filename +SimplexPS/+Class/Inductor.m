% This class defines the model of a single phase inductor for test.

% Author(s): Yitong Li

classdef Inductor < SimplexPS.Class.ModelAdvance
    methods(Static)
        % Set the strings of input, output, state
        function SetString(obj)
            obj.StateString  = {'i'};        % x
            obj.InputString  = {'v'};        % u
            obj.OutputString = {'i'};        % y
        end
        
        % Calculate the equilibrium
        function Equilibrium(obj)
            obj.x_e = 0;
            obj.u_e = 0;
            obj.xi = 0;
        end
        
        % State space model
        function [Output] = StateSpaceEqu(obj,x,u,CallFlag)
            % Get input
            v = u(1);
            % Get state
            i = x(1);
            % Get parameter
            L = obj.Para;
            % State space equations
            if CallFlag == 1
                % State equations: dx/dt = f(x,u)
                di = v/L; 
                f_xu = [di]; 
                Output = f_xu;
            elseif CallFlag == 2
                % Output equations: y = g(x,u)
                g_xu = [i]; 
                Output = g_xu;              
            end
        end
        
    end
end