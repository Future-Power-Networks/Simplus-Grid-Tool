% This class defines the model of a single phase inductor for test.

% Author(s): Yitong Li

classdef Inductor < SimplexPS.Class.ModelAdvance
    
    methods
        % Constructor
        function obj = Inductor(varargin)
            setProperties(obj,nargin,varargin{:});
        end
    end
    
    methods(Static)
        % Set the strings of input, output, state
        function [State,Input,Output] = SignalList(obj)
            StateString  = {'i'};        % x
            InputString  = {'v'};        % u
            OutputString = {'i'};        % y
        end
        
        % Calculate the equilibrium
        function [x_e,u_e,xi] = Equilibrium(obj)
            % Calculate the equilibrium
            x_e = 0;
            u_e = 0;
            xi = 0;
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
         	% dx/dt = f(x,u)
            % y     = g(x,u)
            if CallFlag == 1
                % ### Call state equation: dx/dt = f(x,u)
                di = v/L; 
                f_xu = [di]; 
                Output = f_xu;
                
            elseif CallFlag == 2
                % ### Call output equation: y = g(x,u)
                g_xu = [i]; 
                Output = g_xu;              
            end
        end
        
    end
end