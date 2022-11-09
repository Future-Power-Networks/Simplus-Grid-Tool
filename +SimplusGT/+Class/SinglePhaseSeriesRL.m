% Author(s): Yitong Li

%% Class

classdef SinglePhaseSeriesRL < SimplusGT.Class.ModelAdvanceNetwork
    
    properties
        R;          % Resistor
        L;          % Inductor
    end
    
    methods
        % Constructor
        function obj =SinglePhaseSeriesRL(varargin)
            setProperties(obj,nargin,varargin{:});
        end
    end
    
    methods(Static)
        % Set the strings of input, output, state
        % For an inductor, voltage v is input, current i is output, current
        % i is state.
        function [State,Input,Output] = SignalList(obj)
            State  = {'i'};        	% x
            Input  = {'v'};         % u
            Output = {'i'};       	% y
        end
        
        % Calculate the equilibrium
        % For simplicity, we set the initial v and i of this inductor to i.
        % xi is the initial phase angle for three-phase apparatus and
        % calculated by power flow, which is not useful here.
        function [x_e,u_e] = Equilibrium(obj)
            % Calculate the equilibrium
            x_e = [0];
            u_e = [0];
        end
        
        % State space model
        function [Output] = StateSpaceEqu(obj,x,u,CallFlag)
            % Get input
            v = u(1);
            
            % Get state
            i = x(1);
            
            % Get parameter
            L = obj.L;
            R = obj.R;
            
            % State space equations
         	% dx/dt = f(x,u)
            % y     = g(x,u)
            if CallFlag == 1
                % ### Call state equation: dx/dt = f(x,u)
                di = (v-i*R)/L;
                
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