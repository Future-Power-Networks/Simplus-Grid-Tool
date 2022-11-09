% Author(s): Yitong Li

%% Class

classdef SinglePhaseParallelGC < SimplusGT.Class.ModelAdvanceNetwork
    
    properties
        Conductor;          % Conductor
        Capacitor;          % Capacitor
    end
    
    methods
        % Constructor
        function obj = SinglePhaseParallelGC(varargin)
            setProperties(obj,nargin,varargin{:});
        end
    end
    
    methods(Static)
        % Set the strings of input, output, state
        % For an inductor, voltage v is input, current i is output, current
        % i is state.
        function [State,Input,Output] = SignalList(obj)
            State  = {'v'};       	% x
            Input  = {'i'};         % u
            Output = {'v'};       	% y
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
            i = u(1);

            % Get state
            v = x(1);
            
            % Get parameter
           	C = obj.Capacitor;
            G = obj.Conductor;
            
            % State space equations
         	% dx/dt = f(x,u)
            % y     = g(x,u)
            if CallFlag == 1
                % ### Call state equation: dx/dt = f(x,u)
                dv = (i - v*G)/C;
                
                f_xu = [dv]; 
                Output = f_xu;
                
            elseif CallFlag == 2
                % ### Call output equation: y = g(x,u)
                g_xu = [v]; 
                Output = g_xu;              
            end
            
        end
        
    end
end