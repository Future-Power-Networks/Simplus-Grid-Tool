% Author(s): Yitong Li


%% Class

classdef ThreePhaseStationarySeriesRL < SimplusGT.Class.ModelAdvanceNetwork
    
    methods
        % Constructor
        function obj = ThreePhaseSeriesRLC(varargin)
            setProperties(obj,nargin,varargin{:});
        end
    end
    
    methods(Static)
        % Set the strings of input, output, state
        % For an inductor, voltage v is input, current i is output, current
        % i is state.
        function [State,Input,Output] = SignalList(obj)
            State  = {};        % x
            Input  = {};        % u
            Output = {};        % y
        end
        
        % Calculate the equilibrium
        % For simplicity, we set the initial v and i of this inductor to i.
        % xi is the initial phase angle for three-phase apparatus and
        % calculated by power flow, which is not useful here.
        function [x_e,u_e] = Equilibrium(obj)
            % Calculate the equilibrium
            x_e = 0;
            u_e = 0;
        end
        
        % State space model
        function [Output] = StateSpaceEqu(obj,x,u,CallFlag)
            % Get input
            va1 = u(1);
            vb1 = u(2);
            vc1 = u(3);
            va2 = u(4);
            vb2 = u(5);
            vc2 = u(6);
            % Get state
            i = x(1);
            % Get parameter
            L = obj.Para;
            
            % State space equations
         	% dx/dt = f(x,u)
            % y     = g(x,u)
            if CallFlag == 1
                % ### Call state equation: dx/dt = f(x,u)
                di = v/L;       % di means di/dt
                                % The inductor is in load convention.
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