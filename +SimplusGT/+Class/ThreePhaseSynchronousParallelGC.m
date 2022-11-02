% Author(s): Yitong Li

%% Class

classdef ThreePhaseSynchronousParallelGC < SimplusGT.Class.ModelAdvanceNetwork
    
    properties
        Conductor;          % Conductor
        Capacitor;          % Capacitor
        w;                  % Frequency
    end
    
    methods
        % Constructor
        function obj = ThreePhaseSynchronousParallelGC(varargin)
            setProperties(obj,nargin,varargin{:});
        end
    end
    
    methods(Static)
        % Set the strings of input, output, state
        % For an inductor, voltage v is input, current i is output, current
        % i is state.
        function [State,Input,Output] = SignalList(obj)
            State  = {'vd','vq','theta'};   	% x
            Input  = {'ia','ib','ic'};          % u
            Output = {'va','vb','vc'};        	% y
        end
        
        % Calculate the equilibrium
        % For simplicity, we set the initial v and i of this inductor to i.
        % xi is the initial phase angle for three-phase apparatus and
        % calculated by power flow, which is not useful here.
        function [x_e,u_e] = Equilibrium(obj)
            % Calculate the equilibrium
            x_e = [0;0;0];
            u_e = [0;0;0;0;0;0];
        end
        
        % State space model
        function [Output] = StateSpaceEqu(obj,x,u,CallFlag)
            % Get input
            ia = u(1);
            ib = u(2);
            ic = u(3);

            % Get state
            vd = x(1);
            vq = x(2);
            theta = x(3);

            % Get parameter
            G = obj.Conductor;
            C = obj.Capacitor;
            w = obj.w;
            
            % State space equations
         	% dx/dt = f(x,u)
            % y     = g(x,u)
            if CallFlag == 1
                % ### Call state equation: dx/dt = f(x,u)
                iabc = [ia;ib;ic];
                idq = SimplusGT.abc2dq(iabc,theta);
                id = idq(1);
                iq = idq(2);
                
                dvd = (id - vd*G + w*C*vq)/C;
                dvq = (iq - vq*G - w*C*vd)/C;
                
                dtheta = w;
                
                f_xu = [dvd; dvq; dtheta]; 
                Output = f_xu;
                
            elseif CallFlag == 2
                % ### Call output equation: y = g(x,u)
                vdq = [vd;vq];
                vabc = SimplusGT.dq2abc(vdq,theta);
                va = vabc(1);
                vb = vabc(2);
                vc = vabc(3);
                g_xu = [va; vb; vc]; 
                Output = g_xu;              
            end
            
        end
        
    end
end