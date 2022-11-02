% Author(s): Yitong Li

%% Class

classdef ThreePhaseSynchronousSeriesRL < SimplusGT.Class.ModelAdvanceNetwork
    
    properties
        R;          % Resistor
        L;          % Inductor
        w;          % Frequency
    end
    
    methods
        % Constructor
        function obj = ThreePhaseSynchronousSeriesRL(varargin)
            setProperties(obj,nargin,varargin{:});
        end
    end
    
    methods(Static)
        % Set the strings of input, output, state
        % For an inductor, voltage v is input, current i is output, current
        % i is state.
        function [State,Input,Output] = SignalList(obj)
            State  = {'id','iq','theta'};                         	% x
            Input  = {'va1','vb1','vc1','va2','vb2','vc2'};         % u
            Output = {'ia','ib','ic'};                              % y
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
            va1 = u(1);
            vb1 = u(2);
            vc1 = u(3);
            va2 = u(4);
            vb2 = u(5);
            vc2 = u(6);
            
            % Get state
            id = x(1);
            iq = x(2);
            theta = x(3);

            % Get parameter
            L = obj.L;
            R = obj.R;
            w = obj.w;
            
            % State space equations
         	% dx/dt = f(x,u)
            % y     = g(x,u)
            if CallFlag == 1
                % ### Call state equation: dx/dt = f(x,u)
                va = va1-va2;
                vb = vb1-vb2;
                vc = vc1-vc2;
                vabc = [va;vb;vc];
                vdq = SimplusGT.abc2dq(vabc,theta);
                vd = vdq(1);
                vq = vdq(2);
                
                did = (vd - id*R + w*L*iq)/L;
                diq = (vq - iq*R - w*L*id)/L;
                
                dtheta = w;
                
                f_xu = [did; diq; dtheta]; 
                Output = f_xu;
                
            elseif CallFlag == 2
                % ### Call output equation: y = g(x,u)
                idq = [id;iq];
                iabc = SimplusGT.dq2abc(idq,theta);
                ia = iabc(1);
                ib = iabc(2);
                ic = iabc(3);
                g_xu = [ia; ib; ic]; 
                Output = g_xu;              
            end
            
        end
        
    end
end