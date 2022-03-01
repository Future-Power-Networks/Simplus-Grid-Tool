% This class defines the model of synchronous machine
% in stationary frame

% Author(s): Yitong Li, Yunjie Gu

%% Notes
%
% The model is in load convention, admittance form.

%% Class

classdef SynchronousMachineStationary < SimplusGT.Class.SynchronousMachine
    
    methods(Static)
        
        function [State,Input,Output] = SignalList(obj)
            State	 = {'i_a','i_b','w','theta'};            
            Input	 = {'v_d','v_q','T_m','v_ex'};
            Output = {'i_d','i_q','w','i_ex','theta'};
        end
        
        function [Output] = StateSpaceEqu(obj,x,u,CallFlag)
            % Get states
            i_a   = x(1);
            i_b   = x(2);
            w     = x(3);
            theta = x(4);

            % Get input signals
            v_d  = u(1);
            v_q  = u(2);
            T_m  = u(3);
            v_ex = u(4);

            % Get parameters
            J  = obj.Para(1);
            D  = obj.Para(2);
            L  = obj.Para(3);
            R  = obj.Para(4);
            W0 = obj.Para(5);

            % ??? Temp
            psi_f = obj.psi_f;
            v_ab = SimplusGT.dq2alphabeta([v_d;v_q],theta);
            i_dq = SimplusGT.alphabeta2dq([i_a;i_b],theta);
            v_a  = v_ab(1);
            v_b  = v_ab(2);
            i_d  = i_dq(1);
            i_q  = i_dq(2);   

            % State space equations
          	% dx/dt = f(x,u)
            % y     = g(x,u)
            if CallFlag == 1        
            % ### Call state equation: dx/dt = f(x,u)
                % Auxiliary equation
                Te = psi_f * i_d;
                
                % State equation
                di_a   = (v_a - R*i_a - psi_f*cos(theta)*w)/L;
                di_b   = (v_b - R*i_b - psi_f*sin(theta)*w)/L;
                dw     = (Te - T_m - D*w)/J;
                dtheta = w;

                f_xu = [di_a; di_b; dw; dtheta];
                Output = f_xu;
            elseif CallFlag == 2    
            % ### Call output equation: y = g(x,u)
                i_ex = 0;           % ??? i_ex = f(v_ex,omega,i_d,i_q) will be added later
                g_xu = [i_d; i_q; w; i_ex; theta];
                Output = g_xu;
            end
        end
        
    end
    
    methods(Access = protected)
        function resetImpl(obj)
            if obj.EquiInitial
                i_dq = obj.x_e(1:2);
                theta = obj.x_e(end);
                i_ab = SimplusGT.dq2alphabeta(i_dq,theta);
                obj.x = [i_ab;obj.x_e(3:end)];
            else
                obj.x = obj.x0;   
            end        
        end
    end

end     % End class definition