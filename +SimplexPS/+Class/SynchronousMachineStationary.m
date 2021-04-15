% This class defines the model of synchronous machine

% Author(s): Yitong Li, Yunjie Gu

%% Notes
%
% The model is in load convention, admittance form.

%% Class

classdef SynchronousMachineStationary < SimplexPS.Class.ModelAdvanceStationary
    
    properties(Access = protected)
        psi_f;
    end
    
%     methods
%         % constructor
%         function obj = SynchronousMachine(varargin)
% 
%             % Support name-value pair arguments when constructing object
%             setProperties(obj,nargin,varargin{:});
% 
%         end
%     end
    
    methods(Static)
        
        function [State,Input,Output] = SignalList(obj)
            State	 = {'i_d','i_q','w','theta'};            
            Input	 = {'v_d','v_q','T_m','v_ex'};
            Output = {'i_d','i_q','w','i_ex','theta'};
        end
        
        function [x_e,u_e,xi] = Equilibrium(obj)
            % Get the power PowerFlow values
            P 	= obj.PowerFlow(1);
            Q	= obj.PowerFlow(2);
            V	= obj.PowerFlow(3);
            xi	= obj.PowerFlow(4);
            w   = obj.PowerFlow(5);
            
            % Get parameters
        	D = obj.Para(2);
            L = obj.Para(3);
            R = obj.Para(4);

            % Calculate parameters
            i_D = P/V;
            i_Q = -Q/V;     % Use -Q because S = V*conj(I)
            i_DQ = i_D + 1j*i_Q;
            e_DQ = V - i_DQ * (R + 1j*L*w);
            arg_e = angle(e_DQ);
            abs_e = abs(e_DQ);
            xi = xi + arg_e;

            v_dq = V * exp(-1j*arg_e);
            i_dq = i_DQ * exp(-1j*arg_e);
            v_d = real(v_dq);
            v_q = imag(v_dq);
            i_d = real(i_dq);
            i_q = imag(i_dq);
            psi_f = abs_e/w;
            T_m = psi_f * i_d - D*w;

            % ??? Temp
            obj.psi_f = psi_f;
            v_ex = 0;
            theta = xi;

            % Get equilibrium
            x_e = [i_d; i_q; w; theta];
            u_e = [v_d; v_q; T_m; v_ex];
            xi  = [xi];
        end
        
        function [Output] = StateSpaceEqu(obj,x,u,CallFlag)
            % Get states
            i_d   = x(1);
            i_q   = x(2);
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

            % State space equations
            if CallFlag == 1        % State equation
                % Auxiliary equation
                psi_d = L*i_d;
                psi_q = L*i_q - psi_f;
                Te = psi_f * i_d;
                
                % State equation
                di_d   = (v_d - R*i_d + w * psi_q)/L;
                di_q   = (v_q - R*i_q - w * psi_d)/L;
                dw     = (Te - T_m - D*w)/J;
                dtheta = w;

                f_xu = [di_d; di_q; dw; dtheta];
                Output = f_xu;
            elseif CallFlag == 2    % Output equation
                i_ex = 0;           % ??? i_ex = f(v_ex,omega,i_d,i_q) will be added later
                
%              	  e_d = w*psi_f;
%                 e_q = 0;
%                 e_dq = e_d+1i*e_q;
%                 v_dq = e_dq + (i_d+1i*i_q)*(R+1i*w*L);
%                 theta = angle(v_dq)+theta;
                
                g_xu = [i_d; i_q; w; i_ex; theta];
                Output = g_xu;
            end
        end
        
    end

end     % End class definition