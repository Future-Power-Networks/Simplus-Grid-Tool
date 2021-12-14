% This class defines the model of synchronous machine

% Author(s): Yitong Li, Yunjie Gu

%% Notes
%
% The model is in load convention, admittance form.
%
% Two very important physical relationships:
% w * psi = v
% w * T = P
% The per unit selection of w will influence the results a lot.

%% Class

classdef SynchronousMachine < SimplusGT.Class.ModelAdvance
    
    properties(Access = protected)
        psi_f;
    end
    
    methods
        % constructor
        function obj = SynchronousMachine(varargin)
            % Support name-value pair arguments when constructing object
            setProperties(obj,nargin,varargin{:});
        end
    end
    
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
        	D  = obj.Para(2);
            wL = obj.Para(3);
            R  = obj.Para(4);
            W0 = obj.Para(5);
            
            % Calculate paramters
            D = D/W0^2;
            L = wL/W0;
            
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

            % Notes:
            % Noting that w is Wbase rather than 1 at steady state in this
            % model, i.e., w is not in per unit. This means psi_f, T_m, T_e
            % (or K_S) is also NOT close to 1. This is different from the
            % analysis in Kundur's book, where w is also in per unit
            % value.
            
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
            wL = obj.Para(3);
            R  = obj.Para(4);
            W0 = obj.Para(5);
            
            % Calculate parameters
            J = J*2/W0^2;   % Jpu=J/Pb=[1/2*J*w0^2/Pb]*2/w0^2, [MWs/MW] 
            D = D/W0^2;     % Dpu=dTpu/dw=dPpu/dw/w0=[dP%/dw%]/w0^2, [%/%]
            L = wL/W0;
            
            % Notes:
            % Noting that w is not in per unit system here. So, J and D
            % should be divided by W0^2 rather than W0. In this case,
            % P=T*w0. If P is in per unit and is close to 1, T is close to
            % 1/w0 and is much smaller than 1. For the two forms of swing
            % equations blow:
            % J*dw/dt = Tm - Ks - Kd*w;    (1)
            % J*dw/dt = Pm - Ks - Kd*w;    (2)
            % (1)*w0 is equivalent to (2), which means J, Ks, Kd in (1) are
            % w0 times smaller.
            %
            % This is different from the equations in Kundur's book, where
            % w is also in per unit and P=T.
            
            % State space equations
          	% dx/dt = f(x,u)
            % y     = g(x,u)
            if CallFlag == 1        
            % ### Call state equation: dx/dt = f(x,u)
                % Auxiliary equation
                psi_f = obj.psi_f;
                psi_d = L*i_d;
                psi_q = L*i_q - psi_f;
                if obj.ApparatusType == 0
                    Te = psi_f * i_d;
                elseif obj.ApparatusType == 1
                    Pe = psi_f*W0*i_d;      % Means e_d = psi_f*W0 is constant.
                else
                    error(['Error.']);
                end
                
                % State equation
                if obj.ApparatusType == 0
                di_d   = (v_d - R*i_d + w*psi_q)/L;
                di_q   = (v_q - R*i_q - w*psi_d)/L;
                dw     = (Te - T_m - D*w)/J;
                elseif obj.ApparatusType == 1
             	di_d   = (v_d - R*i_d + w*L*i_q - psi_f*W0)/L;
                di_q   = (v_q - R*i_q - w*L*i_d)/L;
                dw     = (Pe - T_m*W0 - D*w*W0)/(J*W0);
                end
                
             	% Notes:
                %
                % For type 1, we can move the flux inductor outside the SG.
                % This operation is based on a precondition that the Pe is
                % measured from the internal EMF rather than the output
                % electrical terminal.
                % 
                % q-axis is aligned to psi_f, rather than d-axis. q-axis is
                % leading d-axis 90 degrees. These two conventions would be
                % different from conventional SG models for example in
                % Kundur's book.
                
                dtheta = w;
                
                % Notes:
                % Type 0 has more dampings than type 1.

                f_xu = [di_d; di_q; dw; dtheta];
                Output = f_xu;
            elseif CallFlag == 2    
            % ### Call output equation: y = g(x,u)
                i_ex = 0;           % ??? i_ex = f(v_ex,omega,i_d,i_q) will be added later
                
                g_xu = [i_d; i_q; w; i_ex; theta];
                Output = g_xu;
            end
        end
        
    end

end     % End class definition