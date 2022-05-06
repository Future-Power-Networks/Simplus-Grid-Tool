% This class defines the model of synchronous machine. This model has been
% modified to a grid-forming inverter controlled as a virtual synchronous
% generator for the paper.

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
            State	 = {'i_d','i_q','i_ld','i_lq','v_od','v_oq','w','theta'};            
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
            % e_DQ = V - i_DQ * (R + 1j*L*w);
            e_DQ = V - i_DQ * (R + 1j*L*1.1*w);
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

            i_ld = i_d;
            i_lq = i_q;
            
            v_od = v_d;
            v_oq = v_q;
            
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
            x_e = [i_d; i_q; i_ld; i_lq; v_od; v_oq; w; theta];
            u_e = [v_d; v_q; T_m; v_ex];
            xi  = [xi];
        end
        
        function [Output] = StateSpaceEqu(obj,x,u,CallFlag)
            % Get states
            i_d   = x(1);
            i_q   = x(2);
            i_ld  = x(3);
            i_lq  = x(4);
            v_od  = x(5);
            v_oq  = x(6);
            w     = x(7);
            theta = x(8);

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
            
            % For LCL filter
            Cf = L/10;
            Lc = L/10;
            Rc = R/10;
            
            % =================================
            % Step change of damping for faster convergence of initial
            % states
            % =================================
            if (obj.Timer <= 2) && (obj.Timer >= 1e-3)
                D_new = D*10;
                T_m_new = T_m + D*W0 - D_new*W0;
                
                D = D_new;
                T_m = T_m_new;
            end
            
            % State space equations
          	% dx/dt = f(x,u)
            % y     = g(x,u)
            if CallFlag == 1        
            % ### Call state equation: dx/dt = f(x,u)
                % Auxiliary equation
                psi_f = obj.psi_f;
                if obj.ApparatusType == 0
                    Te = psi_f * i_ld;
                elseif obj.ApparatusType == 1
                    Pe = psi_f*W0*i_ld;      % Means e_d = psi_f*W0 is constant.
                else
                    error(['Error.']);
                end
                
                % State equation
                % Lf equation and w equation
                if obj.ApparatusType == 0
                di_ld   = (v_od - R*i_ld + w*L*i_lq - psi_f*w)/L;
                di_lq   = (v_oq - R*i_lq - w*L*i_ld)/L;
                dw     = (Te - T_m - D*w)/J;
                elseif obj.ApparatusType == 1
             	di_ld   = (v_od - R*i_ld + w*L*i_lq - psi_f*W0)/L;
                di_lq   = (v_oq - R*i_lq - w*L*i_ld)/L;
                dw     = (Pe - T_m*W0 - D*w*W0)/(J*W0);
                end
                
                % Cf equation
                % -(i_ld - i_od) = Cf*dv_cd/dt - w*Cf*v_cq
                % -(i_lq - i_oq) = Cf*dv_cq/dt + w*Cf*v_cd
                dv_od = (-(i_ld - i_d) + w*Cf*v_oq)/Cf;
                dv_oq = (-(i_lq - i_q) - w*Cf*v_od)/Cf;

                % Lc equation
                % v_od - v_d = -(Lc*di_od/dt + Rc*i_od - w*Lc*i_oq)
                % v_oq - v_q = -(Lc*di_oq/dt + Rc*i_oq + w*Lc*i_od)
                di_d = (v_d - v_od - Rc*i_d + w*Lc*i_q)/Lc;
                di_q = (v_q - v_oq - Rc*i_q - w*Lc*i_d)/Lc;
               
                % theta equation
                dtheta = w;
                
                f_xu = [di_d; di_q; di_ld; di_lq; dv_od; dv_oq; dw; dtheta];
                Output = f_xu;
            elseif CallFlag == 2    
            % ### Call output equation: y = g(x,u)
                i_ex = 0;
                
                g_xu = [i_d; i_q; w; i_ex; theta];
                Output = g_xu;
            end
        end
        
    end

end     % End class definition