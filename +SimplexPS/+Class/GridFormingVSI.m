% This class defines the model of grid-forming voltage source inverters.

% Author(s): Yitong Li

%% Notes
%
% The model is in load convention, admittance form.

%% Class

classdef GridFormingVSI < SimplexPS.Class.ModelAdvance
    
    methods
        % constructor
        function obj = GridFormingVSI(varargin)

            % Support name-value pair arguments when constructing object
            setProperties(obj,nargin,varargin{:});

        end
    end
    
    methods(Static)
        
        function [State,Input,Output] = SignalList(obj)
         	State  = {'i_ld','i_lq','i_ld_i','i_lq_i','v_od','v_oq','v_od_i','v_oq_i','i_d','i_q','w','v_od_r'};
            Input  = {'v_d','v_q','?','?'};
            Output = {'i_d','i_q','w','?'};
        end
        
        % Calculate the equilibrium
        function Equilibrium(obj)
         	% Get the power PowerFlow values
            P 	= obj.PowerFlow(1);
            Q	= obj.PowerFlow(2);
            V	= obj.PowerFlow(3);
            xi	= obj.PowerFlow(4);
            w   = obj.PowerFlow(5);
            
            % Input
            v_d = 1;
            v_q = 0;
            
            % State
            i_ld = 0;
            i_lq = 0;
            i_ld_i = 0;
            i_lq_i = 0;
            v_od = 1;
            v_oq = 0;
            v_od_i = 0;
            v_oq_i = 0;
            i_d = 0;
            i_q = 0;
            w_e = w;
            v_od_r = 1;
            
            % Get equilibrium
            obj.x_e = [i_ld;i_lq;i_ld_i;i_lq_i;v_od;v_oq;v_od_i;v_oq_i;i_d;i_q;w_e;v_od_r];
            obj.u_e = [v_d;v_q;0;0];
            obj.xi = 0;
        end
        
        % State space model
        function [Output] = StateSpaceEqu(obj,x,u,CallFlag)
            % Get input
            v_d   = u(1);
            v_q   = u(2);
            u(3);
            u(4);

            % Get state
            i_ld   = x(1);
            i_lq   = x(2);
            i_ld_i = x(3);
            i_lq_i = x(4);
            v_od   = x(5);
            v_oq   = x(6); 
            v_od_i = x(6);
            v_oq_i = x(7);
            i_d    = x(8);
            i_q    = x(9);
            w      = x(10);
            v_od_r = x(11);

            % Get parameter
            % LCL filter
            Lf       = obj.Para(1);
            Rf       = obj.Para(2);
            Cf       = obj.Para(3);
            Lc       = obj.Para(4);
            Rc       = obj.Para(5);
            % Droop control
            Dw       = obj.Para(6);     % PF droop gain
            Dv       = obj.Para(7);     % QV droop gain
            Tf       = obj.Para(8);
         	P0       = obj.Para(9);
            Q0       = obj.Para(10);
            w0       = obj.Para(11);
            v_od0    = obj.Para(12);
            v_oq0    = obj.Para(13);
            % Voltage and current loops
            kp_i_ldq = obj.Para(14);
            ki_i_ldq = obj.Para(15);
            kp_v_odq = obj.Para(16);
            ki_v_odq = obj.Para(17);
            Gi_cd    = obj.Para(18);     % Cross-decoupling gain
            Gv_cd    = obj.Para(19);     % Cross-decoupling gain
            Fv       = obj.Para(20);     % Voltage feedforward gain
            Fi       = obj.Para(21);     % Current feedforward gain
            % Virtual impedance control
            Rov      = obj.Para(22);
            Xov      = obj.Para(23);

            % State space equations
            if CallFlag == 1
                % State equations: dx/dt = f(x,u)
                % Power measurement
                p = (v_od*i_d + v_oq*i_q)*(-1);
                q = (-v_od*i_q + v_oq*i_d)*(-1);

                % Droop controller
                % w = w0 + Dw*(P0 - P*LPF)
                % v_od_r = v_od_0 + Dv*(Q0 - Q*LPF)
                % v_oq_r = v_oq_0
                dw = (w0 + Dw*(P0 - p) - w)/Tf;
                dv_od_r = (v_od0 + Dv*(Q0 - p) - v_od_r)/Tf;
                v_oq_r = v_oq0;

                % Voltage control
                delta_v_od = v_od_r - v_od - (i_d*Rov-i_q*Xov)*(-1);
                delta_v_oq = v_oq_r - v_oq - (i_q*Rov+i_d*Xov)*(-1);
                i_ld_r = delta_v_od*kp_v_odq + v_od_i - Gv_cd*w0*Cf*v_oq + Fi*i_ld;
                i_lq_r = delta_v_oq*kp_v_odq + v_oq_i + Gv_cd*w0*Cf*v_od + Fi*i_lq;
                dv_od_i = delta_v_od*ki_v_odq;
                dv_oq_i = delta_v_oq*ki_v_odq;

                % Current control
                e_d = -(i_ld_r-i_ld)*kp_i_ldq + i_ld_i - Gi_cd*w0*Lf*(-i_lq) + Fv*v_od;
                e_q = -(i_lq_r-i_lq)*kp_i_ldq + i_lq_i + Gi_cd*w0*Lf*(-i_ld) + Fv*v_oq;
                di_ld_i = -(i_ld_r - i_ld)*ki_i_ldq;
                di_lq_i = -(i_lq_r - i_lq)*ki_i_ldq;

                % Lf equation
                % e_d - v_od = -(di_ld/dt*Lf + Rf*i_ld - w*Lf*i_lq)
                % e_q - v_oq = -(di_lq/dt*Lf + Rf*i_lq + w*Lf*i_ld)
                di_ld = (v_od - e_d - Rf*i_ld + w*Lf*i_lq)/Lf;
                di_lq = (v_oq - e_q - Rf*i_lq - w*Lf*i_ld)/Lf;

                % Cf equation
                % -(i_od - i_ld) = Cf*dv_od/dt - w*Cf*v_oq
                % -(i_oq - i_lq) = Cf*dv_oq/dt + w*Cf*v_od
                dv_od = -(i_d - i_ld) + w*Cf*v_oq;
                dv_oq = -(i_q - i_lq) - w*Cf*v_od;

                % Lc equation
                % v_od - v_d = -(Lc*di_od/dt + Rc*i_od - w*Lc*i_oq)
                % v_oq - v_q = -(Lc*di_oq/dt + Rc*i_oq + w*Lc*i_od)
                di_d = (v_d - v_od - Rc*i_d + w*Lc*i_q)/Lc;
                di_q = (v_q - v_oq - Rc*i_q - w*Lc*i_d)/Lc;

                f_xu = [di_ld;di_lq;di_ld_i;di_lq_i;dv_od;dv_oq;dv_od_i;dv_oq_i;di_d;di_q;dw;dv_od_r];
                Output = f_xu;
            elseif CallFlag == 2
                % Output equations: y = g(x,u)
                g_xu = [i_d;i_q;w;0];
                Output = g_xu;
            end
            
        end
        
    end
end