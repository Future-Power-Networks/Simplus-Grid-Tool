% This class defines the model of a grid-forming wind turbine generator
% system.

% Author(s): Jianhong Zhao

%% Notes
%
% The model is in load convention, admittance form.
%
% The type of generator is permanent magnet synchronous generator (PMSG).

%% Class

classdef WindTurbineGFM < SimplusGT.Class.ModelAdvance
    
    % For temporary use
    properties(Access = protected)
        v_od_r;
        v_oq_r;
        P0;
        Q0;
    end
    
    methods
        % constructor
        function obj = WindTurbineGFM(varargin)
            % Support name-value pair arguments when constructing object
            setProperties(obj,nargin,varargin{:});
        end
    end
    
    methods(Static)
        
        function [State,Input,Output] = SignalList(obj)
         	State  = {'i_ld','i_lq','i_ld_i','i_lq_i','v_od','v_oq','v_od_i','v_oq_i','i_od','i_oq','v_d_ref','w','theta','v_dc','v_dc_i','i_sd','i_sq','i_sd_i','i_sq_i','w_m','w_m_i','theta_m'};
            Input  = {'v_d','v_q','P0'};
            Output = {'i_d','i_q','w','theta','v_dc','i_sd','i_sq','v_sd','v_sq','T_e','T_m','n','theta_e','beta'};
        end
        
        % Calculate the equilibrium
        function [x_e,u_e,xi] = Equilibrium(obj)
         	% Get the power PowerFlow values
            P 	= obj.PowerFlow(1);
            Q	= obj.PowerFlow(2);
            V	= obj.PowerFlow(3);
            xi	= obj.PowerFlow(4);
            w   = obj.PowerFlow(5);
            
            % Get parameters
            xwLf    = obj.Para(1);
            Rf      = obj.Para(2);
            xwCf    = obj.Para(3);
            xwLc    = obj.Para(4);
            Rc      = obj.Para(5);
            Xov     = obj.Para(6);
            Rov     = 0;
            W0      = obj.Para(11);
            V_dc    = obj.Para(12);
            N       = obj.Para(15);

            % PMSG Parameters
            Ld = 0.8e-4;    % D-axis inductance (pu)
            Lq = 1.2e-4;    % Q-axis inductance (pu)
            Rs = 0.01;      % Stator resistance (pu)
            Psi_f = 0.005;  % PM flux linkage (pu)
            Np = 80;        % Number of pole pairs
           
            % Calculate parameters
            % GSC
            Lf = xwLf/W0;
            Cf = xwCf/W0;
            Lc = xwLc/W0;
            v_gd = V;
            v_gq = 0;
            i_od = P/V;
            i_oq = -Q/V;
            v_gdq = v_gd + 1i*v_gq;
            i_odq = i_od + 1i*i_oq;
            v_odq = v_gdq - i_odq*(Rc + 1i*w*Lc);
            i_cdq = v_odq*(1i*w*Cf);
            i_ldq = i_odq - i_cdq;
            i_ld = real(i_ldq);
            i_lq = imag(i_ldq);
            e_dq  = v_odq - i_ldq*(Rf + 1i*w*Lf);
            e_d = real(e_dq);
            e_q = imag(e_dq);
            i_ld_i = e_d;
            i_lq_i = e_q;
            v_od = real(v_odq);
            v_oq = imag(v_odq);
            v_od_i = -i_ld;
            v_oq_i = -i_lq;
            i_od = real(i_odq);
            i_oq = imag(i_odq);
            P_dc = e_d*i_ld + e_q*i_lq;
            theta = xi;
            obj.P0 = P*(-1);
            obj.Q0 = Q*(-1);
            v_odq_r = v_odq + (Rov + 1i*Xov)*i_odq*(-1);
            v_od_r = real(v_odq_r);
            v_oq_r = imag(v_odq_r);
            obj.v_od_r = v_od_r;
            obj.v_oq_r = v_oq_r;
            % MSC
            w_m = N*pi/30;
            w_e = Np*w_m;
            theta_m = 0;
            e_md = 0;
            e_mq = sqrt(3/2)*w_e*Psi_f;
            i_sd = 0;
            % Calculate i_sq according to:
            % P_dc = (i_sd^2 + i_sq^2)*Rs + e_md*i_sd + e_mq*i_sq
            if Rs == 0
                i_sq = P_dc/e_mq;
            else
                i_sq = (sqrt(e_mq^2 + 4*Rs*P_dc) - e_mq)/(2*Rs);
            end
            v_sd = Rs*i_sd - w_e*Lq*i_sq + e_md;
            v_sq = Rs*i_sq + w_e*Ld*i_sd + e_mq;
            i_sd_i = v_sd;
            i_sq_i = v_sq;
            v_dc_i = i_sq;
            v_dc = V_dc;
            beta = (1-obj.P0)*20;    % Approximation
            if beta < 0
                beta = 0;
            end
            w_m_i = beta;

            % Get equilibrium
            x_e = [i_ld; i_lq; i_ld_i; i_lq_i; v_od; v_oq; v_od_i; v_oq_i; i_od; i_oq; v_od_r; w; theta; v_dc; v_dc_i; i_sd; i_sq; i_sd_i; i_sq_i; w_m; w_m_i; theta_m];
            u_e = [v_gd; v_gq; 0];
        end
        
        % State space model
        function [Output] = StateSpaceEqu(obj,x,u,CallFlag)
            % Get inputs
            v_gd   = u(1);
            v_gq   = u(2);

            % Get states
            i_ld    = x(1);
            i_lq    = x(2);
            i_ld_i  = x(3);
            i_lq_i  = x(4);
            v_od    = x(5);
            v_oq    = x(6); 
            v_od_i  = x(7);
            v_oq_i  = x(8);
            i_od    = x(9);
            i_oq    = x(10);
            v_od_r  = x(11);
            w       = x(12);
            theta   = x(13);
            v_dc    = x(14);
            v_dc_i  = x(15);
            i_sd    = x(16);
            i_sq    = x(17);
            i_sd_i  = x(18);
            i_sq_i  = x(19);
            w_m     = x(20);
            w_m_i   = x(21);
            theta_m = x(22);
            
            % Get parameters
            xwLf    = obj.Para(1);
            Rf      = obj.Para(2);
            xwCf    = obj.Para(3);
            xwLc    = obj.Para(4);
            Rc      = obj.Para(5);
            Xov     = obj.Para(6);
            Rov     = 0;
            xDw     = obj.Para(7);
            Dv      = 0.05;
            xfdroop = obj.Para(8);
            xfvdq   = obj.Para(9);
            xfidq   = obj.Para(10);
            W0      = obj.Para(11);
            v_dc_r  = obj.Para(12);
            C_dc    = obj.Para(13);
            v_wpu   = obj.Para(14);
            n_r     = obj.Para(15);
            f_v_dc  = obj.Para(16);
            f_i_sdq = obj.Para(17);
            f_w_m   = obj.Para(18);

            % Wind Turbine Parameters
            r = 50;             % Wind turbine radius (m)
            beta_limit_H = 25;  % Pitch angle limit high (deg)
            beta_limit_L = 0;   % Pitch angle limit low (deg)
            V_w0 = 12;          % Base Wind speed (m/s)
            Cp_max = 0.48;      % Maximum wind energy utilization coefficient
            % Notes:
            % In grid-forming mode, the pitch angle (beta) is not fixed
            % so that the output power of wind turbine can be adjusted.

            % PMSG parameters
            Ld = 0.8e-4;        % D-axis inductance (pu)
            Lq = 1.2e-4;        % Q-axis inductance (pu)
            Rs = 0.01;          % Stator resistance (pu)
            Psi_f = 0.005;      % PM flux linkage (pu)
            J = 1;              % Rotor inertia (pu)
            D = 0.0001;         % Viscous damping coefficient (pu)
            Np = 80;            % Number of pole pairs
            % Notes:
            % Psi_fpu = Psi_f/Vbase, Jpu = J/Pbase, Dpu = D/Pbase,
            % because w_m, w_e and n are not in per unit system in this
            % model.
            
            % Filter parameter
            Lf = xwLf/W0;
            Cf = xwCf/W0;
            Lc = xwLc/W0;

            % Droop controller parameter
            Dw = xDw*W0;
            wf = xfdroop*2*pi;

            % Current controller parameter
            w_i_ldq = xfidq*2*pi;
            kp_i_ldq = w_i_ldq*Lf;
            ki_i_ldq = w_i_ldq^2*Lf/4;

            % Voltage controller parameter
            w_v_odq = xfvdq*2*pi;
            kp_v_odq = w_v_odq*Cf;
            ki_v_odq = w_v_odq^2*Cf/4*50;

            % Dc link controller parameter
            w_v_dc = f_v_dc*2*pi;
            kp_v_dc	= v_dc_r*C_dc*w_v_dc;
            ki_v_dc	= kp_v_dc*w_v_dc/4;

            % PMSG current controller parameter
            w_i_sdq = f_i_sdq*2*pi;
            kp_i_sdq = w_i_sdq*(Ld+Lq)/2;
            ki_i_sdq = kp_i_sdq*w_i_sdq/4;

            % PMSG speed controller parameter
            w_w_m = f_w_m*2*pi;
            kp_w_m = J*w_w_m*20;
            ki_w_m = kp_w_m*w_w_m/4;

            % Angular speed reference
            w_m_r = n_r*pi/30;

            % Wind speed
            v_w = v_wpu*V_w0;
            

            % Power measurement
            p = (v_od*i_od + v_oq*i_oq)*(-1);
            q = (-v_od*i_oq + v_oq*i_od)*(-1);

            % Droop controller
            P0 = obj.P0;
            Q0 = obj.Q0;
            dw = (W0 + Dw*(P0 - p) - w)*wf;         % P-w droop

            % Phase angle
            dtheta = w;
            
            switch 1
                case 1
                    dv_od_r = (obj.v_od_r + Dv*(Q0 - q) - v_od_r)*wf;   % Q-V droop
                case 2
                    v_od_r = obj.v_od_r + Dv*(Q0 - q);                  % Q-V droop without LPF
                    dv_od_r = 0;
                case 3
                    v_od_r = obj.v_od_r;                                % No Q-V droop
                    dv_od_r = 0;
                otherwise
                    error('Error')
            end
            v_oq_r = obj.v_oq_r;

            % AC voltage control
            error_v_od = v_od_r - v_od - (i_od*Rov-i_oq*Xov)*(-1);
          	error_v_oq = v_oq_r - v_oq - (i_oq*Rov+i_od*Xov)*(-1);
            i_ld_r = -(error_v_od*kp_v_odq + v_od_i);
            i_lq_r = -(error_v_oq*kp_v_odq + v_oq_i);
            dv_od_i = error_v_od*ki_v_odq;
            dv_oq_i = error_v_oq*ki_v_odq;

            % AC current control
            error_i_ld = i_ld_r-i_ld;
            error_i_lq = i_lq_r-i_lq;
            e_d = (-error_i_ld*kp_i_ldq + i_ld_i)/v_dc_r*v_dc;
            e_q = (-error_i_lq*kp_i_ldq + i_lq_i)/v_dc_r*v_dc;
            di_ld_i = -error_i_ld*ki_i_ldq;
            di_lq_i = -error_i_lq*ki_i_ldq;

            % Dc link control
            i_sd_r = 0;
            i_sq_r = -(v_dc_r - v_dc)*kp_v_dc + v_dc_i;
            dv_dc_i = -(v_dc_r - v_dc)*ki_v_dc;
            
            % PMSG current control
            v_sd = ((i_sd_r - i_sd)*kp_i_sdq + i_sd_i)/v_dc_r*v_dc;
            v_sq = ((i_sq_r - i_sq)*kp_i_sdq + i_sq_i)/v_dc_r*v_dc;
            di_sd_i = (i_sd_r - i_sd)*ki_i_sdq;
            di_sq_i = (i_sq_r - i_sq)*ki_i_sdq;

            % PMSG speed control
            beta = -(w_m_r - w_m)*kp_w_m + w_m_i;
            dw_m_i = -(w_m_r - w_m)*ki_w_m;
            % Anti wind-up for speed control
            if (w_m_i >= beta_limit_H && dw_m_i >= 0) || (w_m_i <= beta_limit_L && dw_m_i <= 0)
                dw_m_i = 0;
            end
            % Pitch angle saturation
            beta = min(beta,beta_limit_H);
            beta = max(beta,beta_limit_L);

            % Wind turbine equations
            lambda = w_m*r/v_w;     % Tip-speed ratio
            lambda_i = 1/(1/(lambda + 0.08*beta) - 0.035/(beta^3 + 1));
            Cp = 0.5176*(116/lambda_i - 0.4*beta - 5)*exp(-21/lambda_i) + 0.0068*lambda; % Wind energy utilization coefficient
            P_m = v_wpu^3*(Cp/Cp_max);  % Originally, P_m = (1/2)*rho*A*v_w^3*Cp
            T_m = -P_m/w_m;

            % PMSG mechanical dynamics
            T_e = Np*(sqrt(3/2)*Psi_f*i_sq + (Ld - Lq)*i_sd*i_sq);
            dw_m = (T_e - T_m - D*w_m)/J;
            w_e = Np*w_m;
            n = 30*w_m/pi;
            dtheta_m = w_m;
            theta_e = Np*theta_m;

            % PMSG stator dynamics
            di_sd = (v_sd - Rs*i_sd + w_e*Lq*i_sq)/Ld;
            di_sq = (v_sq - Rs*i_sq - w_e*Ld*i_sd - sqrt(3/2)*w_e*Psi_f)/Lq;
            
            % Dc capacitor
            dv_dc = (e_d*i_ld + e_q*i_lq - v_sd*i_sd - v_sq*i_sq)/v_dc/C_dc;

            % Lf
            di_ld = (v_od - e_d - Rf*i_ld + w*Lf*i_lq)/Lf;
            di_lq = (v_oq - e_q - Rf*i_lq - w*Lf*i_ld)/Lf;

            % Cf
            dv_od = (-(i_ld - i_od) + w*Cf*v_oq)/Cf;
            dv_oq = (-(i_lq - i_oq) - w*Cf*v_od)/Cf;

            % Lc
            di_od = (v_gd - v_od - Rc*i_od + w*Lc*i_oq)/Lc;
            di_oq = (v_gq - v_oq - Rc*i_oq - w*Lc*i_od)/Lc;

            % State space equations
         	% dx/dt = f(x,u)
            % y     = g(x,u)
            if CallFlag == 1    
            % ### Call state equation: dx/dt = f(x,u)
                f_xu = [di_ld; di_lq; di_ld_i; di_lq_i; dv_od; dv_oq; dv_od_i; dv_oq_i; di_od; di_oq; dv_od_r; dw; dtheta; dv_dc; dv_dc_i; di_sd; di_sq; di_sd_i; di_sq_i; dw_m; dw_m_i; dtheta_m];
                Output = f_xu;
            elseif CallFlag == 2     
            % ### Call output equations: y = g(x,u)
                g_xu = [i_od; i_oq; w; theta; v_dc; i_sd; i_sq; v_sd; v_sq; T_e; T_m; n; theta_e; beta];
                Output = g_xu;
            end
        end
    end
end     % End class definition
