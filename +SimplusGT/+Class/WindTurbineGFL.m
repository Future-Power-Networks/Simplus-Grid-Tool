% This class defines the model of a grid-following wind turbine generator
% system.

% Author(s): Jianhong Zhao

%% Notes
%
% The model is in load convention, admittance form.
%
% The type of generator is permanent magnet synchronous generator (PMSG).

%% Class

classdef WindTurbineGFL < SimplusGT.Class.ModelAdvance
    
    % For temporary use
    properties(Access = protected)
        i_q_r;
    end 
    
    methods
        % constructor
        function obj = WindTurbineGFL(varargin)
            % Support name-value pair arguments when constructing object
            setProperties(obj,nargin,varargin{:});
        end
    end

    methods(Static)
        
        function [State,Input,Output] = SignalList(obj)
            State = {'i_d','i_q','i_d_i','i_q_i','w_pll_i','w','theta','v_dc','v_dc_i','i_sd','i_sq','i_sd_i','i_sq_i','w_m','w_m_i','theta_m'};
        	Input = {'v_d','v_q','ang_r'};
            Output = {'i_d','i_q','w','theta','v_dc','i_sd','i_sq','v_sd','v_sq','T_e','T_m','n','theta_e'};
        end
        
        function [x_e,u_e,xi] = Equilibrium(obj)
            % Get the power PowerFlow values
            P 	= obj.PowerFlow(1);
            Q	= obj.PowerFlow(2);
            V	= obj.PowerFlow(3);
            xi	= obj.PowerFlow(4);
            w   = obj.PowerFlow(5);

            % Get parameters
            V_dc = obj.Para(1);
            wLf  = obj.Para(3);
            Rf   = obj.Para(4);
            W0   = obj.Para(9);
            N    = obj.Para(11);
            Lf   = wLf/W0;

            % PMSG Parameters
            Ld = 0.8e-4;    % D-axis inductance (pu)
            Lq = 1.2e-4;    % Q-axis inductance (pu)
            Rs = 0.01;      % Stator resistance (pu)
            Psi_f = 0.005;  % PM flux linkage (pu)
            Np = 80;        % Number of pole pairs

            % Calculate parameters
            % GSC
            i_d = P/V;
            i_q = -Q/V;
            v_d = V;
            v_q = 0;
            i_dq = i_d + 1j*i_q;
            v_dq = v_d + 1j*v_q;
            e_dq = v_dq - i_dq * (Rf + 1j*Lf*w);
            e_d = real(e_dq);
            e_q = imag(e_dq);
            i_d_i = e_d;
            i_q_i = e_q;
            i_q_r = i_q;
            w_pll_i = w;
            v_dc_i = i_d;
            v_dc = V_dc;
            P_dc = e_d*i_d + e_q*i_q;
            ang_r = 0;
            theta = xi;
            obj.i_q_r = i_q_r;
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
            w_m_i = i_sq;

            % Get equilibrium
            x_e = [i_d; i_q; i_d_i; i_q_i; w_pll_i; w; theta; v_dc; v_dc_i; i_sd; i_sq; i_sd_i; i_sq_i; w_m; w_m_i; theta_m];
        	u_e = [v_d; v_q; ang_r];
        end

        function [Output] = StateSpaceEqu(obj,x,u,CallFlag)
           	% Get parameters
            v_dc_r      = obj.Para(1);
            C_dc        = obj.Para(2);
            wLf         = obj.Para(3);
            Rf          = obj.Para(4);
            f_v_dc      = obj.Para(5);
            f_pll       = obj.Para(6);     
            f_tau_pll   = obj.Para(7);
            f_i_dq      = obj.Para(8); 
            W0          = obj.Para(9);
            v_wpu       = obj.Para(10);
            n_r         = obj.Para(11);
            f_i_sdq     = obj.Para(12);
            f_w_m       = obj.Para(13);

            % Wind Turbine Parameters
            r = 50;         % Radius (m)
            beta = 0;       % Pitch angle (deg)
            V_w0 = 12;      % Base wind speed (m/s)
            Cp_max = 0.48;  % Maximum wind energy utilization coefficient
            % Notes:
            % In grid-following mode, the pitch angle (beta) is fixed at
            % zero to capture as much wind energy as possible.

            % PMSG Parameters
            Ld = 0.8e-4;    % D-axis inductance (pu)
            Lq = 1.2e-4;    % Q-axis inductance (pu)
            Rs = 0.01;      % Stator resistance (pu)
            Psi_f = 0.005;  % PM flux linkage (pu)
            J = 1;          % Rotor inertia (pu)
            D = 0.0001;     % Viscous damping coefficient (pu)
            Np = 80;        % Number of pole pairs
            % Notes:
            % Psi_fpu = Psi_f/Vbase, Jpu = J/Pbase, Dpu = D/Pbase,
            % because w_m, w_e and n are not in per unit system in this
            % model.
            
            % Filter inductor
            Lf = wLf/W0;
            
            % Dc link controller parameter
            w_v_dc  = f_v_dc*2*pi;
            kp_v_dc	= v_dc_r*C_dc*w_v_dc;
            ki_v_dc	= kp_v_dc*w_v_dc/4;
            
            % PLL controller parameter
            w_pll     = f_pll*2*pi;
            kp_pll    = w_pll;
            ki_pll    = kp_pll * w_pll/4;
            w_tau_pll = f_tau_pll*2*pi;
            tau_pll   = 1/w_tau_pll;
            
            % Current controller parameter
            w_i_dq  = f_i_dq*2*pi;
            kp_i_dq = Lf * w_i_dq;
            ki_i_dq = Lf * w_i_dq^2 /4;

            % PMSG speed controller parameter
            w_w_m  = f_w_m*2*pi;
            kp_w_m = w_w_m*J/(sqrt(3/2)*Np*Psi_f);
            ki_w_m = kp_w_m*w_w_m/4;

            % PMSG current controller parameter
            w_i_sdq  = f_i_sdq*2*pi;
            kp_i_sdq = w_i_sdq*(Ld+Lq)/2;
            ki_i_sdq = kp_i_sdq*w_i_sdq/4;

            % Angular speed reference
            w_m_r = n_r*pi/30;

            % Wind speed
            v_w = v_wpu*V_w0;
            
            % Get states
      	    i_d   	= x(1);
     	    i_q   	= x(2);
      	    i_d_i  	= x(3);
            i_q_i 	= x(4);
            w_pll_i = x(5);
            w       = x(6);
            theta   = x(7);
            v_dc  	= x(8);
            v_dc_i 	= x(9);
            i_sd    = x(10);
            i_sq    = x(11);
            i_sd_i  = x(12);
            i_sq_i  = x(13);
            w_m     = x(14);
            w_m_i   = x(15);
            theta_m = x(16);

            % Get inputs
        	v_d    = u(1);
            v_q    = u(2);
            ang_r  = u(3);
            
                      
            % Dc link control
            i_d_r = (v_dc_r - v_dc)*kp_v_dc + v_dc_i;
            i_q_r = obj.i_q_r;    % Constant iq control, PQ/PV node in power flow
          	dv_dc_i = (v_dc_r - v_dc)*ki_v_dc;

            % PLL angle measurement
            e_ang = v_q - ang_r;    % vq-PLL
            % Notes:
            % "- ang_r" gives the reference in load convention, like
            % the Tw port.

            % PLL Integral controller 
            dw_pll_i = e_ang*ki_pll;
            
            if 1                                                                          
                dw = (w_pll_i + e_ang*kp_pll - w)/tau_pll;  	% LPF
                % Notes:
                % This introduces an additional state w.
            else
                dw = 0;                                         % No LPF
                w = w_pll_i + e_ang*kp_pll;
            end
            
            % Phase angle
            dtheta = w;
                        
            % Ac current control
            e_d = (-(i_d_r - i_d)*kp_i_dq + i_d_i)/v_dc_r*v_dc;
            e_q = (-(i_q_r - i_q)*kp_i_dq + i_q_i)/v_dc_r*v_dc;
            di_d_i = -(i_d_r - i_d)*ki_i_dq;
            di_q_i = -(i_q_r - i_q)*ki_i_dq;

            % PMSG speed control
            i_sd_r = 0;
            i_sq_r = (w_m_r - w_m)*kp_w_m + w_m_i;
            dw_m_i = (w_m_r - w_m)*ki_w_m;

            % PMSG current control
            v_sd = ((i_sd_r - i_sd)*kp_i_sdq + i_sd_i)/v_dc_r*v_dc;
            v_sq = ((i_sq_r - i_sq)*kp_i_sdq + i_sq_i)/v_dc_r*v_dc;
            di_sd_i = (i_sd_r - i_sd)*ki_i_sdq;
            di_sq_i = (i_sq_r - i_sq)*ki_i_sdq;

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
            dv_dc = (e_d*i_d + e_q*i_q - v_sd*i_sd - v_sq*i_sq)/v_dc/C_dc;
            
            % Ac filter inductor
          	di_d = (v_d - Rf*i_d + w*Lf*i_q - e_d)/Lf;
            di_q = (v_q - Rf*i_q - w*Lf*i_d - e_q)/Lf;
            
            % State space equations
            % dx/dt = f(x,u)
            % y     = g(x,u)
            if CallFlag == 1    
            % ### Call state equation: dx/dt = f(x,u)
                f_xu = [di_d; di_q; di_d_i; di_q_i; dw_pll_i; dw; dtheta; dv_dc; dv_dc_i; di_sd; di_sq; di_sd_i; di_sq_i; dw_m; dw_m_i; dtheta_m];
                Output = f_xu;
                
            elseif CallFlag == 2
          	% ### Call output equation: y = g(x,u)
                g_xu = [i_d; i_q; w; theta; v_dc; i_sd; i_sq; v_sd; v_sq; T_e; T_m; n; theta_e];
                Output = g_xu;
            end
        end
    end
end     % End class definition
