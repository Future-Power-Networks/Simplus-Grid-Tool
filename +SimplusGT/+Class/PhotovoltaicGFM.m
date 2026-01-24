 % This class defines the model of a photovoltaic system with GFM(Matching control) control.

% Author(s): Wenjie Ning
%


%% Notes
%
% The model is in load convention.
% 
% The model is in admittance form.
%
% dw means the derivative of w

%% Class

classdef PhotovoltaicGFM < SimplusGT.Class.ModelAdvance
    
    % For temporary use
    properties(Access = protected)
        v_od_r;
        v_oq_r;
        P0;
        Q0;
    end
    
    methods
        % constructor
        function obj = PhotovoltaicGFM(varargin)
            % Support name-value pair arguments when constructing object
            setProperties(obj,nargin,varargin{:});
        end
    end
    
    methods(Static)
        
        function [State,Input,Output] = SignalList(obj)
            State  = {'v_i','i_i','i_l','v_pv','v_d_i','v_q_i','i_d_i','i_q_i','i_d','i_q','v_d','v_q','i_gd','i_gq','v_o','theta'};  
            Input  = {'v_d','v_q','P0'};
            Output = {'i_d','i_q','w','theta'};
        end
        
        % Calculate the equilibrium
        function [x_e,u_e,xi] = Equilibrium(obj)
            % Get the power PowerFlow values
            P 	= obj.PowerFlow(1);
            Q	= obj.PowerFlow(2);
            V	= obj.PowerFlow(3);
            xi	= obj.PowerFlow(4);
            w   = obj.PowerFlow(5);
            % Get parameter
            xwLf     = obj.Para(1);
            Rf       = obj.Para(2);
            xwCf     = obj.Para(3);
            xwLc     = obj.Para(4);
            Rc       = obj.Para(5);
            Xov      = obj.Para(6);
            Rov      = 0;
            W0       = obj.Para(11);
            v_pv_ref = obj.Para(15);
            v_dc_ref = obj.Para(16);

            % Calculate parameters
            Lf_ac    = xwLf/W0;
            Cf_ac    = xwCf/W0;
            Lc_ac    = xwLc/W0;
            Rf_dc    = 0.01;
            
            v_gd = V;
            v_gq = 0;
            i_gd = P/V;
            i_gq = -Q/V;
            
            v_gdq = v_gd + 1i*v_gq;
            i_gdq = i_gd + 1i*i_gq;
            v_dq  = v_gdq - i_gdq*(Rc + 1i*w*Lc_ac);
            i_cdq = v_dq*(1i*w*Cf_ac);
            i_dq  = i_gdq - i_cdq;
            e_dq  = v_dq - i_dq*(Rf + 1i*w*Lf_ac);
            
            i_d   = real(i_dq);
            i_q   = imag(i_dq);
            i_d_i = -real(e_dq);         
            i_q_i = -imag(e_dq);
            v_d   = real(v_dq);
            v_q   = imag(v_dq);
            v_d_i = -i_d;              
            v_q_i = -i_q;
            i_gd  = real(i_gdq);
            i_gq  = imag(i_gdq);
            theta = xi;

            % dc
            v_o   = v_dc_ref;
            v_pv  = v_pv_ref;
            i_pv     = -P/v_pv;

            i_l   = i_pv;
            v_i   = i_l; 
            ed_dc = v_pv - i_l*Rf_dc;
            i_i   = ed_dc;

            obj.P0 = P*(-1);
            obj.Q0 = Q*(-1);
            
            v_odq_r    = v_dq + (Rov + 1i*Xov)*i_gdq*(-1);
            v_od_r     = real(v_odq_r);
            v_oq_r     = imag(v_odq_r);
            obj.v_od_r = v_od_r;
            obj.v_oq_r = v_oq_r;
            
            % Get equilibrium
            x_e = [v_i; i_i; i_l; v_pv; v_d_i; v_q_i; i_d_i; i_q_i; i_d; i_q; v_d; v_q; i_gd; i_gq; v_o; theta];
            u_e = [v_gd; v_gq; P];
        end
      
        % State space model
        function [Output] = StateSpaceEqu(obj,x,u,CallFlag)
            % Set the invertor control type
            Control_type = 1;  %% 1--Closed loop 0--Open loop

            % Get the power PowerFlow values
            V	= obj.PowerFlow(3);
            
            % Get input
            vgd   = u(1);
            vgq   = u(2);
            P     = u(3);
            % Get state
            v_i    = x(1);
            i_i    = x(2);
            i_l    = x(3);
            v_pv   = x(4);
            v_d_i  = x(5);
            v_q_i  = x(6); 
            i_d_i  = x(7);
            i_q_i  = x(8);
            i_d    = x(9);
            i_q    = x(10);
            v_d    = x(11);
            v_q    = x(12);
            i_gd   = x(13);
            i_gq   = x(14);
            v_o    = x(15);
            theta  = x(16);
            
            % Get parameters
            xwLf     = obj.Para(1);
            Rf       = obj.Para(2);
            xwCf     = obj.Para(3);
            xwLc     = obj.Para(4);
            Rc       = obj.Para(5);
            Xov      = obj.Para(6);
            Rov      = 0;
            N1       = obj.Para(7);        
            R1       = obj.Para(8);            
            xfvdq    = obj.Para(9);
            xfidq    = obj.Para(10);
            W0       = obj.Para(11);
            S_air    = obj.Para(12);
            T_air    = obj.Para(13);
            C_o      = obj.Para(14); 
            v_pv_ref = obj.Para(15);
            v_dc_ref = obj.Para(16);
            xfvdc    = obj.Para(17);
            xfidc    = obj.Para(18);
 
            % Solar para
            T_ref = 25;
            S_ref = 1000;
            I_sc  = 14.880;
            I_m   = 13.88;
            U_m   = 576;
            U_oc  = 708;
            a     = 0.0025;
            b     = 0.00288;
            Rs    = 0.5;
            k     = 0.03;
            C2    = (U_m/U_oc - 1)/(log(1 - I_m/I_sc));
            C1    = (1 - I_m/I_sc)*exp(-U_m/(C2*U_oc));
            T     = T_air + k*S_air;
            dT    = T - T_ref;
            dI    = a * S_air/S_ref * dT + I_sc*(1 - S_ref/S_air);
            dU    = -b * dT - Rs * dI;

            i_pv_ref = -P/v_pv_ref;
            A_mpp    = i_pv_ref/(( I_sc * (1 - C1 * (exp((v_pv_ref*800 + dU)/C2/U_oc) - 1))) + dI);
            % Update paramters
            Lf_ac = xwLf/W0;
            Cf_ac = xwCf/W0;
            Lc_ac = xwLc/W0;

            Rf_dc = 0.01;
            Lf_dc = 0.05/W0;
            Cf_dc = 0.02/W0;

            w_v_dc  = xfvdc*2*pi; 
            kp_v_dc = -Cf_dc*w_v_dc* 20;
            ki_v_dc = -Cf_dc*w_v_dc^2/4* 20;

            w_i_dc  = xfidc*2*pi;
            kp_i_dc = -Lf_dc*w_i_dc;
            ki_i_dc = -Lf_dc*(w_i_dc^2)/4;

            w_v_ac  = xfvdq*2*pi;    % Voltage Controller
            kp_v_ac = Cf_ac*w_v_ac;
            ki_v_ac = Cf_ac*w_v_ac^2/4* 20;

            w_i_ac  = xfidq*2*pi;  % Current Controller
            kp_i_ac = Lf_ac*w_i_ac;
            ki_i_ac = Lf_ac*(w_i_ac^2)/4;
       
  
            % State space equations
            % dx/dt = f(x,u)
            % y     = g(x,u)
            if CallFlag == 1    
            % ### Call state equation: dx/dt = f(x,u)
                i_pv       = (( I_sc * (1 - C1 * (exp((v_pv*800 + dU)/C2/U_oc) - 1))) + dI)*A_mpp;
                error_v_pv = v_pv_ref - v_pv;
                i_r        = kp_v_dc*error_v_pv + v_i;
                dv_i       = ki_v_dc *error_v_pv;

                error_i_l  = i_r - i_l;
                ed_dc      = kp_i_dc*error_i_l + i_i;
                di_i       = ki_i_dc*error_i_l;

                di_l       = (v_pv - ed_dc - Rf_dc*i_l)/Lf_dc;
                dv_pv      = (i_pv - i_l)/Cf_dc;
                
                v_d_r = obj.v_od_r;
                v_q_r = obj.v_oq_r;
                

                if Control_type ==1
                    % AC voltage control
                    error_vd = v_d_r - v_d - (i_gd*Rov - i_gq*Xov) * (-1);
              	    error_vq = v_q_r - v_q - (i_gq*Rov + i_gd*Xov) * (-1);
                    i_d_r    = -(error_vd*kp_v_ac + v_d_i);
                    i_q_r    = -(error_vq*kp_v_ac + v_q_i);
                    dv_d_i   = error_vd*ki_v_ac;
                    dv_q_i   = error_vq*ki_v_ac;
    
                    % AC current control
                    error_id = i_d_r - i_d;
                    error_iq = i_q_r - i_q;
                    e_d      = -(error_id*kp_i_ac + i_d_i);
                    e_q      = -(error_iq*kp_i_ac + i_q_i);
                    di_d_i   = error_id*ki_i_ac;            
                    di_q_i   = error_iq*ki_i_ac;
                else
                    dv_d_i = 0;
                    dv_q_i = 0;
                    di_d_i = 0;
                    di_q_i = 0;
                    e_d    =v_d_r;
                    e_q    =v_q_r;
                end

                P_r = i_l * ed_dc;
                p   = (e_d*i_d + e_q*i_q)*(-1);
                w   = ((P_r - p)/v_o*R1 + (v_o - v_dc_ref))*N1*W0 + W0; 

                % dtheta = w-W0;
                dtheta = w;
                dv_o   = (P_r - p)/v_o/C_o; 

                % Lf equation
                di_d = (v_d - e_d - Rf*i_d + w*Lf_ac*i_q)/Lf_ac;
                di_q = (v_q - e_q - Rf*i_q - w*Lf_ac*i_d)/Lf_ac;

                % Cf equation
                dv_d = (-(i_d - i_gd) + w*Cf_ac*v_q)/Cf_ac;
                dv_q = (-(i_q - i_gq) - w*Cf_ac*v_d)/Cf_ac;

                % Lc equation
                di_gd = (vgd - v_d - Rc*i_gd + w*Lc_ac*i_gq)/Lc_ac;
                di_gq = (vgq - v_q - Rc*i_gq - w*Lc_ac*i_gd)/Lc_ac;
                
                % dx
                f_xu   = [dv_i; di_i; di_l; dv_pv; dv_d_i; dv_q_i; di_d_i; di_q_i; di_d; di_q; dv_d; dv_q; di_gd; di_gq; dv_o; dtheta];
                Output = f_xu;
                
            elseif CallFlag == 2    
                error_v_pv = v_pv_ref - v_pv;
                i_r        = kp_v_dc*error_v_pv + v_i;

                error_i_l  = i_r - i_l;
                ed_dc      = kp_i_dc*error_i_l + i_i;
                
                v_d_r = obj.v_od_r;
                v_q_r = obj.v_oq_r;
                
                if Control_type ==1
                    % AC voltage control
                    error_vd = v_d_r - v_d - (i_gd*Rov-i_gq*Xov) * (-1);
              	    error_vq = v_q_r - v_q - (i_gq*Rov+i_gd*Xov) * (-1);
                    i_d_r = -(error_vd*kp_v_ac + v_d_i);
                    i_q_r = -(error_vq*kp_v_ac + v_q_i);
                    
                    % AC current control
                    error_id = i_d_r - i_d;
                    error_iq = i_q_r - i_q;
                    e_d = -(error_id*kp_i_ac + i_d_i);
                    e_q = -(error_iq*kp_i_ac + i_q_i);
                else
                    e_d   = v_d_r;
                    e_q   = v_q_r;
                end

                % end
                P_r = i_l * ed_dc;
                p   = (e_d*i_d + e_q*i_q)*(-1);
                w   = ((P_r - p)/v_o*R1 + (v_o - v_dc_ref))*N1*W0 + W0; 

                % dx
                g_xu   = [i_gd; i_gq; w; theta];
                Output = g_xu;
            end
            
        end
        
    end
end




