% Model of an ac-dc converter for inter-linking ac and dc grids.

% Author(s): Yitong Li

% This file is NOT ready for use!

%% Notes
%
% The model is in 
% ac-side: load convention, admittance form.
% dc-side: load convention, admittance form.
%
% The model is simply a three phase voltage source inverter, rather than a
% back-to-back topology.

%% Class

classdef InterlinkAcDc < SimplusGT.Class.ModelAdvance
    
    % For temporary use
    properties(Access = protected)
        i_q_r;
    end
    
    methods
        % constructor
        function obj = InterlinkAcDc(varargin)
            % Support name-value pair arguments when constructing object
            setProperties(obj,nargin,varargin{:});
        end
    end

    methods(Static)
        
        function [State,Input,Output] = SignalList(obj)
          	% Notes:
            %
            % The first three inputs must be v_d, v_q, and v; and the first
            % three outputs must be i_d, i_q, and i.
            %
            % v_d, v_q, i_d, i_q are ac-grid electrical port.v, i are the
            % dc-grid electrical port.
            %
            % v_dc is the dc link voltage, there is an inductor between v
            % and v_dc. This inductor makes the system admittance model
            % proper seen from the dc side.
            if obj.ApparatusType==2000 || obj.ApparatusType==2001
                State = {'i_d','i_q','i_d_i','i_q_i','w_pll_i','w','theta','v_dc','v_dc_i','i'};
            else
                error('Error: Invalid ApparatusType.');
            end
        	Input = {'v_d','v_q','v','ang_r'};
            Output = {'i_d','i_q','i','w','v_dc','theta'};
        end
        
        function [x_e,u_e,xi] = Equilibrium(obj)
            % Get the power PowerFlow values
            % The interlink converter has two sets of power flow results
            % because it is connected to two buses: the first one is ac,
            % and the second one is dc.
            P_ac    = obj.PowerFlow(1);
            Q_ac    = obj.PowerFlow(2);
            Vg_ac   = obj.PowerFlow(3);
            xi      = obj.PowerFlow(4);
            w       = obj.PowerFlow(5);
            
            P_dc    = obj.PowerFlow(6);
            Vg_dc   = obj.PowerFlow(8);
            
            % Get parameters
            wL_ac = obj.Para(2);
            W0 = obj.Para(9);
            L_ac  = wL_ac/W0;
            R_ac  = obj.Para(3);
            R_dc  = obj.Para(5);

            % Calculate paramters
            i_d = P_ac/Vg_ac;
            i_q = -Q_ac/Vg_ac;     % Because of conjugate "i"
            v_d = Vg_ac;
            v_q = 0;
            i_dq = i_d + 1j*i_q;
            v_dq = v_d + 1j*v_q;
            e_dq = v_dq - i_dq * (R_ac + 1j*L_ac*w);
            e_d = real(e_dq);
            e_q = imag(e_dq);
            i_d_i = e_d;
            i_q_i = e_q;
            i_d_r = i_d;
            i_q_r = i_q;
            w_pll_i = w;
            v_dc_i = i_d_r;

            theta = xi;
            
            v = Vg_dc;
            i = P_dc/Vg_dc;
            
            v_dc = Vg_dc + i*R_dc;
            
            ang_r = 0;
            
            % ??? Temp
            obj.i_q_r = i_q_r;

            % Get equilibrium
        	x_e = [i_d; i_q; i_d_i; i_q_i; w_pll_i; w; theta; v_dc; v_dc_i; i];
        	u_e = [v_d; v_q; v; ang_r];
        end

        function [Output] = StateSpaceEqu(obj,x,u,CallFlag)
            % Get power flow
            P_ac    = obj.PowerFlow(1);
            Vg_ac   = obj.PowerFlow(3);
            Vg_dc   = obj.PowerFlow(8);
            
           	% Get parameters
            xC_dc = obj.Para(1);
            xwL_ac= obj.Para(2);
            xR_ac= obj.Para(3);
            xwL_dc= obj.Para(4);
            xR_dc= obj.Para(5);
            xfidq= obj.Para(6);
            xfvdc= obj.Para(7);
            xfpll= obj.Para(8);
            W0= obj.Para(9);
            
            V_dc = 1;
            w_vdc = xfvdc*2*pi; 	% (rad/s) bandwidth, vdc
            w_pll     = xfpll*2*pi;  	% (rad/s) bandwidth, pll
            w_i     = xfidq*2*pi; 	% (rad/s) bandwidth, idq
            w_tau_pll = 200*2*pi;   % (rad/s) PLL filter bandwidth
            
            L_ac    = xwL_ac/W0;
            R_ac    = xR_ac;
        	L_dc    = xwL_dc/W0;
            R_dc    = xR_dc;
            C_dc    = xC_dc;
            kp_v_dc = V_dc*xC_dc*w_vdc;
            ki_v_dc = kp_v_dc*w_vdc/4;
            kp_i_dq = L_ac * w_i;
            ki_i_dq = L_ac * w_i^2 /4;
            kp_pll = w_pll;
            ki_pll = kp_pll * w_pll/4;
            tau_pll = 1/w_tau_pll;
            
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
            i       = x(10);

            % Get input
        	v_d    = u(1);
            v_q    = u(2);
            v      = u(3);
            ang_r  = u(4);
            
            % State space equations
            % dx/dt = f(x,u)
            % y     = g(x,u)
            if CallFlag == 1    
            % ### Call state equation: dx/dt = f(x,u)
                
                % Dc-link voltage control
                v_dc_r = Vg_dc;
                if obj.ApparatusType==2000
                    dv_dc_i = 0;
                elseif obj.ApparatusType==2001
                    dv_dc_i = (v_dc_r - v_dc)*ki_v_dc;
                end
                
                % Ac current control
                if obj.ApparatusType==2000
                    i_d_r = P_ac/Vg_ac;
                elseif obj.ApparatusType==2001
                    i_d_r = (v_dc_r - v_dc)*kp_v_dc + v_dc_i;
                end
                i_q_r = obj.i_q_r;                  % Constant iq
                di_d_i = -(i_d_r - i_d)*ki_i_dq;
                di_q_i = -(i_q_r - i_q)*ki_i_dq;
                
              	% Ac voltage (duty cycle*v_dc)
                e_d = -(i_d_r - i_d)*kp_i_dq + i_d_i;
                e_q = -(i_q_r - i_q)*kp_i_dq + i_q_i;
                
                % Ac-side L filter
                di_d = (v_d - R_ac*i_d + w*L_ac*i_q - e_d)/L_ac;
                di_q = (v_q - R_ac*i_q - w*L_ac*i_d - e_q)/L_ac;
                
                % PLL
             	e_ang = v_q - ang_r;
                dw_pll_i = e_ang*ki_pll;
                dw = (w_pll_i + e_ang*kp_pll - w)/tau_pll;
                dtheta = w;
                
            	% Dc-side inductor
                d_i = (v - v_dc - R_dc*i)/L_dc;
                
                % Dc-side capacitor
                dv_dc = ((e_d*i_d + e_q*i_q)/v_dc + i)/C_dc;
                
                % Output state
            	f_xu = [di_d; di_q; di_d_i; di_q_i; dw_pll_i; dw; dtheta; dv_dc; dv_dc_i; d_i];
                Output = f_xu;
                
            elseif CallFlag == 2
          	% ### Call output equation: y = g(x,u)
                g_xu = [i_d; i_q; i; w; v_dc; theta];
                Output = g_xu;
            end
        end

    end

end     % End class definition