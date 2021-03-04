% Model of an ac-dc converter for inter-linking ac and dc grids.

% Author(s): Yitong Li

% This file is NOT ready for use!

%% Notes
%
% The model is in 
% ac-side: load convention, admittance form.
% dc-side: load convention, admittance form.

%% Class

classdef InterlinkAcDc < SimplexPS.Class.ModelAdvance
    
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
            % For the interlink ac-dc converter.
            % The first three inputs must be v_d, v_q, and v; and the first
            % three outputs must be i_d, i_q, and i.
            if (obj.DeviceType == 10) || (obj.DeviceType == 12)
                State = {'i_d','i_q','i_d_i','i_q_i','w_pll_i','w','theta','v_dc','v_dc_i'};
            else
                error('Error: Invalid DeviceType.');
            end
        	Input = {'v_d','v_q','v'};
            Output = {'i_d','i_q','i','w','v_dc','theta'};
        end
        
        function [x_e,u_e,xi] = Equilibrium(obj)
            % Get the power PowerFlow values
            P_ac = obj.PowerFlow{1}(1);
            Q_ac = obj.PowerFlow{1}(2);
            V_ac = obj.PowerFlow{1}(3);
            xi   = obj.PowerFlow{1}(4);
            w    = obj.PowerFlow{1}(5);
            
            P_dc = obj.PowerFlow{2}(1);
            V_dc = obj.PowerFlow{2}(3);
            
            % Get parameters
            C_dc = obj.Para(1);
            V_dc = obj.Para(2);
            L_ac = obj.Para(11);
            R_ac = obj.Para(12);
            
            L_dc =
            R_dc =

            % Calculate paramters
            i_d = P_ac/V_ac;
            i_q = -Q_ac/V_ac;     % Because of conjugate "i"
            v_d = V_ac;
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
            v_dc = V_dc;

            theta = xi;
            
            v = V_dc;
            i = P_dc/V_dc;
            
            % ??? Temp
            obj.i_q_r = i_q_r;

            % Get equilibrium
        	x_e = [i_d; i_q; i_d_i; i_q_i; w_pll_i; w; theta; v_dc; v_dc_i];
        	u_e = [v_d; v_q; v];
        end

        function [Output] = StateSpaceEqu(obj,x,u,CallFlag)
            % Get the power PowerFlow values
            P 	= obj.PowerFlow(1);
            Q	= obj.PowerFlow(2);
            V	= obj.PowerFlow(3);
            xi	= obj.PowerFlow(4);
            w   = obj.PowerFlow(5);
            
           	% Get parameters
            C_dc    = obj.Para(1);
            v_dc_r  = obj.Para(2);
            kp_v_dc = obj.Para(3);
            ki_v_dc = obj.Para(4);
            kp_pll  = obj.Para(5);
            ki_pll  = obj.Para(6);
            tau_pll = obj.Para(7);
            kp_i_dq = obj.Para(8);
            ki_i_dq = obj.Para(9);
            L_ac       = obj.Para(11);
            R_ac       = obj.Para(12);
            L_dc
            R_dc
            
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

            % Get input
        	v_d    = u(1);
            v_q    = u(2);
            v      = u(3);
            
            % State space equations
            % dx/dt = f(x,u)
            % y     = g(x,u)
            if CallFlag == 1    
            % ### Call state equation: dx/dt = f(x,u)
                
                % Current reference
                i_d_r = (v_dc_r - v_dc)*kp_v_dc + v_dc_i;
                i_q_r = obj.i_q_r;    % Constant iq control
               
                % Ac voltage (duty cycle*v_dc)
                e_d = -(i_d_r - i_d)*kp_i_dq + i_d_i;
                e_q = -(i_q_r - i_q)*kp_i_dq + i_q_i;
                
                i_dc = P_dc/v_dc_r;
                
                % Dc-side capacitor
                dv_dc = ((e_d*i_d + e_q*i_q)/v_dc - i_dc)/C_dc; 	% C_dc
                dv_dc_i = (v_dc_r - v_dc)*ki_v_dc;                  % v_dc I
                
                % Current controller
                di_d_i = -(i_d_r - i_d)*ki_i_dq;               	% i_d I
                di_q_i = -(i_q_r - i_q)*ki_i_dq;             	% i_q I
                
                % Ac-side L filter
                di_d = (v_d - R_ac*i_d + w*L_ac*i_q - e_d)/L_ac;      	% L
                di_q = (v_q - R_ac*i_q - w*L_ac*i_d - e_q)/L_ac;      	% L
                
                % PLL
             	e_ang = v_q;   % vq-PLL
                dw_pll_i = e_ang*ki_pll;                    	% PLL I
                dw = (w_pll_i + e_ang*kp_pll - w)/tau_pll;      % PLL tau
                dtheta = w;
                
                % Output state
            	f_xu = [di_d; di_q; di_d_i; di_q_i; dw_pll_i; dw; dtheta; dv_dc; dv_dc_i];
                Output = f_xu;
                
            elseif CallFlag == 2
          	% ### Call output equation: y = g(x,u)
                g_xu = [i_d; i_q; i; w; v_dc; theta];
                Output = g_xu;
            end
        end

    end

end     % End class definition