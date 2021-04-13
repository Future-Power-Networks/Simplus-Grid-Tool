% This class defines the model of grid-following VSI

% Author(s): Yitong Li, Yunjie Gu

%% Notes
%
% The model is in 
% ac-side: load convention, admittance form.
% dc-side: source convention, impedance form.

%% Class

classdef GridFollowingInverterStationary < SimplexPS.Class.ModelAdvance
    
    % For temporary use
    properties(Access = protected)
        i_q_r;
        w_temp;
    end
    
    methods
        % constructor
        function obj = GridFollowingVSI(varargin)
            % Support name-value pair arguments when constructing object
            setProperties(obj,nargin,varargin{:});
        end
    end

    methods(Static)
        
        function [State,Input,Output] = SignalList(obj)
          	% Notes:
            % P_dc is the output power to dc side.
            % The "_i" in "i_d_i", "i_q_i", "v_dc_i" means integral. These
            % states appear because PI controllers are used
            % correspondingly.
            if (obj.DeviceType == 19)
                State = {'i_d','i_q','i_d_i','i_q_i','w_pll_i','w','theta'};
            else
                error('Error: Invalid DeviceType.');
            end
        	Input = {'v_d','v_q','ang_r','P_dc'};
            Output = {'i_d','i_q','w','v_dc','theta'};
        end
        
        function [x_e,u_e,xi] = Equilibrium(obj)
            % Get the power PowerFlow values
            P 	= obj.PowerFlow(1);
            Q	= obj.PowerFlow(2);
            V	= obj.PowerFlow(3);
            xi	= obj.PowerFlow(4);
            w   = obj.PowerFlow(5);

            % Get parameters
            C_dc    = obj.Para(1);
            V_dc    = obj.Para(2);
            L       = obj.Para(11);
            R       = obj.Para(12);
            W0      = obj.Para(13);
            Gi_cd   = obj.Para(14);

            % Calculate paramters
            i_d = P/V;
            i_q = -Q/V;     % Because of conjugate "i"
            v_d = V;
            v_q = 0;
            i_dq = i_d + 1j*i_q;
            v_dq = v_d + 1j*v_q;
            e_dq = v_dq - i_dq * (R + 1j*L*w);
            e_d = real(e_dq);
            e_q = imag(e_dq);
            i_d_i = e_d + Gi_cd*W0*L*(-i_q);
            i_q_i = e_q - Gi_cd*W0*L*(-i_d);
            i_d_r = i_d;
            i_q_r = i_q;
            w_pll_i = w;
            v_dc_i = i_d;
            v_dc = V_dc;
            P_dc = e_d*i_d + e_q*i_q;
            ang_r = 0;
            theta = xi;
            
            % ??? Temp
            obj.i_q_r = i_q_r;
            obj.w_temp = W0;

            % Get equilibrium
            x_e = [i_d; i_q; i_d_i; i_q_i; w_pll_i; w; theta];
        	u_e = [v_d; v_q; ang_r; P_dc];
        end

        function [Output] = StateSpaceEqu(obj,x,u,CallFlag)
            % Get the power PowerFlow values
            P 	= obj.PowerFlow(1);
            Q	= obj.PowerFlow(2);
            V	= obj.PowerFlow(3);
            xi	= obj.PowerFlow(4);
            w   = obj.PowerFlow(5);
            
           	% Get parameters
            kp_pll  = obj.Para(5);      % PLL, P
            ki_pll  = obj.Para(6);      % PLL, I
            kp_i    = obj.Para(8);      % i_dq, P
            ki_i    = obj.Para(9);      % i_dq, I
            L       = obj.Para(11);     % L filter
            R       = obj.Para(12);     % L filter's inner resistance
            W0      = obj.Para(13);
            
            % Get states
          	i_al   	= x(1);
         	i_be   	= x(2);
          	i_d_i  	= x(3);
            i_q_i 	= x(4);
            w_pll_i = x(5);
            w       = x(6);
            theta   = x(7);

            % Get input
        	v_d    = u(1);
            v_q    = u(2);
            ang_r  = u(3);
            P_dc   = u(4);
            
            % Frame transformation
            v_albe = SimplexPS.dq2alphabeta([v_d;v_q],theta);
            i_dq = SimplexPS.alphabeta2dq([i_al;i_be],theta);
            v_al = v_albe(1);
            v_be = v_albe(2);
            i_d = i_dq(1);
            i_q = i_dq(2); 
            
            % State space equations
            % dx/dt = f(x,u)
            % y     = g(x,u)
            if CallFlag == 1    
            % ### Call state equation: dx/dt = f(x,u)
                
                % Get current reference                       
              	i_d_r = P/V;
                i_q_r = obj.i_q_r;    % Constant iq control, PQ/PV node in power flow
                i_albe_r = SimplexPS.alphabeta2dq([i_al;i_be],theta);
                i_al_r = i_albe_r(1);
                i_be_r = i_albe_r(2);
                
             	% Current controller
                di_d_i = -(i_d_r - i_d)*ki_i;               	% i_d I
                di_q_i = -(i_q_r - i_q)*ki_i;                  % i_q I
                
                % Ac voltage (duty cycle*v_dc)
                e_d = -(i_d_r - i_d)*kp_i + i_d_i;
                e_q = -(i_q_r - i_q)*kp_i + i_q_i;

                % Filter inductor
                di_d = (v_d - R*i_d + w*L*i_q - e_d)/L;      	% L
                di_q = (v_q - R*i_q - w*L*i_d - e_q)/L;      	% L
              	
                % PLL
              	% Q = v_q*i_d - v_d*i_q;              % Q-PLL
                Q = e_q*i_d - e_d*i_q;
                if i_d<0
                    e_ang = - Q - ang_r;
                else
                    e_ang = Q - ang_r;
                end
                dw_pll_i = e_ang*ki_pll;                            % integral control of PLL
                if 0                                                                        % ??? 
                    % With PLL LPF
                    dw = (w_pll_i + e_ang*kp_pll - w)/tau_pll;      % PLL filter
                else
                    % Without PLL LPF
                    dw = 0;
                    w = w_pll_i + e_ang*kp_pll;
                end
                obj.w_temp = w;
                dtheta = w;
                
                % Output state
                f_xu = [di_d; di_q; di_d_i; di_q_i; dw_pll_i; dw; dtheta];
                Output = f_xu;
                
            elseif CallFlag == 2
          	% ### Call output equation: y = g(x,u)
                w = obj.w_temp;                                                             % ???
                g_xu = [i_d; i_q; w; v_dc; theta];
                Output = g_xu;
            end
        end

    end

end     % End class definition