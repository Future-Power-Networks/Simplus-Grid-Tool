% This class defines the model of grid-following VSI

% Author(s): Yitong Li, Yunjie Gu

%% Notes
%
% The model is in 
% ac-side: load convention, admittance form.
% dc-side: source convention, impedance form.

%% Class

classdef GridFollowingInverterStationary < SimplusGT.Class.ModelAdvance
    
    % For temporary use
    properties(Access = protected)
        i_q_r;
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
            if (obj.ApparatusType == 19)
                State = {'i_al','i_be','i_al_i','i_be_i','i_al_ii','i_be_ii','w_pll_i','w','theta'};
            else
                error('Error: Invalid ApparatusType.');
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
            ki_i    = obj.Para(9);      % i, I
            L       = obj.Para(11);
            R       = obj.Para(12);
            W0      = obj.Para(13);

            % Calculate paramters
            % dq frame
            i_d = P/V;
            i_q = -Q/V;     % Because of conjugate "i"
            v_d = V;
            v_q = 0;
            i_dq = i_d + 1j*i_q;
            obj.i_q_r = i_q;
            v_dq = v_d + 1j*v_q;
            e_dq = v_dq - i_dq * (R + 1j*L*w);
            e_d = real(e_dq);
            e_q = imag(e_dq);
            w_pll_i = w;
            P_dc = e_d*i_d + e_q*i_q;
            theta = xi;
            ang_r = 0;
            
            % alpha/beta frame
            i_albe = SimplusGT.dq2alphabeta([i_d;i_q],theta);
            i_al = i_albe(1);
            i_be = i_albe(2);
         	e_albe = SimplusGT.dq2alphabeta([e_d;e_q],theta);
            e_al = e_albe(1);
            e_be = e_albe(2);
            
            i_al_ii = e_al/(W0*ki_i);
            i_be_ii = e_be/(W0*ki_i);
            
            % Get equilibrium
            x_e = [i_al; i_be; i_al_i; i_be_i; i_al_ii; i_be_ii; w_pll_i; w; theta];
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
            kp_i    = obj.Para(8);      % i, P
            ki_i    = obj.Para(9);      % i, I
            L       = obj.Para(11);     % L filter
            R       = obj.Para(12);     % L filter's inner resistance
            W0      = obj.Para(13);
            
            % Get states
          	i_al   	= x(1);
         	i_be   	= x(2);
          	i_al_i 	= x(3);
            i_be_i 	= x(4);
            i_al_ii = x(5);
            i_be_ii = x(6);
            w_pll_i = x(7);
            w       = x(8);
            theta   = x(9);

            % Get input
        	v_d    = u(1);
            v_q    = u(2);
            ang_r  = u(3);
            P_dc   = u(4);
            
            % Frame transformation
            v_albe = SimplusGT.dq2alphabeta([v_d;v_q],theta);
            v_al = v_albe(1);
            v_be = v_albe(2);
            i_dq = SimplusGT.alphabeta2dq([i_al;i_be],theta);
            i_d = i_dq(1);
            i_q = i_dq(2); 
            
          	% Get current reference                       
            i_d_r = P/V;
            i_q_r = obj.i_q_r;    % Constant iq control, PQ/PV node in power flow
            i_albe_r = SimplusGT.alphabeta2dq([i_al;i_be],theta);
            i_al_r = i_albe_r(1);
            i_be_r = i_albe_r(2);

            % Current controller
            di_al_ii = i_al_i;
            di_be_ii = i_be_i;
            di_al_i = -(i_al_r - i_al) + W0^2*i_al_ii;
            di_be_i = -(i_be_r - i_be) + W0^2*i_be_ii;

            % Ac voltage (duty cycle*v_dc)
            e_al = -(i_al_r - i_al)*kp_i + ki_i*W0*i_al_ii;
            e_be = -(i_be_r - i_be)*kp_i + ki_i*W0*i_be_ii;
            e_dq = SimplusGT.alphabeta2dq([e_al;e_be],theta);
            e_d = e_dq(1);
            e_q = e_dq(2);

            % Filter inductor
            di_al = (v_al - R*i_d - e_al)/L;      	% L
            di_be = (v_be - R*i_q - e_be)/L;      	% L

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
            dtheta = w;
            
            % State space equations
            % dx/dt = f(x,u)
            % y     = g(x,u)
            if CallFlag == 1    
            % ### Call state equation: dx/dt = f(x,u)
                f_xu = [di_al; di_be; di_al_i; di_be_i; di_al_ii; di_be_ii; dw_pll_i; dw; dtheta];
                Output = f_xu;
            elseif CallFlag == 2
          	% ### Call output equation: y = g(x,u)
                g_xu = [i_d; i_q; w; v_dc; theta];
                Output = g_xu;
            end
        end

    end

end     % End class definition