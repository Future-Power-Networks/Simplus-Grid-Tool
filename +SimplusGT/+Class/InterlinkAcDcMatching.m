% Model of an ac-dc converter for inter-linking ac and dc grids.

% Author(s): Yitong Li

%% Notes
%
% The model is in 
% ac-side: load convention, admittance form.
% dc-side: load convention, admittance form.
%
% The model is simply a three phase voltage source inverter, rather than a
% back-to-back topology.

%% Class

classdef InterlinkAcDcMatching < SimplusGT.Class.ModelAdvance
    
    % For temporary use
    properties(Access = protected)
         e_dq;
    end
    
    methods
        % constructor
        function obj = InterlinkAcDcMatching(varargin)
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
            if obj.ApparatusType==2100 || obj.ApparatusType==2101
                State = {'i_d','i_q','theta','v_dc','i'};
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
            W0 = obj.Para(6);
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

            theta = xi;
            
            v = Vg_dc;
            i = P_dc/Vg_dc;
            
            v_dc = Vg_dc - i*R_dc;
            ang_r = 0;
                     
            p_ac = e_d*i_d + e_q*i_q;
            p_dc = v_dc*i;

            % Protected Temp
            obj.e_dq = e_dq;

            % Get equilibrium
        	x_e = [i_d; i_q; theta; v_dc; i];
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
            W0= obj.Para(6);

            R = obj.Para(7);
            K = obj.Para(8);     
            N = obj.Para(9);   
            
           
            L_ac    = xwL_ac/W0;
            R_ac    = xR_ac;
        	L_dc    = xwL_dc/W0;
            R_dc    = xR_dc;
            C_dc    = xC_dc;
            
            % Get states
          	i_d   	= x(1);
         	i_q   	= x(2);

            theta   = x(3);
            v_dc  	= x(4);

            i       = x(5);

            % Get input
        	v_d    = u(1);
            v_q    = u(2);
            v      = u(3);
            ang_r  = u(4);

            % State space equations
            % dx/dt = f(x,u)
            % y     = g(x,u)
            
            V_dc0 = Vg_dc;
            Np = W0/V_dc0;

            if CallFlag == 1    
            % ### Call state equation: dx/dt = f(x,u)
                              
              	% % Ac voltage (duty cycle*v_dc)
                e_dq = obj.e_dq;
                e_d = real(e_dq);
                e_q = imag(e_dq);

                % Dc-link voltage control
                w = (R*(e_d*i_d + e_q*i_q + v_dc*i)/v_dc + K*(v_dc - V_dc0))*N*Np + W0;

                % Ac-side L filter
                di_d = (v_d - R_ac*i_d + w*L_ac*i_q - e_d)/L_ac;
                di_q = (v_q - R_ac*i_q - w*L_ac*i_d - e_q)/L_ac;
                
                dtheta = w;
                
            	% Dc-side inductor
                d_i = (v - v_dc - R_dc*i)/L_dc;
                
                % Dc-side capacitor
                dv_dc = ((e_d*i_d + e_q*i_q)/v_dc + i)/C_dc;
                                
                % Output state
            	f_xu = [di_d; di_q; dtheta; dv_dc; d_i];
                Output = f_xu;
                
            elseif CallFlag == 2
          	% ### Call output equation: y = g(x,u)
                e_dq = obj.e_dq;
                e_d = real(e_dq);
                e_q = imag(e_dq);

                w = (R*(e_d*i_d + e_q*i_q + v_dc*i)/v_dc + K*(v_dc - V_dc0))*N*Np + W0;

                g_xu = [i_d; i_q; i; w; v_dc; theta];
                Output = g_xu;
            end
        end

    end

end     % End class definition