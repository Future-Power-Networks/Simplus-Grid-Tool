% Model of a grid-following buck dc-dc converter.

% Author(s): Yitong Li

%% Notes
%
% The model is in 
% dc-grid-side: load convention, admittance form.
% dc-load/source-side: source convention, impedance form.
%
% For the buck dc-dc converter, the dc generation/load side suffers
% discontinous current, and the dc grid side is in continous current
% because of the output inductor.

%% Class

classdef GridFollowingBuck < SimplexPS.Class.ModelAdvance

    methods
        % constructor
        function obj = GridFollowingBuck(varargin)
            % Support name-value pair arguments when constructing object
            setProperties(obj,nargin,varargin{:});
        end
    end

    methods(Static)
        
        function [State,Input,Output] = SignalList(obj)
            if (obj.DeviceType == 10)
                State = {'i','i_i','v_dc','v_dc_i'};
            elseif obj.DeviceType == 11
                State = {'i','i_i'};
            else
                error('Error: Invalid DeviceType.');
            end
        	Input = {'v','P_dc'};
            Output = {'i','v_dc','w','theta'};
        end
        
        function [x_e,u_e,xi] = Equilibrium(obj)
            % Get the power PowerFlow values
            P 	= obj.PowerFlow(1);
            V	= obj.PowerFlow(3);
            xi	= 0;

            % Get parameters
            V_dc = obj.Para(2);
            R    = obj.Para(8);

            % Calculate paramters
            i = P/V;
            v = V;
            e = v - i * R;
            i_i = e;
            v_dc_i = i;
            v_dc = V_dc;
            P_dc = e*i;

            % Get equilibrium
            x_e_1 = [i; i_i];
            if obj.DeviceType == 1000
                x_e = x_e_1;
            elseif obj.DeviceType == 1001
                x_e = [x_e_1; v_dc; v_dc_i];
            end
        	u_e = [v_d; v_q; P_dc];
        end

        function [Output] = StateSpaceEqu(obj,x,u,CallFlag)
            % Get the power PowerFlow values
            P 	= obj.PowerFlow(1);
            V	= obj.PowerFlow(3);
            xi	= 0;
            w   = 0;
            
           	% Get parameters
            C_dc    = obj.Para(1);
            v_dc_r  = obj.Para(2);
            kp_v_dc = obj.Para(3);      % v_dc, P
            ki_v_dc = obj.Para(4);      % v_dc, I
            kp_i    = obj.Para(5);      % i_dq, P
            ki_i    = obj.Para(6);      % i_dq, I
            L       = obj.Para(7);      % L filter
            R       = obj.Para(8);      % L filter's inner resistance
            W0      = obj.Para(9);
            
            % Get states
          	i   	= x(1);
          	i_i  	= x(2);
            if obj.DeviceType == 1011
                v_dc  	= x(3);
                v_dc_i 	= x(4);
            elseif obj.DeviceType == 1010
                v_dc    = v_dc_r;
                v_dc_i  = 0;
            end

            % Get input
        	v    = u(1);
            P_dc = u(2);
            
            % State space equations
            % dx/dt = f(x,u)
            % y     = g(x,u)
            if CallFlag == 1    
            % ### Call state equation: dx/dt = f(x,u)
                
                % Get current reference
               	if obj.DeviceType == 1011
                    % DC-link control
                    i_r = (v_dc_r - v_dc)*kp_v_dc + v_dc_i;
                elseif obj.DeviceType == 1010
                    % Direct current control                                           
                    i_r = P/V;
                end
                
                % Converter leg node voltage (duty cycle*v_dc)
                e = -(i_r - i)*kp_i + i_i;
                        
                % State equations
                if obj.DeviceType == 1011
                    dv_dc = (e*i - P_dc)/v_dc/C_dc; 	% C_dc
                    dv_dc_i = (v_dc_r - v_dc)*ki_v_dc;             	% v_dc I
                elseif obj.DeviceType == 1010
                    % No dc link control
                end
                di_i = -(i_r - i_d)*ki_i;
                di = (v - R*i - e)/L;

                % Output state
                f_xu_1 = [di; di_i];
                if obj.DeviceType == 1011
                    f_xu = [f_xu_1; dv_dc; dv_dc_i];
                elseif obj.DeviceType == 1010
                    f_xu = f_xu_1;
                end
                Output = f_xu;
                
            elseif CallFlag == 2
          	% ### Call output equation: y = g(x,u)
                w = 0;
                theta = 0;
                g_xu = [i; v_dc;  w; theta];
                Output = g_xu;
            end
        end

    end

end     % End class definition