% Model of a grid-following buck dc-dc converter.

% Author(s): Yitong Li

%% Notes
%
% Comprehensive comments will be added for this model here with the first
% priority, so that this model can be a specific example of a dc apparatus.
%
% The model is in 
% dc-grid-side: load convention, admittance form.
% dc-load/source-side: source convention, impedance form.
%
% For the buck dc-dc converter, the dc generation/load side suffers
% discontinous current, and the dc grid side is in continous current
% because of the output inductor.

%% Class

classdef GridFeedingBuck < SimplusGT.Class.ModelAdvance

    methods
        % constructor
        function obj = GridFeedingBuck(varargin)
            % Support name-value pair arguments when constructing object
            setProperties(obj,nargin,varargin{:});
        end
    end

    methods(Static)
        
        function [State,Input,Output] = SignalList(obj)
        	if obj.ApparatusType == 1010
                State = {'i','i_i'};
            elseif obj.ApparatusType == 1011
                State = {'i','i_i','v_dc','v_dc_i'};
            else
                error('Error: Invalid ApparatusType.');
            end
        	Input = {'v','P_dc'};
            Output = {'i','v_dc'};
        end
        
        function [x_e,u_e,xi] = Equilibrium(obj)
            % Get the power PowerFlow values
            P 	= obj.PowerFlow(1);
            V	= obj.PowerFlow(3);
            xi	= 0;

            % Get parameters
                                  
            V_dc = obj.Para(1);
            R    = obj.Para(4);

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
            if obj.ApparatusType == 1010
                x_e = x_e_1;
            elseif obj.ApparatusType == 1011
                x_e = [x_e_1; v_dc; v_dc_i];
            end
        	u_e = [v; P_dc];
        end

        function [Output] = StateSpaceEqu(obj,x,u,CallFlag)
            % Get the power PowerFlow values
            P 	= obj.PowerFlow(1);
            V	= obj.PowerFlow(3);
            
           	% Get parameters
            
            xVdc  = obj.Para(1);
            xCdc  = obj.Para(2);
            xwL   = obj.Para(3);
            xR    = obj.Para(4);
            xfi   = obj.Para(5);
            xfvdc = obj.Para(6);
            W0    = obj.Para(7);
            
            w_vdc	= xfvdc*2*pi; 	% (rad/s) bandwidth, vdc
            w_i     = xfi*2*pi;     % (rad/s) bandwidth, i
            L       = xwL/W0;     	% L filter
            R       = xR;           % L filter's inner resistance
            C_dc    = xCdc;
            v_dc_r  = xVdc;
            kp_v_dc = xVdc*xCdc*w_vdc;      % v_dc, P
            ki_v_dc = kp_v_dc*w_vdc/4;      % v_dc, I
            kp_i    = L * w_i;              % i_dq, P
            ki_i    = L * w_i^2 /4;         % i_dq, I   
            
            % Get states
          	i   	= x(1);
          	i_i  	= x(2);
            if obj.ApparatusType == 1011
                v_dc  	= x(3);
                v_dc_i 	= x(4);
            elseif obj.ApparatusType == 1010
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
               	if obj.ApparatusType == 1011
                    % DC-link control
                    i_r = (v_dc_r - v_dc)*kp_v_dc + v_dc_i;
                elseif obj.ApparatusType == 1010
                    % Direct current control                                           
                    i_r = P/V;
                end
                
                % Converter leg node voltage (duty cycle*v_dc)
                e = -(i_r - i)*kp_i + i_i;
                        
                % State equations
                if obj.ApparatusType == 1011
                    dv_dc = (e*i - P_dc)/v_dc/C_dc;         % C_dc
                    dv_dc_i = (v_dc_r - v_dc)*ki_v_dc;     	% v_dc I
                elseif obj.ApparatusType == 1010
                    % No dc link control
                end
                di_i = -(i_r - i)*ki_i;
                di = (v - R*i - e)/L;

                % Output state
                f_xu_1 = [di; di_i];
                if obj.ApparatusType == 1011
                    f_xu = [f_xu_1; dv_dc; dv_dc_i];
                elseif obj.ApparatusType == 1010
                    f_xu = f_xu_1;
                end
                Output = f_xu;
                
            elseif CallFlag == 2
          	% ### Call output equation: y = g(x,u)
                g_xu = [i; v_dc];
                Output = g_xu;
            end
        end

    end

end     % End class definition