% This class defines the model of a single phase inductor for test.

% Author(s): Yitong Li

classdef PassiveLoad < SimplusGT.Class.ModelAdvance
    properties(Access = protected)
        Connection;
        W0;
      	R;
        L;
    end
    
    methods(Static)
        % Load the parameters foe different types of loads
        function LoadPara(obj)
            % Load basic parameters
            obj.W0 = obj.Para(1);
            obj.Connection = obj.Para(2);
            
            % PQ passive load
            if obj.ApparatusType == 90
                % Load power flow data
                % Notes: P and Q are in load convention
               	P 	= obj.PowerFlow(1);
                Q	= obj.PowerFlow(2);
                V	= obj.PowerFlow(3);
                xi	= obj.PowerFlow(4);
                w   = obj.PowerFlow(5);
                
                % Error check
                if P < 0
                    error(['Error: wrong power flow setting for load, should absorb active power']);
                elseif Q < 0
                    error(['Error: wrong power flow setting for load, should absorb reactive power']);
                end
                
                % Calculate the equivalent passive RL load
                % Series RL connection
                if obj.Connection == 1
                    S = P+j*Q;
                    I = conj(S/V);
                    Z = V/I;
                    obj.R = real(Z);
                    obj.L = imag(Z)/w;
                % Parallel RL connection
                elseif obj.Connection == 2
                    obj.R = V^2/P;
                    obj.L = V^2/Q/w;
                end
            % RL load    
            elseif obj.ApparatusType == 91
                obj.R = obj.Para(2);
                obj.L = obj.Para(3);
            end
        end
        
        % Set the strings of input, output, state
        function SetString(obj)
            obj.LoadPara(obj);
         	obj.InputString  = {'v_d','v_q'};  	% u
        	obj.OutputString = {'i_d','i_q'};  	% y
            if (obj.L ~= 0) 
              	obj.StateString = {'i_Ld','i_Lq'};
            else
                obj.StateString  = {};                  % x
            end
        end
        
        % Calculate the equilibrium
        function Equilibrium(obj)
          	% Get the power PowerFlow values
            P 	= obj.PowerFlow(1);
            Q	= obj.PowerFlow(2);
            V	= obj.PowerFlow(3);
            xi	= obj.PowerFlow(4);
            w   = obj.PowerFlow(5);

            % Calculate
          	i_d = P/V;
            i_q = -Q/V;     % Use -Q because S = V*conj(I)
            v_d = V;
            v_q = 0;
            
            % Get the equilibrium
            if obj.L~=0
                obj.x_e = [i_d;i_q];    % Only for RL series connection
            else
                obj.x_e = [];
            end
            obj.u_e = [v_d;v_q];
            obj.xi = xi;
        end
        
        % State space model
        function [Output] = StateSpaceEqu(obj,x,u,CallFlag)
            
        	% Get state
            if obj.L~= 0
                i_Ld = x(1);
                i_Lq = x(2);
            end
            
            % Get input
            v_d = u(1);
            v_q = u(2);
            
            % Get parameters
            W0 = obj.W0;
            R = obj.R;
            L = obj.L;
            w = W0;

            % State space equations
            if CallFlag == 1
                % State equations: dx/dt = f(x,u)
                if obj.L ~= 0
                    if obj.Connection == 1
                        % v_d = R*i_d + w*L*di_d/dt - w*L*i_q
                        % v_q = R*i_q + w*Ldi_q/dt + w*L*i_d
                        di_Ld = (v_d - R*i_Ld + w*L*i_Lq)/L;
                        di_Lq = (v_q - R*i_Lq - w*L*i_Ld)/L;
                    elseif obj.Connection == 2
                        % v_d = w*L*di_Ld/dt - w*L*i_Lq
                        % v_q = w*Ldi_Lq/dt + w*L*i_Ld
                     	di_Ld = (v_d + w*L*i_Lq)/L;
                        di_Lq = (v_q - w*L*i_Ld)/L;
                    end
                    f_xu = [di_Ld; di_Lq];
                else
                    f_xu = [];
                end
                Output = f_xu;
            elseif CallFlag == 2
                % Output equations: y = g(x,u)
                if obj.L~= 0
                    if obj.Connection == 1
                        i_d = i_Ld;
                        i_q = i_Lq;
                    elseif obj.Connection == 2
                        i_d = i_Ld+v_d/R;
                        i_q = i_Lq+v_q/R;
                    end
                else
                    i_d = v_d/R;
                    i_q = v_q/R;
                end
                g_xu = [i_d; i_q];
                Output = g_xu;              
            end
        end
        
    end
end