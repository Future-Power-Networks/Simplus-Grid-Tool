% This class defines the model of a grid-forming voltage source inverter.

% Author(s): Yitong Li

%% Notes
%
% The model is in load convention.
% 
% The model is in admittance form.
%
% dw means the derivative of w

%% Class

classdef GridFormingVSI < SimplexPS.Class.ModelAdvance
    
    % For temporary use
    properties(Access = protected)
        v_od_r;
        v_oq_r;
        W0;
    end
    
    methods
        % constructor
        function obj = GridFormingVSI(varargin)
            % Support name-value pair arguments when constructing object
            setProperties(obj,nargin,varargin{:});
        end
    end
    
    methods(Static)
        
        function [State,Input,Output] = SignalList(obj)
         	State  = {'i_ld','i_lq','i_ld_i','i_lq_i','v_od','v_oq','v_od_i','v_oq_i','i_od','i_oq','w','theta'};
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
                        xwLf = obj.Para(1);
                        xRf = obj.Para(2);
                        xwCf = obj.Para(3);
                        xwLc = obj.Para(4);
                        xRc = obj.Para(5);
                        xXov = obj.Para(6);
                        xDw = obj.Para(7);
                        xfdroop = obj.Para(8);
                        xfvdc = obj.Para(9);
                        xfidq = obj.Para(10);
                        W0=w;
                        obj.W0=W0;
                        
                        Lf = xwLf/W0;
                        Rf = xRf;
                        Cf = xwCf/W0;
                        Lc = xwLc/W0;
                        Rc = xRc/W0;
                        Xov= xXov;
                        Dw = xDw*W0;
                        wf = xfdroop*2*pi;
                        w_v_odq = xfvdc*2*pi;
                        w_i_ldq = xfidq*2*pi;                   
                                   
%             case 1;  CellPara{row(i)}.Lf      = user_value/W0;
%           	case 2;  CellPara{row(i)}.Rf      = user_value;
%           	case 3;  CellPara{row(i)}.Cf      = user_value/W0;
%            	case 4;  CellPara{row(i)}.Lc  	  = user_value/W0;
%          	case 5;  CellPara{row(i)}.Rc  	  = user_value/W0;
%            	case 6;  CellPara{row(i)}.Xov 	  = user_value;
%             case 7;  CellPara{row(i)}.Dw      = user_value*W0;
%             case 8;  CellPara{row(i)}.wf      = user_value*2*pi;
%           	case 9;  CellPara{row(i)}.w_v_odq = user_value*2*pi;
%           	case 10; CellPara{row(i)}.w_i_ldq = user_value*2*pi;
            % Calculate parameters
            v_gd = V;
            v_gq = 0;
            i_od = P/V;
            i_oq = -Q/V;
            
            v_gdq = v_gd + 1i*v_gq;
            i_odq = i_od + 1i*i_oq;
            v_odq = v_gdq - i_odq*(Rc + 1i*w*Lc);
            i_cdq = v_odq*(1i*w*Cf);
            i_ldq = i_odq - i_cdq;
            e_dq = v_odq - i_ldq*(Rf + 1i*w*Lf);
            
            i_ld = real(i_ldq);
            i_lq = imag(i_ldq);
            i_ld_i = real(e_dq);
            i_lq_i = imag(e_dq);
            v_od = real(v_odq);
            v_oq = imag(v_odq);
            v_od_i = i_ld;
            v_oq_i = i_lq;
            i_od = real(i_od);
            i_oq = imag(i_oq);
            theta = xi;
            
            P0 = 0;
            
            % ??? Temp
            obj.v_od_r = real(v_odq);
            obj.v_oq_r = imag(v_odq);
            
            % Get equilibrium
            x_e = [i_ld; i_lq; i_ld_i; i_lq_i; v_od; v_oq; v_od_i; v_oq_i; i_od; i_oq; w; theta];
            u_e = [v_gd; v_gq; P0];
        end
        
        % State space model
        function [Output] = StateSpaceEqu(obj,x,u,CallFlag)
        	% Get the power PowerFlow values
            V	= obj.PowerFlow(3);
            
            % Get input
            v_gd   = u(1);
            v_gq   = u(2);
            P0     = u(3);

            % Get state
            i_ld   = x(1);
            i_lq   = x(2);
            i_ld_i = x(3);
            i_lq_i = x(4);
            v_od   = x(5);
            v_oq   = x(6); 
            v_od_i = x(7);
            v_oq_i = x(8);
            i_od   = x(9);
            i_oq   = x(10);
            w      = x(11);
            theta  = x(12);
            
            % Get parameter
                        xwLf = obj.Para(1);
                        xRf = obj.Para(2);
                        xwCf = obj.Para(3);
                        xwLc = obj.Para(4);
                        xRc = obj.Para(5);
                        xXov = obj.Para(6);
                        xDw = obj.Para(7);
                        xfdroop = obj.Para(8);
                        xfvdc = obj.Para(9);
                        xfidq = obj.Para(10);
                        w0=obj.W0;
                        W0=obj.W0;
                        Lf = xwLf/W0;
                        Rf = xRf;
                        Cf = xwCf/W0;
                        Lc = xwLc/W0;
                        Rc = xRc/W0;
                        Xov= xXov;
                        Dw = xDw*W0;
                        wf = xfdroop*2*pi;
                        w_v_odq = xfvdc*2*pi;
                        w_i_ldq = xfidq*2*pi; 
%             Lf       = obj.Para(1);
%             Rf       = obj.Para(2);
%             Cf       = obj.Para(3);
%             Lc       = obj.Para(4);
%             Rc       = obj.Para(5);
%             Xov      = obj.Para(6);
%             Dw       = obj.Para(7);
%             wf       = obj.Para(8);
%             w_v_odq  = obj.Para(9);
%             w_i_ldq  = obj.Para(10);
%             w0       = obj.Para(11);
                        
            Rov      = 0;
            kp_i_ldq = w_i_ldq*Lf;
            ki_i_ldq = w_i_ldq^2*Lf/4;
            kp_v_odq = w_v_odq*Cf;
            ki_v_odq = w_v_odq^2*Cf/4*100;
                % This is a different way of setting voltage PI
                % kp_v_odq = 1/(16*w_i_ldq*Lf);
                % ki_v_odq = 1/(4*Lf);
            
            % State space equations
         	% dx/dt = f(x,u)
            % y     = g(x,u)
            if CallFlag == 1    
            % ### Call state equation: dx/dt = f(x,u)
            
                % Power measurement
                p = (v_od*i_od + v_oq*i_oq)*(-1);   % (-1) appears because the model is in load convention
                q = (-v_od*i_oq + v_oq*i_od)*(-1);

                % Droop controller
                % Standard droop equations:
                % PF droop: w = w0 + Dw*(P0 - P*LPF)
                % QV droop: v_od_r = v_od_0 + Dv*(Q0 - Q*LPF)
                %           v_oq_r = v_oq_0
                if 1
                    dw = (w0 + Dw*(P0 - p) - w)*wf;         % P-w droop
                else
                    dw = (w0 - Dw*(P0/V - i_od) - w)*wf; 	% id-w droop
                end
                v_od_r = obj.v_od_r;
                v_oq_r = obj.v_oq_r;
                  	% dv_od_r = (v_od0 + Dv*(Q0 - p) - v_od_r)/Tf;
                    % v_oq_r = v_oq0;
                
                % AC voltage control
                error_v_od = v_od_r - v_od - (i_od*Rov-i_oq*Xov)*(-1);
              	error_v_oq = v_oq_r - v_oq - (i_oq*Rov+i_od*Xov)*(-1);
                i_ld_r = -(error_v_od*kp_v_odq + v_od_i);
                i_lq_r = -(error_v_oq*kp_v_odq + v_oq_i);
                dv_od_i = error_v_od*ki_v_odq;
                dv_oq_i = error_v_oq*ki_v_odq;

                % AC current control
                error_i_ld = i_ld_r-i_ld;
                error_i_lq = i_lq_r-i_lq;
                e_d = -error_i_ld*kp_i_ldq + i_ld_i;
                e_q = -error_i_lq*kp_i_ldq + i_lq_i;
                di_ld_i = -error_i_ld*ki_i_ldq;
                di_lq_i = -error_i_lq*ki_i_ldq;

                % Lf equation
                % e_d - v_od = -(di_ld/dt*Lf + Rf*i_ld - w*Lf*i_lq)
                % e_q - v_oq = -(di_lq/dt*Lf + Rf*i_lq + w*Lf*i_ld)
                di_ld = (v_od - e_d - Rf*i_ld + w*Lf*i_lq)/Lf;
                di_lq = (v_oq - e_q - Rf*i_lq - w*Lf*i_ld)/Lf;

                % Cf equation
                % -(i_ld - i_od) = Cf*dv_cd/dt - w*Cf*v_cq
                % -(i_lq - i_oq) = Cf*dv_cq/dt + w*Cf*v_cd
                dv_od = (-(i_ld - i_od) + w*Cf*v_oq)/Cf;
                dv_oq = (-(i_lq - i_oq) - w*Cf*v_od)/Cf;

                % Lc equation
                % v_od - v_d = -(Lc*di_od/dt + Rc*i_od - w*Lc*i_oq)
                % v_oq - v_q = -(Lc*di_oq/dt + Rc*i_oq + w*Lc*i_od)
                di_od = (v_gd - v_od - Rc*i_od + w*Lc*i_oq)/Lc;
                di_oq = (v_gq - v_oq - Rc*i_oq - w*Lc*i_od)/Lc;

                % Phase angle
                dtheta = w;
                
                % dx
                f_xu = [di_ld; di_lq; di_ld_i; di_lq_i; dv_od; dv_oq; dv_od_i; dv_oq_i; di_od; di_oq; dw; dtheta];
                Output = f_xu;
                
            elseif CallFlag == 2     
            % ### Call output equations: y = g(x,u)
                g_xu = [i_od; i_oq; w; theta];
                Output = g_xu;
            end
            
        end
        
    end
end