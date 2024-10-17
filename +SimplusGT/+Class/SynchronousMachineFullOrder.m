% This class defines the model of synchronous machine with full order.

% Author(s): Duange Guo

%% Notes
%
% The model is in admittance form.
%
% dw means the derivative of w
%
% This model contains 6 windings equations:
% stator: d and q axis
% rotor: d-axis field winding (fd), d-axis damper (1d), 2 q-axis dampers (gq, kq), 
%
% 1. q-axis leads the d-axis (IEEE standard)
% 2. rotor angle with regard to q-axis
% 3. the whole model is in load convention, 
% but the equations below are in generator convention 
% 4. 0-axis is omitted
% 5. w is not in p.u.
%

%% Class

classdef SynchronousMachineFullOrder < SimplusGT.Class.ModelAdvance
    
  	methods
        % Constructor
        function obj = SynchronousMachineFullOrder(varargin)
            setProperties(obj,nargin,varargin{:});
        end
    end
    
    methods(Static)
        
        function [State,Input,Output] = SignalList(obj)
        	State  = {'psi_d','psi_q',...    % d-q axis of stator windings
                      'psi_fd','psi_1d',...  % d-axis of rotor windings
                      'psi_gq','psi_kq',...  % q-axis of rotor windings
                      'w','theta',...        % rotation states
                      'T_m',...
                      'V_f','R_f','E_fd',... % exciter related
                      'Px_1','Px_2','Px_3'}; % PSS related
            Input  = {'v_d','v_q',...        % input voltage
                      'T_0',...              % input torque
                      'w_ref',...            % governor speed ref
                      'U_ref'};              % exciter voltage ref
            Output = {'i_d','i_q',...        % output current
                      'w',...                % angular rotation speed (rad/s)
                      'theta',...            % rotor angle
                      'T_e',...              % EM torque
                      'U_mag',...            % voltage magnitude
                      'psi_d',...
                      'psi_q',...
                      'psi_fd','psi_1d',...
                      'psi_gq','psi_kq'}; 
        end


        % Calculate the equilibrium
        % The equilibrium is determined by the power flow data and apparatus's
        % own paramters.This function will be called once, at the
        % beginneing of simulation.
        function [x_e,u_e,xi] = Equilibrium(obj)
         	% Get the power PowerFlow values
            P 	= obj.PowerFlow(1);
            Q	= obj.PowerFlow(2);
            V	= obj.PowerFlow(3);
            xi	= obj.PowerFlow(4);
            w   = obj.PowerFlow(5);
            
          	% Get parameters in pu
            x_d   = obj.Para(1);   %
            x_dd  = obj.Para(2);   % X_dd represents Xd'
            x_ddd = obj.Para(3);   %
            x_q   = obj.Para(4);   %
            x_qd  = obj.Para(5);   % X_qd represents Xq'
            x_qdd = obj.Para(6);   %
            x_l   = obj.Para(7);   % or x_2, x_ls in some practical cases
            R_a   = obj.Para(8);   % 
            T_dd  = obj.Para(9);   % T_dd represents Td'
            T_ddd = obj.Para(10);  %
            T_qd  = obj.Para(11);  % T_qd represents Tq'
            T_qdd = obj.Para(12);  %
            H     = obj.Para(13);  % Inertial constant
            S_A   = obj.Para(14);  % Saturation parameters
            S_B   = obj.Para(15);
            S_N   = obj.Para(16);
            D     = obj.Para(17);            % Damp constant
            SpeedPara   = obj.Para(18);      % simplified governor parameters
            K_A   = obj.Para(19);            % simplified exciter parameters
            EnableSaturation = obj.Para(20); % Saturation setting
            wb    = obj.Para(21);
            
            T_E = 0.314;  % Time constant of exciter
            T_A = 0.02;   % Time constant of AV
            T_F = 0.35;   % Time constant of Feedback
            K_E = 1.0;    % Exciter Gain 
            K_F = 0.063;  % Feedback Gain
            % K_A = 20;   % K_A can be adjusted

            % Calculated parameters
            xd_diff = x_d - x_dd;
            xdd_diff = x_dd - x_ddd;
            xdl_diff = x_dd-x_l;
            xddl_diff = x_ddd-x_l;
            
            xq_diff = x_q-x_qd;
            xqd_diff = x_qd-x_qdd;
            xql_diff = x_qd-x_l;
            xqql_diff = x_qdd-x_l;

            xd_frac = xdd_diff/xdl_diff^2;
            xq_frac = xqd_diff/xql_diff^2;
            
            % Calculate equilibrium
            % We calculate the equilibrium by solving the algebra equations
            % directly. By doing so, we need to import symbol variable, which
            % can only be applied in the simulink when the interpretable
            % execution is used! (setting in the MATLAB system)

            S_D0 = -P/V;  % The equations below are in generator convention
            S_Q0 = -Q/V;
            S_DQ0 = S_D0 + 1j*S_Q0;
            i_abs = abs(S_DQ0);
            ui_argdiff = angle(S_DQ0);
            i_arg = xi - ui_argdiff;

            i_DQ = i_abs * exp(1i * i_arg);
            i_Q  = imag(i_DQ);
            i_D = real(i_DQ);

            v_DQ = V * exp(1i * xi);
            v_D = real(v_DQ);
            v_Q = imag(v_DQ);

            syms theta
            T_rotate = [cos(pi/2-theta) -sin(pi/2-theta); sin(pi/2-theta) cos(pi/2-theta)];
            v_dq = T_rotate*[v_D; v_Q];
            i_dq = T_rotate*[i_D; i_Q];

            v_d = v_dq(1);
            v_q = v_dq(2);
            i_d = i_dq(1);
            i_q = i_dq(2);

            syms T_0 V_ref E_0 
            syms psi_d psi_q psi_fd psi_1d psi_gq psi_kq w theta T_m V_f T_dd T_qd E_fd R_f
            
            T_e = psi_d*i_q - psi_q*i_d;
            U_mag = sqrt((R_a*i_d+(w/wb)*psi_q)^2+(-R_a*i_q+(w/wb)*psi_d)^2);

            % The equations below come from the differential equations.
            % The theta will be governed by current constraints: eqn12 and eqn13,
            % bringing two solutions for us. But only one is small-signal stable.

            eqn1 = -wb*(-v_d - (w/wb)*psi_q - R_a*i_d) == 0; % state1
            eqn2 = -wb*(-v_q + (w/wb)*psi_d - R_a*i_q) == 0; % state2
            eqn3 = (1/T_dd)*(E_fd - psi_fd...
                                             -xd_diff*(i_d-...
                                                      xd_frac*...
                                                      (psi_1d + xdl_diff*i_d - psi_fd))) == 0; % state3
            eqn4 = (1/T_ddd)*(psi_fd - psi_1d - xdl_diff*i_d); % state4
            eqn5 = (1/T_qd)*(-psi_gq +...
                                          xq_diff*(i_q-...
                                                xq_frac*...
                                                (psi_kq + xql_diff*i_q+psi_gq))) == 0; % state5
            eqn6 = (1/T_qdd)*(-psi_gq - psi_kq - xql_diff*i_q) == 0; % state6
            eqn7 = (wb/(2*H))*(T_m - T_e - D*w) == 0; % state7
            eqn9 = (1/T_A)*(K_A*(V_ref - U_mag) - V_f + K_A*R_f - E_fd*(K_A*K_F/T_F)) == 0;
            eqn10 = (1/T_F)*(-R_f + E_fd*(K_F/T_F)) == 0;
            eqn11 = (1/T_E)*(V_f - E_fd*K_E) == 0;
            eqn12 = i_d + (1/x_ddd)*(psi_d +...
                             -psi_1d*(xdd_diff/xdl_diff) +...
                             -psi_fd*(xddl_diff/xdl_diff)) == 0;
            eqn13 = i_q + (1/x_qdd)*(psi_q -...
                             -psi_kq*(xqd_diff/xql_diff) +...
                             psi_gq*(xqql_diff/xql_diff)) == 0;
            % The theta will be governed by current constraints: eqn12 and eqn13.
            w = wb;
            [psi_d psi_q psi_fd psi_1d psi_gq psi_kq T_m R_f V_f E_fd V_ref theta]...
                      = solve(eval([eqn1, eqn2, eqn3, eqn4, eqn5, eqn6, eqn7,...
                      eqn9, eqn10, eqn11, eqn12, eqn13]),...
                      [psi_d psi_q psi_fd psi_1d psi_gq psi_kq T_m R_f V_f E_fd V_ref theta]);
            sol = 2;
            % We have two solutions here, sol = 1 is unstable.
            psi_d  = vpa(psi_d(sol), 16);  % 16 digits remain here
            psi_q  = vpa(psi_q(sol), 16);
            psi_fd = vpa(psi_fd(sol), 16);
            psi_1d = vpa(psi_1d(sol), 16);
            psi_gq = vpa(psi_gq(sol), 16);
            psi_kq = vpa(psi_kq(sol), 16);
            theta  = vpa(theta(sol), 16);
            T_m    = vpa(T_m(sol), 16);
            V_f    = vpa(V_f(sol), 16);
            R_f    = vpa(R_f(sol), 16);
            E_fd   = vpa(E_fd(sol), 16);
            V_ref  = vpa(V_ref(sol), 16);
            Px_1   = 0;
            Px_2   = 0;
            Px_3   = 0;
            T_0    = T_m;
            w_ref  = wb;

            v_d = eval(v_d);
            v_q = eval(v_q);

            % Set equilibrium
            x_e = eval([psi_d; psi_q; psi_fd; psi_1d; psi_gq; psi_kq; w; theta; T_m; V_f; R_f; E_fd; Px_1; Px_2; Px_3]);
            u_e = eval([v_d; v_q; T_0; w_ref; V_ref]);
            xi  = xi;
        end
        
    	% State space model
        function [Output] = StateSpaceEqu(obj,x,u,CallFlag)
          	
            % Get parameters in pu
            x_d   = obj.Para(1);   %
            x_dd  = obj.Para(2);   % X_dd represents Xd'
            x_ddd = obj.Para(3);   %
            x_q   = obj.Para(4);   %
            x_qd  = obj.Para(5);   % X_qd represents Xq'
            x_qdd = obj.Para(6);   %
            x_l   = obj.Para(7);   % or x_2, x_ls in some practical cases
            R_a   = obj.Para(8);   % 
            T_dd  = obj.Para(9);   % T_dd represents Td'
            T_ddd = obj.Para(10);  %
            T_qd  = obj.Para(11);  % T_qd represents Tq'
            T_qdd = obj.Para(12);  %
            H     = obj.Para(13);  % Inertial constant
            S_A   = obj.Para(14);  % Saturation parameters
            S_B   = obj.Para(15);
            S_N   = obj.Para(16);
            D     = obj.Para(17);       % Damp constant
            SpeedPara   = obj.Para(18); % simplified governor parameters
            K_A = obj.Para(19);         % simplified exciter parameters
            EnableSaturation = obj.Para(20); % Saturation setting
            wb = obj.Para(21);

            T_speed = 0.01; % Time constant of governer

            T_E = 0.314;    % Time constant of exciter
            T_A = 0.02;     % Time constant of AV
            T_F = 0.35;     % Time constant of Feedback
            K_E = 1.0;      % Exciter Gain 
            K_F = 0.063;    % Feedback Gain

            T_1 = 0.5;      % Phase 1
            T_2 = 0.05;     % Phase Compensator 1
            T_3 = 0.5;      % Phase 2
            T_4 = 0.05;     % Phase Compensator 2
            T_w = 10;       % Washout
            K_s = 2;        % PSS Gain

%           you can change the paramters in the simulation, to observe some
%           influences from the parameters.
%
%             if obj.Timer>5
%             % Get parameters in pu
%             x_d   = 0.2;   %
%             x_dd  = 0.033;   % X_dd represents Xd'
%             x_ddd = 0.023;   %
%             x_q   = 0.19;   %
%             x_qd  = 0.061;   % X_qd represents Xq'
%             x_qdd = 0.051;   %
%             x_l   = 0.022;   % or x_2, x_ls in some practical cases
%             H=0.5;
%             end

            % Calculated parameters
            xd_diff = x_d - x_dd;
            xdd_diff = x_dd - x_ddd;
            xdl_diff = x_dd-x_l;
            xddl_diff = x_ddd-x_l; 

            xq_diff = x_q - x_qd;
            xqd_diff = x_qd - x_qdd;
            xql_diff = x_qd - x_l;
            xqql_diff = x_qdd - x_l;

            xd_frac = xdd_diff/xdl_diff^2;
            xq_frac = xqd_diff/xql_diff^2;

        	% Get state
            psi_d  = x(1);
            psi_q  = x(2);
            psi_fd = x(3);
            psi_1d = x(4);
            psi_gq = x(5);
            psi_kq = x(6);
            w      = x(7);
            theta  = x(8);
            T_m    = x(9);
            V_f    = x(10);
            R_f    = x(11);
            E_fd   = x(12);
            Px_1   = x(13);
            Px_2   = x(14);
            Px_3   = x(15);
            
            % Get input
            v_d   = u(1);
            v_q   = u(2);
            T_0   = u(3);
            w_ref = u(4);
            V_ref = u(5);

            % saturation setting
            V_f_limit_H = 6;
            V_f_limit_L = -6;

%             if obj.Timer>10
%                 V_ref =V_ref+0.05;
%                 T_0 = T_0+0.05;
%                 w_ref = w_ref+10;
%             end
%             
            % Dynamic equations
          	% dx/dt = f(x,u)
            % y     = g(x,u)
            if CallFlag == 1
                % ### Call dynamic equation: dx/dt = f(x,u)
                i_d = -(1/x_ddd)*(psi_d +...
                                 -psi_1d*(xdd_diff/xdl_diff) +...
                                 -psi_fd*(xddl_diff/xdl_diff));
                i_q = -(1/x_qdd)*(psi_q -...
                                 -psi_kq*(xqd_diff/xql_diff) +...
                                 psi_gq*(xqql_diff/xql_diff));

                dpsi_d = -wb*(-v_d - (w/wb)*psi_q - R_a*i_d);  % state1
                dpsi_q = -wb*(-v_q + (w/wb)*psi_d - R_a*i_q);  % state2  
                dpsi_fd = (1/T_dd)*(E_fd - psi_fd...
                                 -xd_diff*(i_d-...
                                          xd_frac*...
                                          (psi_1d + xdl_diff*i_d - psi_fd)));   % state3
                dpsi_1d = (1/T_ddd)*(psi_fd - psi_1d - xdl_diff*i_d);           % state4
                dpsi_gq = (1/T_qd)*(-psi_gq +...
                              xq_diff*(i_q-...
                                    xq_frac*...
                                    (psi_kq + xql_diff*i_q+psi_gq)));          % state5
                dpsi_kq = (1/T_qdd)*(-psi_gq - psi_kq - xql_diff*i_q);         % state6
                
                T_e = psi_d*i_q - psi_q*i_d;
                dw = (wb/(2*H))*(T_m - T_e - D*w);      % state7
                dtheta = w;                             % state8
                
                U_mag = abs(-(R_a*i_d+(w/wb)*psi_q)+...
                        (-R_a*i_q+(w/wb)*psi_d)*1i);

                if 0    % simplifed governer and exciter with 1 order
                    % T_0 comes from the prime motor
                    dT_m  = (1/T_speed)*(T_0 + (SpeedPara*(w_ref - w)/wb) - T_m); % state9
                    dV_f  = (1/T_A)*(K_A*(V_ref - U_mag  + K_s*Px_3) - V_f);      % state10
                    dR_f =  0;
                    dE_fd = 0;
                    dPx_1 = 0;
                    dPx_2 = 0;
                    dPx_3 = 0;

                    if 1  % exciter saturation can be considered
                        V_f = min(V_f,V_f_limit_H);
                        V_f = max(V_f,V_f_limit_L);
                    end  
                end

                if 0    % other possible governer and exciter
                    dT_m  = (1/T_speed)*(T_0 + (SpeedPara*(w_ref - w)/wb) - T_m);
                    dV_f  = (1/T_A)*(K_A*(V_ref - U_mag + K_s*Px_3) - V_f + K_A*R_f - E_fd*(K_A*K_F/T_F));
                    dR_f =  (1/T_F)*(-R_f + E_fd*(K_F/T_F));
                    dE_fd = (1/T_E)*(V_f - E_fd*K_E);
                    dPx_1 = 0;
                    dPx_2 = 0;
                    dPx_3 = 0;

                    if 1  % exciter saturation can be considered
                        V_f = min(V_f,V_f_limit_H);
                        V_f = max(V_f,V_f_limit_L);
                    end 
                end

                if 1    % governer, and IEEE DC Type1 exciter with \Delta w PSS
                    dT_m  = (1/T_speed)*(T_0 + (SpeedPara*(w_ref - w)/wb) - T_m);
                    dV_f  = (1/T_A)*(K_A*(V_ref - U_mag + K_s*Px_3) - V_f + K_A*R_f - E_fd*(K_A*K_F/T_F));
                    dR_f =  (1/T_F)*(-R_f + E_fd*(K_F/T_F));
                    dE_fd = (1/T_E)*(V_f - E_fd*K_E);
                    dPx_1 = (-1/T_w)*Px_1 + (1/T_w)*(w - wb);
                    dPx_2 = (-1/T_1)*Px_2 + (1/T_1)*(Px_1 + T_2*dPx_1);
                    dPx_3 = (-1/T_3)*Px_3 + (1/T_3)*(Px_2 + T_4*dPx_2);

                    if 1  % exciter saturation can be considered
                        V_f = min(V_f,V_f_limit_H);
                        V_f = max(V_f,V_f_limit_L);
                    end 

                end

                % Flux saturation considered in the future
                if EnableSaturation
                    % pass
                end

                f_xu = [dpsi_d; dpsi_q; dpsi_fd; dpsi_1d; dpsi_gq; dpsi_kq;...
                        dw; dtheta; dT_m; dV_f; dR_f; dE_fd; dPx_1; dPx_2; dPx_3];
                Output = f_xu;

            elseif CallFlag == 2
                % ### Call output equation: y = g(x,u)
                % Just be careful that the input and output must coincide
                % with the signal sequence as well as the ports in simulink model.
                
                i_d = -(1/x_ddd)*(psi_d +...
                                 -psi_1d*(xdd_diff/xdl_diff) +...
                                 -psi_fd*(xddl_diff/xdl_diff));
                i_q = -(1/x_qdd)*(psi_q -...
                                 -psi_kq*(xqd_diff/xql_diff) +...
                                 psi_gq*(xqql_diff/xql_diff));
                T_e = psi_d*i_q - psi_q*i_d;
                U_mag = abs(-(R_a*i_d+(w/wb)*psi_q)+...
                        (-R_a*i_q+(w/wb)*psi_d)*1i);
                % We need to change the direction of the current to make it a
                % load convention.
                g_xu = [-i_d; -i_q; w; theta-pi/2; T_e; U_mag; psi_d; psi_q; psi_fd; psi_1d; psi_gq; psi_kq];
                Output = g_xu;

            end
        end
        
    end
end
