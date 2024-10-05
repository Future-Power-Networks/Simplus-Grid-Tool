% This class defines the model of doubly-fed induction generator (DFIG).

% Author(s): Duange Guo

%% Notes
%
% This model is in SI, not p.u.
% In load convention
% The wind d
% Some of outer PI controller is omitted
%
%
%
%
%

%% Class

classdef DoublyFedInductionGenerator < SimplusGT.Class.ModelAdvance
    
  	methods
        % Constructor
        function obj = DoublyFedInductionGenerator(varargin)
            setProperties(obj,nargin,varargin{:});
        end
    end
    
    methods(Static)
        
        % Set the strings of state, input, output
        %
        % These strings are mainly for printing output and searching the
        % corresponding ports.
        %
        % For ac apparatuses, the first two inputs and outputs should be
        % {'v_d','v_q'} and {'i_d','i_q'}, and {'w'} should be output. For
        % dc apparatuses, the first input and output should be {'v'} and {'i'}.
        % No specific requirementFor other inputs, outputs, and states.
        %
        % The dimensions of x, u, y must be coinsistent in following three
        % functions: SignalList, Equilibrium, and StateSpaceEqu.
        function [State,Input,Output] = SignalList(obj)
        	State  = {'Psi_ds' 'Psi_dr' 'Psi_qs' 'Psi_qr',...
                      'Phi_ra', 'Phi_rb', 'Phi_rc',... % RSC control
                      'igd', 'igq',...
                      'Phi_ga', 'Phi_gb', 'Phi_gc',... % GSC control
                      'V_dc',...
                      'Phipllr', 'theta_pllr',...
                      'Phipllg', 'theta_pllg',...
                      'w_m', 'Turbine_w', 'shaft_w'}; 	% x, state
            Input  = {'u_d', 'u_q', 'wm_ref', 'Vdc_ref'};   % u, input
            Output = {'i_d', 'i_q', 'w_m', 'V_dc', 'T_e'};  % y, output
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
            
            % Get parameters

            % Base Value
            Ps = ;
            Vs = ;
            Is = ;
            

            % Induction Motor - 9
            J   = obj.Para(1);
            p   = obj.Para(1);
            Rs  = obj.Para(1);
            Rr  = obj.Para(1);
            Lm  = obj.Para(1);
            Ls  = obj.Para(1);
            Lr  = obj.Para(1);
            Lsi = obj.Para(1);
            Lri = obj.Para(1);
            % AC Filter
            Lg  = obj.Para(1);
            Rg  = obj.Para(1);
            % DC Capacitor
            Cbus = obj.Para(1);
            % PLL - 4
            rPLLp = obj.Para(1);
            rPLLi = obj.Para(1);
            gPLLp = obj.Para(1);
            gPLLi = obj.Para(1);
            % RSC - 6
            KPn = obj.Para(1);
            KIn = obj.Para(1);
            KPd = obj.Para(1);
            KId = obj.Para(1);
            KPq = obj.Para(1);
            KIq = obj.Para(1);
            % GSC - 6
            KP_V = obj.Para(1);
            KI_V = obj.Para(1);
            KPdg = obj.Para(1);
            KIdg = obj.Para(1);
            KPqg = obj.Para(1);
            KIqg = obj.Para(1);
            % Wind Turbine
            Ksh = ;
            D_mutual = ;
            H_WT = ;




            % Calculated parameters
            
            
            % Calculate equilibrium
            


            % Set equilibrium
            x_e = [];
            u_e = [];
            xi  = xi;
        end
        
    	% State space model
        % This function defines the state space model of this apparatus,
        % and is the core part for capturing the dynamics of this apparatus.
        %
        % This function will be called at each step, i.e., Ts, during the
        % whole precedure of the discrete simulation.
        %
        % The state space model should be a large-signal model rather than
        % a small-signal model. The linearized model will be calculated by
        % functions in the parent class and the linearization point (i.e.
        % equilibrium) is calculated above.
        function [Output] = StateSpaceEqu(obj,x,u,CallFlag)
          	% Get parameter
            obj.Para(1);
            

            % default parameters
            

            KPdg = 2*wnig*Lg - Rg;
            KIdg = Lg*wnig^2;
            KPqg = KPdg;
            KIqg = KIdg;
            KP_V = 10*0.01;
            KI_V = 110*0.0001;

            % Saturation and limit
            pitch_max  = 27;
            pitch_rate = 10;

        	% Get state
            Psi_ds     = x(1);
            Psi_qs     = x(2);
            Psi_dr     = x(3);
            Psi_qr     = x(4);
            Phi_ra     = x(5);
            Phi_rc     = x(6);
            Phi_rb     = x(7);
            i_gd       = x(8);
            i_gq       = x(9);
            Phi_ga     = x(10);
            Phi_gc     = x(11);
            Phi_gb     = x(12);
            Vdc        = x(13);
            Phipllr    = x(14);
            theta_pllr = x(15);
            Phipllg    = x(16);
            theta_pllg = x(17);
            w_m        = x(18);
            Turbine_w  = x(19);
            shaft_w    = x(20);


            
            % Get input
            u_d = u(1);
            u_q = u(2);

            
            % State space equations
          	% dx/dt = f(x,u)
            % y     = g(x,u)
            if CallFlag == 1
                % ### Call state equation: dx/dt = f(x,u)

                DQdq_Rotation = [cos(u(3)) -sin(u(3)); sin(u(3)) cos(u(3))];

                Psi2I = inv([Lm Ls 0 0;Lr Lm 0 0;0 0 Lm Ls;0 0 Lr Lm]);
                I = Psi2I*[Psids;Psidr;Psiqs;Psiqr];
                Ird = I(1); Isd = I(2); Irq = I(3); Isq = I(4);
                Ptot = usq*Isq + usd*Isd + urq*Irq + urd*Ird; % output1
                Te = p*Lm*(Ird*Isq - Isd*Irq);                % output2

                % Induction
                w_r = w_g - w_m;
                dPsi_ds = u_sd - Rs*Isd + w_g*Psi_qs; % state1
                dPsi_qs = u_sq - Rs*Isq - w_g*Psi_ds; % state2
                dPsi_dr = u_rd - Rr*Ird + w_r*Psi_qr; % state3
                dPsi_qr = u_rq - Rr*Irq - w_r*Psi_dr; % state4

                alpha1 = -Lm/Ls; alpha2 = Lr + Lm*alpha1;
                PHI = sqrt(2/3)*sqrt(Psids^2 + Psiqs^2);
                
                compensatord = wr*Irq*(-alpha2);
                compensatorq = wr*PHI*alpha1 + wr*Ird*alpha2;

                % RSC
                dPhi_ra = KId*(Ird_ref - I_rd); % state5
                u_rd = KPd*(Ird_ref - Ird) + Phi_ra + compensatord;
                
                dPhi_rc = KIn*(wm_ref - w_m); % state6
                Irq_ref = (-1)*(KPn*(wm_ref - w_m) + Phi_rc);
                
                dPhi_rb = KIq*(Irq_ref - I_rq); % state7
                u_rq = KPq*(Irq_ref - I_rq) + Phi_rb + compensatorq;

                % AC Filter
                di_gd = (1/Lg)*(u_gd - Rg*i_gd - vrefd + Lg*w_g*i_gq); % state8
                di_gq = (1/Lg)*(u_gq - Rg*i_gq - vrefq - Lg*w_g*i_gd); % state9

                % GSC
                compensatorgd = ugd + w_g*Lg*igq;
                compensatorgq = ugq - w_g*Lg*igd;
                
                dPhi_ga = KIqg*(Igq_ref - i_gq);                  % state10
                vrefq = compensatorgq - (KPqg*(Igq_ref - i_gq)+ Phi_ga);
                
                dPhi_gc = KI_V*(Vdc_ref - V_dc);                  % state11
                Igd_ref = KP_V*(Vdc_ref - V_dc) + Phi_gc;
                
                dPhi_gb = KIdg*(Igd_ref - igd);                   % state12
                vrefd = compensatorgd - (KPdg*(Igd_ref - i_gd)+ Phi_gb);

                % DC Capacitor
                Prsc = u_rd*I_rd + u_rq*Irq;
                Pgsc = u_gd*i_gd + u_gq*i_gq;
                dVdc = (1/(C_bus*V_dc))*(P_gsc-P_rsc);            % state13

                % PLL
                dPhipllr = ki_pll*(u_sd - usd_ref);               % state14
                dtheta_pllr = Phipllr + kp_pll*(u_sd-usd_ref);    % state15

                dPhipllg = gPLLi*(u_gq - ugq_ref)                 % state16
                dtheta_pllg =  Phipllg + gPLLp*(u_gq - ugq_ref);  % state17

                % Wind Turbine
                T_shaft = Ksh*shaft_w + D_mutual*(Turbine_w - (1/Dw_g)*w_m);
                T_wt = ()/();
                Pitch = 200*((1/Dw_g)*w_m-1.2);
                if 1
                    pass
                end
                dw_m = (p/J)*(T_e-T_L);                           % state18
                dTurbine_w = (1/(2*H_WT))*(T_wt - T_shaft);       % state19 
                dshaft_w = wbase*(Turbine_w - (1/Dw_g)*w_m);      % state20
                

                f_xu = [dPsi_ds; dPsi_qs; dPsi_dr; dPsi_qr; dPhi_ra; dPhi_rc; dPhi_rb; ...
                        di_gd; di_gq; dPhi_ga; dPhi_gc; dPhi_gb; dVdc; dPhipllr; dtheta_pllr; ...
                        dPhipllg; dtheta_pllg; dw_m; dTurbine_w; dshaft_w];
                Output = f_xu;

            elseif CallFlag == 2
                % ### Call output equation: y = g(x,u)
                % Just be careful that the input and output must coincide
                % with the signal sequence as well as the ports in simulink model.

                g_xu = [i_d, i_q];
                Output = g_xu;
            end
        end
        
    end
end
