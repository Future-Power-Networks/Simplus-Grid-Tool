% This class defines the model of doubly-fed induction generator (DFIG).

% Author(s): Duange Guo

%% Notes
%
% This is just a template. For practical examples, please see
% "Inductor.m", which is a single-phase inductor, and see
% "SynchronousMachine.m", which is a three-phase synchronous machine.

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
        	State  = {'w_m',...
                      'Psi_ds' 'Psi_dr' 'Psi_qs' 'Psi_qr',...
                      'Phi_ra', 'Phi_rb', 'Phi_rc',... % RSC control
                      'igd', 'igq',...
                      'Phi_ga', 'Phi_gb', 'Phi_gc',... % GSC control
                      'V_dc',...
                      'Phipllr', 'theta_pllr',...
                      'Phipllg', 'theta_pllg'}; 	% x, state
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
            obj.Para(1);
            obj.Para(1);
            obj.Para(1);
            obj.Para(1);
            obj.Para(1);
            obj.Para(1);
            obj.Para(1);
            obj.Para(1);
            obj.Para(1);
            obj.Para(1);
            obj.Para(1);
            obj.Para(1);
            obj.Para(1);
            obj.Para(1);
            obj.Para(1);
            obj.Para(1);

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

        	% Get state
            w_m  = x(1);
            Psi_ds  = x(2);
            Psi_qs = x(3);
            Psi_dr = x(4);
            Psi_qr = x(5);
            Phi_ra = x(6);
            Phi_rc      = x(7);
            Phi_rb  = x(8);
            i_gd    = x(9);
            i_gq    = x(10);
            Phi_ga    = x(11);
            Phi_gc   = x(12);
            Phi_gb   = x(13);
            Vdc   = x(14);
            Px_3   = x(15);
            
            % Get input
            u(1);
            
            % State space equations
          	% dx/dt = f(x,u)
            % y     = g(x,u)
            if CallFlag == 1
                % ### Call state equation: dx/dt = f(x,u)

                Psi2I = inv([Lm Ls 0 0;Lr Lm 0 0;0 0 Lm Ls;0 0 Lr Lm]);
                I = Psi2I*[Psids;Psidr;Psiqs;Psiqr];
                Ird = I(1); Isd = I(2); Irq = I(3); Isq = I(4);
                Ptot = usq*Isq + usd*Isd + urq*Irq + urd*Ird; % output1
                Te = p*Lm*(Ird*Isq - Isd*Irq); % output2

                dw_m = (p/J)*(T_e-T_L); % state1
                w_r = w_g - w_m;
                dPsi_ds = u_sd - Rs*Isd + w_g*Psi_qs; % state2
                dPsi_qs = u_sq - Rs*Isq - w_g*Psi_ds; % state3
                dPsi_dr = u_rd - Rr*Ird + w_r*Psi_qr; % state4
                dPsi_qr = u_rq - Rr*Irq - w_r*Psi_dr; % state5

                alpha1 = -Lm/Ls; alpha2 = Lr + Lm*alpha1;
                PHI = sqrt(2/3)*sqrt(Psids^2 + Psiqs^2);
                
                compensatord = wr*Irq*(-alpha2);
                compensatorq = wr*PHI*alpha1 + wr*Ird*alpha2;
                
                dPhi_ra = KId*(Ird_ref - I_rd); % state6
                u_rd = KPd*(Ird_ref - Ird) + Phi_ra + compensatord;
                
                dPhi_rc = KIn*(wm_ref - w_m); % state7
                Irq_ref = (-1)*(KPn*(wm_ref - w_m) + Phi_rc);
                
                dPhi_rb = KIq*(Irq_ref - I_rq); % state8
                u_rq = KPq*(Irq_ref - I_rq) + Phi_rb + compensatorq;

                di_gd = (1/Lg)*(u_gd - Rg*i_gd - vrefd + Lg*w_g*i_gq); % state9
                di_gq = (1/Lg)*(u_gq - Rg*i_gq - vrefq - Lg*w_g*i_gd); % state10

                compensatorgd = ugd + w_g*Lg*igq;
                compensatorgq = ugq - w_g*Lg*igd;
                
                dPhi_ga = KIqg*(Igq_ref - i_gq);              % state11
                vrefq = compensatorgq - (KPqg*(Igq_ref - i_gq)+ Phi_ga);
                
                dPhi_gc = KI_V*(Vdc_ref - V_dc);              % state12
                Igd_ref = KP_V*(Vdc_ref - V_dc) + Phi_gc;
                
                dPhi_gb = KIdg*(Igd_ref - igd);               % state13
                vrefd = compensatorgd - (KPdg*(Igd_ref - i_gd)+ Phi_gb);

                Prsc = u_rd*I_rd + u_rq*Irq;
                Pgsc = u_gd*i_gd + u_gq*i_gq;
                dVdc = (1/(C_bus*V_dc))*(P_gsc-P_rsc);        % state14

                dPhiplla = ki_pll*(u_sd-usd_ref);             % state15
                dtheta_pll = Phiplla + kp_pll*(u_sd-usd_ref); % state16

                f_xu = [];
                Output = f_xu;
            elseif CallFlag == 2
                % ### Call output equation: y = g(x,u)
                % Just be careful that the input and output must coincide
                % with the signal sequence as well as the ports in simulink model.

                g_xu = [];
                Output = g_xu;
            end
        end
        
    end
end
