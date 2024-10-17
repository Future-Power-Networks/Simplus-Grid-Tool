% This class defines the model of doubly-fed induction generator (DFIG).

% Author(s): Duange Guo

%% Notes
%
% This model is in SI, not p.u.
% Asynchrnous Machine is in load convention
% and GSC is in Generator convention
% Control schemes like MPPT and Low voltage protection are neglected
% Outer PI controllers for Q are also neglected
% The rotor quantities are referred to the stator
% The system dynamic is in global steady frame, which is different from
% other apparatus.
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

        function [State,Input,Output] = SignalList(obj)
        	State  = {'w_m',...                        % Induction
                      'Psi_ds' 'Psi_dr' 'Psi_qs' 'Psi_qr',...
                      'Turbine_w', 'shaft_w'...        % Axis
                      'Phi_ra', 'Phi_rc','Phi_rb',...  % RSC control
                      'V_dc',...
                      'i_gd', 'i_gq',...               % Filter
                      'Phi_ga', 'Phi_gc', 'Phi_gb',... % GSC control
                      'Phipllr', 'theta_pllr',...      % PLL_RSC
                      'Phipllg', 'theta_pllg'};
            Input  = {'v_d', 'v_q',...                 % Global voltage
                      'wm_ref', 'Vdc_ref',...
                      'Ird_ref', 'Igq_ref'};
            Output = {'i_d', 'i_q',...                 % Global Current
                      'w_m', 'V_dc','T_e',...
                      'theta_r', 'theta_g',...
                      'u_rd', 'u_rq', 'I_rd', 'I_rq',...
                      'vrefd', 'vrefq', 'i_gd', 'i_gq'};                
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
                      
            % Base Value
            f   = obj.Para(1);
            p   = 2;
            w_g = 2*pi*f;
            Sbase = obj.Para(2);
            Vbase = obj.Para(3);
            Ibase = Sbase/(sqrt(3)*Vbase);
            wbase = w_g;
            Zbase = Vbase^2/Sbase;
            Lbase = Zbase/w_g;
            Tbase = Sbase/(w_g/p);
            
            % Reference
            Vdc_ref = obj.Para(4);
            Igq_ref = obj.Para(5);

            % Default parameters:
            % Induction Motor
            Hg  = 0.8;
            J   = (2*Hg*p^2*Sbase)/w_g^2;
            Rs  = 0.011*Zbase;
            Rr  = 1.1*Rs;
            Lm  = 4*Lbase;
            Ls  = 1.01*Lm;
            Lr  = 1.01*Lm;
            % Wind Turbine
            rho = 1.225;
            R = 44;
            Pitch = 0;
            N = 100;
            Ksh = 1.2543e+04;
            D_mutual = 2.5e3;
            H_WT = 0.4;
            wind = 12;
            % RSC
            sigma = 1- Lm^2/(Ls*Lr); 
            tau_i = (sigma*Lr)/Rr;
            tau_n = Hg/16;
            wni = 1/tau_i;
            wnn = 1/tau_n;
            KPn = (2*wnn*J)/p;
            KIn = ((wnn^2)*J)/p;
            KPq = (2*wni*sigma*Lr)-Rr;
            KIq = (wni^2)*Lr*sigma;
            KPd = (2*wni*sigma*Lr)-Rr;
            KId = (wni^2)*Lr*sigma;
            % AC Filter
            Lg  = 483e-6;
            Rg  = 20e-6;
            % DC Capacitor
            C_bus = 80e-3;
            % GSC
            K_pg = 1/(1.5*Vbase); 
            K_qg = -K_pg;
            KP_V = -1000;
            KI_V = -300000;
            KPdg = 2*(f*2*pi)*Lg - Rg;
            KIdg = Lg*(f*2*pi)^2;
            KPqg = KPdg;
            KIqg = KIdg;
            % PLL - 4
            rPLLp = obj.Para(6);
            rPLLi = obj.Para(7);
            gPLLp = obj.Para(8);
            gPLLi = obj.Para(9);

            % Calculate equilibrium
            S_D0 = P/V;
            S_Q0 = Q/V;
            S_DQ0 = S_D0 + 1j*S_Q0;
            i_abs = abs(S_DQ0);
            ui_argdiff = angle(S_DQ0);
            i_arg = xi - ui_argdiff;
            i_DQ = i_abs * exp(1i * i_arg);
            i_Q  = imag(i_DQ); % the global current should be output
            i_D = real(i_DQ);
            v_DQ = V * exp(1i * xi);
            v_d = real(v_DQ);  % the global voltage should be input
            v_q = imag(v_DQ);
            
            syms Phipllr theta_pllr Phipllg theta_pllg w_m...
                Psi_ds Psi_dr Psi_qs Psi_qr...
                Turbine_w shaft_w...
                Phi_ra Phi_rc Phi_rb...
                V_dc i_gd i_gq...
                Phi_ga Phi_gc Phi_gb...
                Ird_ref wm_ref

            theta_r = theta_pllr - pi;
            theta_g = theta_pllg;
            Tr_DQdq = [cos(theta_r) -sin(theta_r); sin(theta_r) cos(theta_r)];
            vr_dq_real = Vbase*Tr_DQdq*[v_d;v_q]; % RSC PLL
            Tg_DQdq = [cos(theta_g) -sin(theta_g); sin(theta_g) cos(theta_g)];
            vg_dq_real = Vbase*Tg_DQdq*[v_d;v_q]; % GSC PLL
            Tr_dqDQ = [cos(-theta_r) -sin(-theta_r); sin(-theta_r) cos(-theta_r)];
            Tg_dqDQ = [cos(-theta_g) -sin(-theta_g); sin(-theta_g) cos(-theta_g)];
            u_sd = vr_dq_real(1);
            u_sq = vr_dq_real(2);
            u_gd = vg_dq_real(1);
            u_gq = vg_dq_real(2);
            eqnd = u_sq/Vbase - V == 0;          % state14
            eqnq = u_gd/Vbase - V == 0;          % state16
            [pllr pllg] = solve(eval([eqnd, eqnq]),[theta_pllr theta_pllg]);
            theta_pllr = real(vpa(pllr(2)));
            theta_pllg = real(vpa(pllg(2)));
            Phipllr = wbase;                     % state15
            Phipllg = wbase;                     % state17
            
            u_sd = round(eval(u_sd));
            u_sq = round(eval(u_sq));
            u_gd = round(eval(u_gd));
            u_gq = round(eval(u_gq));
            
            % Turbine and Axis
            lambda = Turbine_w/N * R / wind;
            lambda_i = 1/((1/(lambda-0.02*Pitch)+(0.003/(Pitch^3+1))));
            Cp = 0.73*(151/lambda_i-0.58*Pitch-0.002*Pitch^2.14-13.2)*(exp(-18.4/lambda_i));
            P_wt = 0.5*rho*pi*(R)^2*(wind)^3*Cp;
            T_wt = -P_wt/(Turbine_w/N)/N;
            T_shaft = Ksh*shaft_w + D_mutual*(Turbine_w - w_m/p);
            % Induction
            Psi2I = inv([Lm Ls 0 0;Lr Lm 0 0;0 0 Lm Ls;0 0 Lr Lm]);
            I = Psi2I*[Psi_ds;Psi_dr;Psi_qs;Psi_qr];
            I_rd = I(1); I_sd = I(2); I_rq = I(3); I_sq = I(4);
            w_r = w_g - w_m;
            T_e  = 1.5*p*Lm*(I_rd*I_sq - I_sd*I_rq);
            % RSC
            alpha1 = Lm/Ls; alpha2 = Lr - Lm*alpha1;
            PHI = sqrt(Psi_ds^2 + Psi_qs^2);
            compensatord = w_r*I_rq*(-alpha2);
            compensatorq = w_r*PHI*alpha1 + w_r*I_rd*alpha2;
            u_rd = (KPd*(Ird_ref - I_rd) + Phi_ra + compensatord);
            Irq_ref = (KPn*(wm_ref - w_m) + Phi_rc)/(-1.5*p*alpha1*PHI);
            u_rq = (KPq*(Irq_ref - I_rq) + Phi_rb + compensatorq);
            % GSC
            compensatorgd = - w_g*Lg*i_gq;
            compensatorgq = + w_g*Lg*i_gd;
            vrefq = (KPqg*(Igq_ref - i_gq)+ Phi_ga) + compensatorgq;
            Igd_ref = (KP_V*(Vdc_ref - V_dc) + Phi_gc)*K_pg;
            vrefd = (KPdg*(Igd_ref - i_gd)+ Phi_gb) - compensatorgd;
            % Output
            P_rsc = 1.5*(u_rd*I_rd + u_rq*I_rq);
            P_gsc = 1.5*(u_gd*i_gd + u_gq*i_gq);
            is_dq_pu = (1/Ibase)*Tr_dqDQ*[I_sd;I_sq];
            ig_dq_pu = (1/Ibase)*Tg_dqDQ*[i_gd;i_gq];
            i_d = is_dq_pu(1)-ig_dq_pu(1);
            i_q = is_dq_pu(2)-ig_dq_pu(2);

            eqn1 = (p/J)*(T_e - T_shaft) == 0; % state1
            eqn2 = (1/(2*H_WT))*(T_wt - T_shaft) == 0;% state20
            eqn3 = Turbine_w - w_m/p;  % state21
            eqn4 = u_sd - Rs*I_sd + w_g*Psi_qs == 0; % state2
            eqn5 = u_sq - Rs*I_sq - w_g*Psi_ds == 0; % state3
            eqn6 = u_rd - Rr*I_rd + w_r*Psi_qr == 0; % state4
            eqn7 = u_rq - Rr*I_rq - w_r*Psi_dr == 0; % state5
            eqn8 = KId*(Ird_ref - I_rd) == 0;
            eqn9 = KIn*(wm_ref - w_m) == 0;
            eqn10 = KIq*(Irq_ref - I_rq) == 0;
            eqn11 = (1/(C_bus*V_dc))*(-P_gsc-P_rsc) == 0;
            eqn12 = (1/Lg)*(vrefd - u_gd - Rg*i_gd + Lg*w_g*i_gq) == 0;
            eqn13 = (1/Lg)*(vrefq - u_gq - Rg*i_gq - Lg*w_g*i_gd) == 0;
            eqn14 = KIqg*(Igq_ref - i_gq) == 0;
            eqn15 = KI_V*(Vdc_ref - V_dc) == 0;
            eqn16 = KIdg*(Igd_ref - i_gd) == 0;
            eqn17 = i_d == i_D;
            eqn18 = i_q == i_Q;

            [Turbine_w w_m Psi_ds Psi_qs Psi_dr Psi_qr shaft_w...
                Phi_ra Phi_rc Phi_rb...
                V_dc i_gd i_gq Phi_ga Phi_gc Phi_gb wm_ref Ird_ref]...
                = vpasolve(eval([eqn1, eqn2, eqn3, eqn4, eqn5, eqn6, eqn7,...
                                 eqn8, eqn9, eqn10, eqn11, eqn12, eqn13, eqn14,...
                                 eqn15, eqn16, eqn17, eqn18]),...
                                [Turbine_w w_m Psi_ds Psi_qs Psi_dr Psi_qr shaft_w...
                                 Phi_ra Phi_rc Phi_rb...
                                 V_dc i_gd i_gq Phi_ga Phi_gc Phi_gb wm_ref Ird_ref],...
                                [40,300;-2*w_g,2*w_g;-inf,inf;-inf,inf;-inf,inf;-inf,inf;-inf,inf;...
                                 -inf,inf;-inf,inf;-inf,inf;-inf,inf;-inf,inf;-inf,inf;...
                                 -inf,inf;-inf,inf;-inf,inf;-2*w_g,2*w_g;-inf,inf]);

            % Set equilibrium
            x_e = eval([w_m;...
                       Psi_ds; Psi_qs; Psi_dr; Psi_qr;...
                       Turbine_w; shaft_w;...
                       Phi_ra; Phi_rc; Phi_rb;...
                       V_dc;...
                       i_gd; i_gq;...
                       Phi_ga; Phi_gc; Phi_gb;...
                       Phipllr; theta_pllr;...
                       Phipllg; theta_pllg]);
            u_e = eval([v_d; v_q; wm_ref; Vdc_ref; Ird_ref; Igq_ref]);
            xi  = xi;
        end
        
    	% State space model
        function [Output] = StateSpaceEqu(obj,x,u,CallFlag)
          	% Get parameter
            % Base Value
            f   = obj.Para(1);
            p   = 2;
            w_g = 2*pi*f;
            Sbase = obj.Para(2);
            Vbase = obj.Para(3);
            Ibase = Sbase/(sqrt(3)*Vbase);
            wbase = w_g;
            Zbase = Vbase^2/Sbase;
            Lbase = Zbase/w_g;
            Tbase = Sbase/(w_g/p);
            
            % Reference
            Vdc_ref = 1200;
            Igq_ref = obj.Para(3);
            usd_ref = 0;
            ugq_ref = 0;

            % Default parameters:
            % Induction Motor
            Hg  = 0.8;
            J   = (2*Hg*p^2*Sbase)/w_g^2;
            Rs  = 0.011*Zbase;
            Rr  = 1.1*Rs;
            Lm  = 4*Lbase;
            Ls  = 1.01*Lm;
            Lr  = 1.01*Lm;
            % Wind Turbine
            rho = 1.225;
            R = 44;
            Pitch = 0;
            N = 100;
            Ksh = 1.2543e+04;
            D_mutual = 2.5e3;
            H_WT = 0.4;
            wind = 12;
            % RSC
            sigma = 1- Lm^2/(Ls*Lr); 
            tau_i = (sigma*Lr)/Rr;
            tau_n = Hg/16;
            wni = 1/tau_i;
            wnn = 1/tau_n;
            KPn = (2*wnn*J)/p;
            KIn = ((wnn^2)*J)/p;
            KPq = (2*wni*sigma*Lr)-Rr;
            KIq = (wni^2)*Lr*sigma;
            KPd = (2*wni*sigma*Lr)-Rr;
            KId = (wni^2)*Lr*sigma;
            % AC Filter
            Lg  = 483e-6;
            Rg  = 20e-6;
            % DC Capacitor
            C_bus = 80e-3;
            % GSC
            K_pg = 1/(1.5*Vbase); 
            K_qg = -K_pg;
            KP_V = -1000;
            KI_V = -300000;
            KPdg = 2*(f*2*pi)*Lg - Rg;
            KIdg = Lg*(f*2*pi)^2;
            KPqg = KPdg;
            KIqg = KIdg;
            % PLL - 4
            rPLLp = obj.Para(6);
            rPLLi = obj.Para(7);
            gPLLp = obj.Para(8);
            gPLLi = obj.Para(9);

            % Saturation and limit
            pitch_max  = 27;
            pitch_rate = 10;

        	% Get state
            w_m        = x(1);
            Psi_ds     = x(2);
            Psi_qs     = x(3);
            Psi_dr     = x(4);
            Psi_qr     = x(5);
            Turbine_w  = x(6);
            shaft_w    = x(7);
            Phi_ra     = x(8);
            Phi_rc     = x(9);
            Phi_rb     = x(10);
            V_dc       = x(11);
            i_gd       = x(12);
            i_gq       = x(13);
            Phi_ga     = x(14);
            Phi_gc     = x(15);
            Phi_gb     = x(16);
            Phipllr    = x(17);
            theta_pllr = x(18);
            Phipllg    = x(19);
            theta_pllg = x(20);

            % Get input
            v_d = u(1); 
            v_q = u(2); 
            wm_ref = u(3); 
            Vdc_ref = u(4); 
            Ird_ref = u(5); 
            Igq_ref = u(6);

            % State space equations
          	% dx/dt = f(x,u)
            % y     = g(x,u)
            if CallFlag == 1
                % ### Call state equation: dx/dt = f(x,u)
                % PLL
                theta_r = theta_pllr - pi;
                theta_g = theta_pllg;
                Tr_DQdq = [cos(theta_r) -sin(theta_r); sin(theta_r) cos(theta_r)];
                vr_dq_real = Vbase*Tr_DQdq*[v_d;v_q]; % RSC PLL
                Tg_DQdq = [cos(theta_g) -sin(theta_g); sin(theta_g) cos(theta_g)];
                vg_dq_real = Vbase*Tg_DQdq*[v_d;v_q]; % GSC PLL
                u_sd = vr_dq_real(1);
                u_sq = vr_dq_real(2);
                u_gd = vg_dq_real(1);
                u_gq = -vg_dq_real(2);
                
                Pitch = 0;
                lambda = Turbine_w/N * R / wind;
                lambda_i = 1/((1/(lambda-0.02*Pitch)+(0.003/(Pitch^3+1))));
                Cp = 0.73*(151/lambda_i-0.58*Pitch-0.002*Pitch^2.14-13.2)*(exp(-18.4/lambda_i));
                P_wt = 0.5*rho*pi*(R)^2*(wind)^3*Cp;
                T_wt = -P_wt/(Turbine_w/N)/N;
                T_shaft = Ksh*shaft_w + D_mutual*(Turbine_w - w_m/p);
                
                Psi2I = inv([Lm Ls 0 0;Lr Lm 0 0;0 0 Lm Ls;0 0 Lr Lm]);
                I = Psi2I*[Psi_ds;Psi_dr;Psi_qs;Psi_qr];
                I_rd = I(1); I_sd = I(2); I_rq = I(3); I_sq = I(4);
                w_r = w_g - w_m;
                T_e  = 1.5*p*Lm*(I_rd*I_sq - I_sd*I_rq);
                
                alpha1 = Lm/Ls; alpha2 = Lr - Lm*alpha1;
                PHI = sqrt(Psi_ds^2 + Psi_qs^2);
                compensatord = w_r*I_rq*(-alpha2);
                compensatorq = w_r*PHI*alpha1 + w_r*I_rd*alpha2;
                u_rd = (KPd*(Ird_ref - I_rd) + Phi_ra + compensatord);
                Irq_ref = (KPn*(wm_ref - w_m) + Phi_rc)/(-1.5*p*alpha1*PHI);
                u_rq = (KPq*(Irq_ref - I_rq) + Phi_rb + compensatorq);
                
                compensatorgd = - w_g*Lg*i_gq;
                compensatorgq = + w_g*Lg*i_gd;
                vrefq = (KPqg*(Igq_ref - i_gq)+ Phi_ga) + compensatorgq;
                Igd_ref = (KP_V*(Vdc_ref - V_dc) + Phi_gc)*K_pg;
                vrefd = (KPdg*(Igd_ref - i_gd)+ Phi_gb) - compensatorgd;
                
                P_rsc = 1.5*(u_rd*I_rd + u_rq*I_rq);
                P_gsc = 1.5*(u_gd*i_gd + u_gq*i_gq);
                
                dw_m        = (p/J)*(T_e - T_shaft);
                dPsi_ds     = u_sd - Rs*I_sd + w_g*Psi_qs;
                dPsi_qs     = u_sq - Rs*I_sq - w_g*Psi_ds;
                dPsi_dr     = u_rd - Rr*I_rd + w_r*Psi_qr;
                dPsi_qr     = u_rq - Rr*I_rq - w_r*Psi_dr;
                dTurbine_w  = (1/(2*H_WT))*(T_wt - T_shaft); 
                dshaft_w    = Turbine_w - w_m/p;
                dPhi_ra     = KId*(Ird_ref - I_rd); 
                dPhi_rc     = KIn*(wm_ref - w_m); 
                dPhi_rb     = KIq*(Irq_ref - I_rq);
                dVdc        = (1/(C_bus*V_dc))*(-P_gsc-P_rsc);
                di_gd       = (1/Lg)*(vrefd - u_gd - Rg*i_gd + Lg*w_g*i_gq); 
                di_gq       = (1/Lg)*(vrefq - u_gq - Rg*i_gq - Lg*w_g*i_gd);
                dPhi_ga     = KIqg*(Igq_ref - i_gq); 
                dPhi_gc     = KI_V*(Vdc_ref - V_dc); 
                dPhi_gb     = KIdg*(Igd_ref - i_gd);
                dPhipllr    = rPLLi*(u_sd/Vbase - usd_ref);
                dtheta_pllr = Phipllr + rPLLp*(u_sd/Vbase-usd_ref);
                dPhipllg    = gPLLi*(u_gq/Vbase - ugq_ref); 
                dtheta_pllg = Phipllg + gPLLp*(u_gq/Vbase - ugq_ref);

                f_xu = [dw_m; ...
                        dPsi_ds; dPsi_qs; dPsi_dr; dPsi_qr;...
                        dTurbine_w; dshaft_w;...
                        dPhi_ra; dPhi_rc; dPhi_rb;...
                        dVdc;...
                        di_gd; di_gq;...
                        dPhi_ga; dPhi_gc; dPhi_gb;...
                        dPhipllr; dtheta_pllr;...
                        dPhipllg; dtheta_pllg];
                Output = f_xu;

            elseif CallFlag == 2
                % ### Call output equation: y = g(x,u)
                % Just be careful that the input and output must coincide
                % with the signal sequence as well as the ports in simulink model.
                
                Psi2I = inv([Lm Ls 0 0;Lr Lm 0 0;0 0 Lm Ls;0 0 Lr Lm]);
                I = Psi2I*[Psi_ds;Psi_dr;Psi_qs;Psi_qr];
                I_rd = I(1); I_sd = I(2); I_rq = I(3); I_sq = I(4);
                w_r = w_g - w_m;

                T_e  = 1.5*p*Lm*(I_rd*I_sq - I_sd*I_rq);

                alpha1 = Lm/Ls; alpha2 = Lr - Lm*alpha1;
                PHI = sqrt(Psi_ds^2 + Psi_qs^2);
                compensatord = w_r*I_rq*(-alpha2);
                compensatorq = w_r*PHI*alpha1 + w_r*I_rd*alpha2;
                u_rd = (KPd*(Ird_ref - I_rd) + Phi_ra + compensatord);
                Irq_ref = (KPn*(wm_ref - w_m) + Phi_rc)/(-1.5*p*alpha1*PHI);
                u_rq = (KPq*(Irq_ref - I_rq) + Phi_rb + compensatorq);
                
                compensatorgd = - w_g*Lg*i_gq;
                compensatorgq = + w_g*Lg*i_gd;
                vrefq = (KPqg*(Igq_ref - i_gq)+ Phi_ga) + compensatorgq;
                Igd_ref = (KP_V*(Vdc_ref - V_dc) + Phi_gc)*K_pg;
                vrefd = (KPdg*(Igd_ref - i_gd)+ Phi_gb) - compensatorgd;
                
                theta_r = theta_pllr - pi;
                theta_g = theta_pllg;

                Tr_dqDQ = [cos(-theta_r) -sin(-theta_r); sin(-theta_r) cos(-theta_r)];
                Tg_dqDQ = [cos(-theta_g) -sin(-theta_g); sin(-theta_g) cos(-theta_g)];

                is_dq_pu = (1/Ibase)*Tr_dqDQ*[I_sd;I_sq];
                ig_dq_pu = (1/Ibase)*Tg_dqDQ*[i_gd;i_gq];
                
                i_d = is_dq_pu(1) - ig_dq_pu(1);
                i_q = is_dq_pu(2) - ig_dq_pu(2);

                g_xu = [i_d; i_q; w_m; V_dc; T_e; theta_r; theta_g;...
                        u_rd; u_rq; I_rd; I_rq; vrefd; vrefq; i_gd; i_gq];
                Output = g_xu;
            end
        end
        
    end
end
