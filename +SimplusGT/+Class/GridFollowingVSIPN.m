% This class defines a grid-following VSI with positive- and
% negative-sequence control.
%
% Author(s): Wanjun Gao

%% Notes
%
% The model uses load convention and admittance form. DDSRF states are
% represented in the positive frame:
%   dx_p = wf*(u - x_n - x_p)
%   dx_n = wf*(u - x_p - x_n) - 2*w*J*x_n
% Here x_n is the negative-sequence filter state mapped into the positive
% controller frame. With J = [0,-1;1,0] and dq rotation R, for
% z_n = R(-2*theta)*z_n_minus,
%   dz_n/dt = R(-2*theta)*dz_n_minus/dt - 2*w*J*z_n.

%% Class

classdef GridFollowingVSIPN < SimplusGT.Class.ModelAdvance

    % Cached equilibrium values
    properties(Access = protected)
        Vod0_ = 0; Voq0_ = 0; Iod0_ = 0; Ioq0_ = 0;
        Ild0_ = 0; Ilq0_ = 0; P0_ = 0; Q0_ = 0;
    end

    methods
        % Constructor
        function obj = GridFollowingVSIPN(varargin)
            % The MATLAB System block applies ApparatusType after construction.
            setProperties(obj,nargin,varargin{:});
        end
    end

    methods(Static)

        % Define state, input, and output signals
        function [State,Input,Output] = SignalList(~)
            State = {'i_ld','i_lq','i_ld_p_i','i_lq_p_i', ...
                     'i_ld_n_i','i_lq_n_i','v_cd','v_cq','i_od','i_oq', ...
                     'v_od_p_f','v_oq_p_f','v_od_n_f','v_oq_n_f', ...
                     'i_ld_p_f','i_lq_p_f','i_ld_n_f','i_lq_n_f', ...
                     'i_od_p_f','i_oq_p_f','i_od_n_f','i_oq_n_f', ...
                     'P_p_i','Q_p_i','w_pll_i', ...
                     'x_del_d','x_del_q','w','theta'};
            % P0 keeps the shared Simulink wrapper input nonempty.
            Input = {'v_d','v_q','P0'};
            Output = {'i_d','i_q','w','theta'};
        end

        % Calculate the equilibrium
        function [x_e,u_e,xi] = Equilibrium(obj)
            P = obj.PowerFlow(1);
            Q = obj.PowerFlow(2);
            V = obj.PowerFlow(3);
            xi = obj.PowerFlow(4);
            w = obj.PowerFlow(5);

            p = obj.GetPara(obj);
            % Align capacitor-terminal voltage with the PLL d axis.
            v_gd = V;
            v_oq = 0;
            i_od = P/V;
            i_oq = -Q/V;
            v_gq = p.Rc*i_oq+w*p.Lc*i_od;

            v_od = v_gd-p.Rc*i_od+w*p.Lc*i_oq;
            v_odq = v_od+1i*v_oq;
            i_cdq = v_odq/(p.Rd+1/(1i*w*p.Cf));
            i_ldq = (i_od+1i*i_oq)-i_cdq;
            i_ld = real(i_ldq);
            i_lq = imag(i_ldq);
            v_cdq = v_odq-p.Rd*i_cdq;
            e_dq_ss = v_odq-i_ldq*(p.Rf+1i*w*p.Lf);
            e_d_ss = real(e_dq_ss);
            e_q_ss = imag(e_dq_ss);

            obj.Vod0_ = v_od; obj.Voq0_ = v_oq;
            obj.Iod0_ = i_od; obj.Ioq0_ = i_oq;
            obj.Ild0_ = i_ld; obj.Ilq0_ = i_lq;
            obj.P0_ = v_od*i_od+v_oq*i_oq;
            obj.Q0_ = v_oq*i_od-v_od*i_oq;

            x_e = zeros(29,1);
            x_e(1:2) = [i_ld;i_lq];
            x_e(3:4) = [e_d_ss;e_q_ss];
            x_e(7:8) = [real(v_cdq);imag(v_cdq)];
            x_e(9:10) = [i_od;i_oq];

            % At balanced equilibrium, all negative-sequence states are zero.
            x_e(11:12) = [v_od;v_oq];
            x_e(15:16) = [i_ld;i_lq];
            x_e(19:20) = [i_od;i_oq];

            x_e(23:24) = [i_ld;i_lq];
            x_e(25) = w;
            x_e(26:27) = [e_d_ss;e_q_ss]*p.Td/2;
            x_e(28) = w;
            x_e(29) = xi;

            u_e = [v_gd;v_gq;0];
        end

        % State-space model
        function Output = StateSpaceEqu(obj,x,u,CallFlag)
            p = obj.GetPara(obj);

            v_gd = u(1); v_gq = u(2);
            i_ld = x(1); i_lq = x(2);
            i_ld_p_i = x(3); i_lq_p_i = x(4);
            i_ld_n_i = x(5); i_lq_n_i = x(6);
            v_cd = x(7); v_cq = x(8);
            i_od = x(9); i_oq = x(10);

            v_od = v_cd+p.Rd*(i_od-i_ld);
            v_oq = v_cq+p.Rd*(i_oq-i_lq);

            % DDSRF filter states for terminal voltage, converter-side
            % current, and grid-side current. Both sequence channels are
            % stored in the positive controller frame. Each signal uses two
            % states per sequence, giving 12 DDSRF states in total.
            x_v_p = x(11:12); x_v_n = x(13:14);
            x_il_p = x(15:16); x_il_n = x(17:18);
            x_io_p = x(19:20); x_io_n = x(21:22);
            P_p_i = x(23); Q_p_i = x(24);
            w_pll_i = x(25);
            x_del_d = x(26); x_del_q = x(27);
            w = x(28); theta = x(29);

            if CallFlag == 1
                % State equations: dx/dt = f(x,u)
                [dx_v_p,dx_v_n,v_p,~] = obj.DDSRF( ...
                    obj,[v_od;v_oq],x_v_p,x_v_n,w,p.w_pll);
                [dx_il_p,dx_il_n,il_p,il_n] = obj.DDSRF( ...
                    obj,[i_ld;i_lq],x_il_p,x_il_n,w,p.w_pll);
                [dx_io_p,dx_io_n,io_p,~] = obj.DDSRF( ...
                    obj,[i_od;i_oq],x_io_p,x_io_n,w,p.w_pll);

                v_od_p_f = v_p(1); v_oq_p_f = v_p(2);
                i_ld_p_f = il_p(1); i_lq_p_f = il_p(2);
                i_ld_n_f = il_n(1); i_lq_n_f = il_n(2);
                i_od_p_f = io_p(1); i_oq_p_f = io_p(2);

                % Power PI loops.
                P_p_cal = v_od_p_f*i_od_p_f+v_oq_p_f*i_oq_p_f;
                Q_p_cal = v_oq_p_f*i_od_p_f-v_od_p_f*i_oq_p_f;
                i_ld_r_p = P_p_i+(p.P0-P_p_cal)*p.kp;
                i_lq_r_p = Q_p_i-(p.Q0-Q_p_cal)*p.kp;
                dP_p_i = (p.P0-P_p_cal)*p.ki;
                dQ_p_i = -(p.Q0-Q_p_cal)*p.ki;

                % q-axis PLL and filtered frequency state.
                e_ang = v_oq_p_f;
                w_raw = w_pll_i+e_ang*p.kp_pll;
                dw_pll_i = e_ang*p.ki_pll;
                dw = (w_raw-w)*p.w_pll;
                dtheta = w;

                % Positive-sequence current controller.
                e_d_p = i_ld_p_i+(i_ld_p_f-i_ld_r_p)*p.kpi_p;
                e_q_p = i_lq_p_i+(i_lq_p_f-i_lq_r_p)*p.kpi_p;
                di_ld_p_i = (i_ld_p_f-i_ld_r_p)*p.kii_p;
                di_lq_p_i = (i_lq_p_f-i_lq_r_p)*p.kii_p;

                % Map the negative-sequence PI integrator into the positive
                % frame. Differentiating R(-2*theta)*z_n_minus produces the
                % additional -2*w*J*z_n term.
                J = [0,-1;1,0];
                e_n_pos = [i_ld_n_i;i_lq_n_i]+ ...
                    p.kpi_n*[i_ld_n_f;i_lq_n_f];
                di_n_i = p.kii_n*[i_ld_n_f;i_lq_n_f] ...
                    -2*w*J*[i_ld_n_i;i_lq_n_i];

                e_d_ref = e_d_p+e_n_pos(1);
                e_q_ref = e_q_p+e_n_pos(2);
                dx_del_d = -(2/p.Td)*x_del_d+e_d_ref;
                dx_del_q = -(2/p.Td)*x_del_q+e_q_ref;
                e_d = (4/p.Td)*x_del_d-e_d_ref;
                e_q = (4/p.Td)*x_del_q-e_q_ref;

                % LCL plant.
                di_ld = (v_od-e_d-p.Rf*i_ld+w*p.Lf*i_lq)/p.Lf;
                di_lq = (v_oq-e_q-p.Rf*i_lq-w*p.Lf*i_ld)/p.Lf;
                dv_cd = ((i_od-i_ld)+w*p.Cf*v_cq)/p.Cf;
                dv_cq = ((i_oq-i_lq)-w*p.Cf*v_cd)/p.Cf;
                di_od = (v_gd-v_od-p.Rc*i_od+w*p.Lc*i_oq)/p.Lc;
                di_oq = (v_gq-v_oq-p.Rc*i_oq-w*p.Lc*i_od)/p.Lc;

                Output = [di_ld;di_lq;di_ld_p_i;di_lq_p_i; ...
                          di_n_i;dv_cd;dv_cq;di_od;di_oq; ...
                          dx_v_p;dx_v_n;dx_il_p;dx_il_n;dx_io_p;dx_io_n; ...
                          dP_p_i;dQ_p_i;dw_pll_i; ...
                          dx_del_d;dx_del_q;dw;dtheta];
            elseif CallFlag == 2
                % Output equations: y = g(x,u)
                Output = [i_od;i_oq;w;theta];
            else
                error('GridFollowingVSIPN:StateSpaceEqu:InvalidCallFlag', ...
                    'CallFlag must be 1 (f) or 2 (g).');
            end
        end

        % Convert the parameter vector to physical quantities
        function p = GetPara(obj)
            q = obj.Para;
            w0 = q(16);

            % Build all fields at once to satisfy MATLAB Coder.
            p = struct( ...
                'w0',       w0, ...
                'Td',       1.5*obj.Ts/2, ...
                'Lf',       q(1)/w0, ...
                'Rf',       q(2), ...
                'Cf',       q(3)/w0, ...
                'Rd',       q(4), ...
                'Lc',       q(5)/w0, ...
                'Rc',       q(6), ...
                'w_pll',    q(7)*2*pi, ...
                'kp',       q(8), ...
                'ki',       q(9), ...
                'kp_pll',   q(10), ...
                'ki_pll',   q(11), ...
                'kpi_p',    q(12), ...
                'kii_p',    q(13), ...
                'kpi_n',    q(14), ...
                'kii_n',    q(15), ...
                'Vod0',     obj.Vod0_, ...
                'Voq0',     obj.Voq0_, ...
                'Iod0',     obj.Iod0_, ...
                'Ioq0',     obj.Ioq0_, ...
                'Ild0',     obj.Ild0_, ...
                'Ilq0',     obj.Ilq0_, ...
                'P0',       obj.P0_, ...
                'Q0',       obj.Q0_);
        end

        % DDSRF in positive-frame coordinates.
        % Let Tpn = R(2*theta) and Tnp = Tpn.'. With x_n_minus denoting the
        % negative-frame LPF state, the conventional DDSRF equations are
        %   dx_p/dt       = wf*(u - Tnp*x_n_minus - x_p)
        %   dx_n_minus/dt = wf*(Tpn*u - Tpn*x_p - x_n_minus).
        % Define z_n = Tnp*x_n_minus to store the negative-sequence state in
        % the positive frame. Since
        %   dTnp/dt = -2*w*J*Tnp,  J = [0,-1;1,0],
        % the equivalent equations become
        %   dx_p/dt = wf*(u - z_n - x_p)
        %   dz_n/dt = wf*(u - x_p - z_n) - 2*w*J*z_n.
        % In this function x_n denotes z_n. Both decoupled outputs are in the
        % positive frame: y_p = u - x_n and y_n = u - x_p.
        function [dx_p,dx_n,y_p,y_n] = DDSRF(~,u,x_p,x_n,w,wf)
            J = [0,-1;1,0];
            dx_p = wf*(u-x_n-x_p);
            dx_n = wf*(u-x_p-x_n)-2*w*J*x_n;
            y_p = u-x_n;
            y_n = u-x_p;
        end
    end
end
