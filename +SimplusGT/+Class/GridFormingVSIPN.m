% This class defines a grid-forming VSI with positive- and
% negative-sequence control.
%
% Author(s): Wanjun Gao

%% Notes
%
% The model uses load convention and admittance form. DDSRF states are
% represented in the positive frame:
%   dx_p = wf * (u - x_n - x_p)
%   dx_n = wf * (u - x_p - x_n) - 2*w*J*x_n
% Here x_n is the negative-sequence filter state mapped into the positive
% controller frame. With J = [0,-1;1,0] and dq rotation R, for
% z_n = R(-2*theta)*z_n_minus,
%   dz_n/dt = R(-2*theta)*dz_n_minus/dt - 2*w*J*z_n.

%% Class

classdef GridFormingVSIPN < SimplusGT.Class.ModelAdvance

    % Cached equilibrium values
    properties(Access = protected)
        Vod0_ = 0; Voq0_ = 0; Iod0_ = 0; Ioq0_ = 0;
        Ild0_ = 0; Ilq0_ = 0; Uvd0_ = 0; Uvq0_ = 0;
        P0_ = 0; Q0_ = 0; U0_ = 0;
    end

    methods
        % Constructor
        function obj = GridFormingVSIPN(varargin)
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
                     'V_mag_lpf','Q_p_lpf','u_d_ref_p_i', ...
                     'x_del_d','x_del_q','w','theta'};
            % P0 is retained for compatibility with the Simulink template.
            Input  = {'v_d','v_q','P0'};
            Output = {'i_d','i_q','w','theta'};
        end

        % Calculate the equilibrium
        function [x_e,u_e,xi] = Equilibrium(obj)
            P  = obj.PowerFlow(1);
            Q  = obj.PowerFlow(2);
            V  = obj.PowerFlow(3);
            xi = obj.PowerFlow(4);
            w  = obj.PowerFlow(5);

            p = obj.GetPara(obj);
            v_gdq = V;
            i_od = P/V;
            i_oq = -Q/V;

            v_odq = v_gdq ...
                    -(i_od+1i*i_oq)*(p.Rc+1i*w*p.Lc);
            i_cdq = v_odq/(p.Rd+1/(1i*w*p.Cf));
            i_ldq = (i_od+1i*i_oq)-i_cdq;
            v_cdq = v_odq-p.Rd*i_cdq;
            e_dq_ss = v_odq-i_ldq*(p.Rf+1i*w*p.Lf);


            u_vdq = v_odq-(p.Rv+1i*p.Xv)*i_ldq;
            % Retain both virtual-voltage components in the global dq frame.
            i_odq = i_od+1i*i_oq;

            v_gd = real(v_gdq); v_gq = imag(v_gdq);
            v_od = real(v_odq); v_oq = imag(v_odq);
            i_od = real(i_odq); i_oq = imag(i_odq);
            i_ld = real(i_ldq); i_lq = imag(i_ldq);
            v_cd = real(v_cdq); v_cq = imag(v_cdq);
            e_d_ss = real(e_dq_ss);
            e_q_ss = imag(e_dq_ss);


            obj.Vod0_ = v_od; obj.Voq0_ = v_oq;
            obj.Iod0_ = i_od; obj.Ioq0_ = i_oq;
            obj.Ild0_ = i_ld; obj.Ilq0_ = i_lq;
            obj.Uvd0_ = real(u_vdq);
            obj.Uvq0_ = imag(u_vdq);
            obj.P0_ = -(v_od*i_od+v_oq*i_oq);
            obj.Q0_ = -(v_oq*i_od-v_od*i_oq);
            obj.U0_ = sqrt(v_od^2+v_oq^2);


            den_v = p.Rv^2+p.Xv^2;
            Y_R = p.Rv/den_v;
            Y_X = p.Xv/den_v;
            i_ld_r_p = -(Y_R*(obj.Uvd0_-v_od)+Y_X*(obj.Uvq0_-v_oq));
            i_lq_r_p = -(Y_R*(obj.Uvq0_-v_oq)-Y_X*(obj.Uvd0_-v_od));

            x_e = zeros(29,1);
            x_e(1) = i_ld; x_e(2) = i_lq;
            x_e(3) = e_d_ss-v_od-p.kpi*(i_ld_r_p-i_ld);
            x_e(4) = e_q_ss-v_oq-p.kpi*(i_lq_r_p-i_lq);
            x_e(7) = v_cd; x_e(8) = v_cq;
            x_e(9) = i_od; x_e(10) = i_oq;

            % At balanced equilibrium, all negative-sequence states are zero.
            x_e(11:12) = [v_od;v_oq];
            x_e(15:16) = [i_ld;i_lq];
            x_e(19:20) = [i_od;i_oq];

            x_e(23) = sqrt(v_od^2+v_oq^2);
            x_e(24) = v_od*i_oq-v_oq*i_od;
            % The Q-loop regulates virtual-voltage magnitude.
            Uv0_mag = sqrt(obj.Uvd0_^2+obj.Uvq0_^2);
            x_e(25) = Uv0_mag-obj.U0_;
            x_e(26) = e_d_ss*p.Td/2;
            x_e(27) = e_q_ss*p.Td/2;
            x_e(28) = w;
            x_e(29) = xi;


            u_e = [v_gd;v_gq;0];
        end

        % State-space model
        function Output = StateSpaceEqu(obj,x,u,CallFlag)
            p = obj.GetPara(obj);

            v_gd = u(1);
            v_gq = u(2);

            i_ld = x(1);     i_lq = x(2);
            i_ld_p_i = x(3); i_lq_p_i = x(4);
            i_ld_n_i = x(5); i_lq_n_i = x(6);
            v_cd = x(7);     v_cq = x(8);
            i_od = x(9);     i_oq = x(10);

            v_od = v_cd + p.Rd*(i_od-i_ld);
            v_oq = v_cq + p.Rd*(i_oq-i_lq);

            % DDSRF filter states for terminal voltage, converter-side
            % current, and grid-side current. Both sequence channels are
            % stored in the positive controller frame. Each signal uses two
            % states per sequence, giving 12 DDSRF states in total.
            x_v_p  = x(11:12); x_v_n  = x(13:14);
            x_il_p = x(15:16); x_il_n = x(17:18);
            x_io_p = x(19:20); x_io_n = x(21:22);

            V_mag_lpf = x(23);  Q_p_lpf = x(24);
            u_d_ref_p_i = x(25);
            x_del_d = x(26); x_del_q = x(27);
            w = x(28); theta = x(29);

            if CallFlag == 1
                % State equations: dx/dt = f(x,u)
                w0 = p.w0;

                [dx_v_p,dx_v_n,v_p,~] = obj.DDSRF( ...
                    obj,[v_od;v_oq],x_v_p,x_v_n,w,p.w_lpf);
                [dx_il_p,dx_il_n,il_p,il_n] = obj.DDSRF( ...
                    obj,[i_ld;i_lq],x_il_p,x_il_n,w,p.w_lpf);
                [dx_io_p,dx_io_n,io_p,~] = obj.DDSRF( ...
                    obj,[i_od;i_oq],x_io_p,x_io_n,w,p.w_lpf);

                v_od_p_f = v_p(1); v_oq_p_f = v_p(2);
                i_ld_p_f = il_p(1); i_lq_p_f = il_p(2);
                i_od_p_f = io_p(1); i_oq_p_f = io_p(2);

                V_mag = sqrt(v_od_p_f^2+v_oq_p_f^2);
                dV_mag_lpf = p.w_amp*(V_mag-V_mag_lpf);
                Q_p_cal = -(v_oq_p_f*i_od_p_f-v_od_p_f*i_oq_p_f);
                P_p_cal = -(v_od_p_f*i_od_p_f+v_oq_p_f*i_oq_p_f);
                dQ_p_lpf = p.w_lpf*(Q_p_cal-Q_p_lpf);

                dw = (w0+p.m*(p.P0-P_p_cal)-w)*p.w_lpf;
                dtheta = w;

                Q_ref = p.Q0+p.kp*(p.U0-V_mag_lpf);
                e_Q = Q_ref-Q_p_lpf;
                u_mag_ref_p = p.U0+p.kpq*e_Q+u_d_ref_p_i;
                du_d_ref_p_i = p.kiq*e_Q;

                % Apply the scalar command along the steady virtual-voltage direction.
                Uv0_mag = sqrt(p.Uvd0^2+p.Uvq0^2);
                u_d_ref_p = (p.Uvd0/Uv0_mag)*u_mag_ref_p;
                u_q_ref_p = (p.Uvq0/Uv0_mag)*u_mag_ref_p;

                den_v = p.Rv^2+p.Xv^2;
                Y_R = p.Rv/den_v; Y_X = p.Xv/den_v;
                i_ld_r_p = -(Y_R*(u_d_ref_p-v_od_p_f) ...
                              +Y_X*(u_q_ref_p-v_oq_p_f));
                i_lq_r_p = -(Y_R*(u_q_ref_p-v_oq_p_f) ...
                              -Y_X*(u_d_ref_p-v_od_p_f));

                % Positive-sequence current controller.
                e_d_p = i_ld_p_i+(i_ld_p_f-i_ld_r_p)*p.kpi ...
                        +p.miu*x_v_p(1);
                e_q_p = i_lq_p_i+(i_lq_p_f-i_lq_r_p)*p.kpi ...
                        +p.miu*x_v_p(2);
                di_ld_p_i = (i_ld_p_f-i_ld_r_p)*p.kii;
                di_lq_p_i = (i_lq_p_f-i_lq_r_p)*p.kii;

                % Map the negative-sequence PI integrator into the positive
                % frame. Differentiating R(-2*theta)*z_n_minus produces the
                % additional -2*w*J*z_n term.
                J = [0,-1;1,0];
                e_n_pos = [i_ld_n_i;i_lq_n_i]+p.kpi*il_n+p.miu*x_v_n;
                di_n_i = p.kii*il_n-2*w*J*[i_ld_n_i;i_lq_n_i];
                di_ld_n_i = di_n_i(1);
                di_lq_n_i = di_n_i(2);

                e_d_ref = e_d_p+e_n_pos(1)+(1-p.miu)*v_od;
                e_q_ref = e_q_p+e_n_pos(2)+(1-p.miu)*v_oq;

                dx_del_d = -(2/p.Td)*x_del_d+e_d_ref;
                dx_del_q = -(2/p.Td)*x_del_q+e_q_ref;
                e_d = (4/p.Td)*x_del_d-e_d_ref;
                e_q = (4/p.Td)*x_del_q-e_q_ref;

                di_ld = (v_od-e_d-p.Rf*i_ld+w*p.Lf*i_lq)/p.Lf;
                di_lq = (v_oq-e_q-p.Rf*i_lq-w*p.Lf*i_ld)/p.Lf;
                dv_cd = ((i_od-i_ld)+w*p.Cf*v_cq)/p.Cf;
                dv_cq = ((i_oq-i_lq)-w*p.Cf*v_cd)/p.Cf;
                di_od = (v_gd-v_od-p.Rc*i_od+w*p.Lc*i_oq)/p.Lc;
                di_oq = (v_gq-v_oq-p.Rc*i_oq-w*p.Lc*i_od)/p.Lc;

                Output = [di_ld;di_lq;di_ld_p_i;di_lq_p_i; ...
                          di_ld_n_i;di_lq_n_i;dv_cd;dv_cq;di_od;di_oq; ...
                          dx_v_p;dx_v_n;dx_il_p;dx_il_n;dx_io_p;dx_io_n; ...
                          dV_mag_lpf;dQ_p_lpf;du_d_ref_p_i; ...
                          dx_del_d;dx_del_q;dw;dtheta];
            elseif CallFlag == 2
                % Output equations: y = g(x,u)
                Output = [i_od;i_oq;w;theta];
            else
                error('GridFormingVSIPN:StateSpaceEqu:InvalidCallFlag', ...
                      'CallFlag must be 1 (f) or 2 (g).');
            end
        end

        % Convert the parameter vector to physical quantities
        function p = GetPara(obj)
            q = obj.Para;
            w0 = q(18);

            p = struct( ...
                'w0',      w0, ...
                'Td',      1.5*obj.Ts/2, ...
                'Lf',      q(1)/w0, ...
                'Rf',      q(2), ...
                'Cf',      q(3)/w0, ...
                'Rd',      q(4), ...
                'Lc',      q(5)/w0, ...
                'Rc',      q(6), ...
                'm',       q(7)*w0, ...
                'w_lpf',   q(8)*2*pi, ...
                'kp',      q(9), ...
                'kpq',     q(10), ...
                'kiq',     q(11), ...
                'w_amp',   q(12)*2*pi, ...
                'kpi',     q(13), ...
                'kii',     q(14), ...
                'miu',     q(15), ...
                'Rv',      q(16), ...
                'Xv',      q(17), ...
                'Vod0',    obj.Vod0_, ...
                'Voq0',    obj.Voq0_, ...
                'Iod0',    obj.Iod0_, ...
                'Ioq0',    obj.Ioq0_, ...
                'Ild0',    obj.Ild0_, ...
                'Ilq0',    obj.Ilq0_, ...
                'Uvd0',    obj.Uvd0_, ...
                'Uvq0',    obj.Uvq0_, ...
                'P0',      obj.P0_, ...
                'Q0',      obj.Q0_, ...
                'U0',      obj.U0_);
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
