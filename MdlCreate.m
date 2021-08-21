% create a linearized model for a converter or a machine
% output: transfer function, state space, symbolic real, symbolic complex

function [Gtf, Gss, Gc, Gr] = MdlCreate(varargin) 

% ### motor convention, d lagging q, with d aligned to voltage
% ### model can use either be perunit or SI unit, but time related
% variables (omega) should never be in perunit, which means:
% current = Ib; voltage = Vb; frequency(omega) = 1;
% resistance = Vb/Ib; inductance = Vb/Ib; capacitance = Ib/Vb;
% power = torque = Ib*Vb; flux linkage = Vb;
% damping torque = inertia = Ib*Vb;

%% load arguments and common symbols
for n = 1:length(varargin)
    if(strcmpi(varargin{n},'para'))
        para = varargin{n+1};
    elseif(strcmpi(varargin{n},'flow'))
        flow = varargin{n+1};
    elseif(strcmpi(varargin{n},'type'))
        type = varargin{n+1};
    elseif(strcmpi(varargin{n},'freq'))
        w0 = 2*pi*varargin{n+1};
    end
end

try 
    w0;
catch
    w0 = 2*pi*50; %default base frequency
end

try 
    type;
catch
    type = 0; %0-9 generator %10-19 converter
end

if floor(type/10) == 0
    try 
        para;
    catch
        % default parameter in perunit
        para.J = (3.5)*2/w0^2;  % Jpu=J/Pb=[1/2*J*w0^2/Pb]*2/w0^2, [MWs/MW] 
        para.D = (1)/w0^2;   % Dpu=dTpu/dw=dPpu/dw/w0=[dP%/dw%]/w0^2, [%/%]
        para.L = (0.05)/w0;
        para.R = (0.01);
    end
elseif floor(type/10) == 1
    try 
        para;
    catch
        % default parameter in perunit
        para.V_dc = 2.5;
        para.C_dc = 2*0.1*para.V_dc^2;
        para.kp_v_dc = para.V_dc*para.C_dc*(10*2*pi);
        para.ki_v_dc = para.kp_v_dc*(10*2*pi)/4;
        para.kp_pll = 2*2*pi;
        para.ki_pll = para.kp_pll * (2*2*pi)/4; 
        para.tau_pll = 1/(2*pi*200);
        para.k_pf = 0;
        para.L = 0.05/(w0);
        para.R = 0.01;
        para.kp_i_dq = para.L * (500*2*pi);
        para.ki_i_dq = para.kp_i_dq *(500*2*pi)/4;
    end
end

try 
    flow;
catch
    flow = [-1,0,1,0,w0];  %[P Q V xi omega]
    % note the frequency in flow can be different to 'freq' in parameter
    % the frequency in flow is steady-state frequency
    % 'freq' is the nominal frequency only used for default parameters
    % 'freq' is useless if 'para' and 'flow' are set by users
end

s = sym('s');
omega = sym('omega');


%% state space for a synchronous generator
if floor(type/10) == 0  %0-9

    % ### state variables
    i_d = sym('i_d'); 
    i_q = sym('i_q');
    i_ex = sym('i_ex'); %excitation (field) current

    % ### input and output
    v_d = sym('v_d');
    v_q = sym('v_q');
    v_ex= sym('v_ex');  %excitation (field) voltage
    T_m = sym('T_m');

    % ### parameters
    psi_f = sym('psi_f');
    J = sym('J');
    D = sym('D');
    L = sym('L');
    R = sym('R');

    % ### ancillary equations
    % empty

    % ### output equations
    % empty

    % ### state equations
    di_d = (v_d - R*i_d + omega * (L*i_q-psi_f))/L;
    di_q = (v_q - R*i_q - omega * (L*i_d)      )/L;
    domega = (psi_f * i_d - T_m - D*omega)/J;

    % ### Jacobian
    x = [i_d i_q omega].';
    f = [di_d di_q domega].';
    u = [T_m   v_ex v_d v_q].';
    y = [omega i_ex i_d i_q].';

    As = jacobian(f,x);
    Bs = jacobian(f,u);
    Cs = jacobian(y,x);
    Ds = jacobian(y,u);

    % ### get parameters
    J = para.J;
    D = para.D;
    L = para.L;
    R = para.R;

    % ### get operating points
    P = flow(1);
    Q = flow(2);
    V = flow(3);
    xi = flow(4);
    omega = flow(5);

    i_D = P/V;
    i_Q = -Q/V;
    i_DQ = i_D + 1j*i_Q;
    e_DQ = V - i_DQ * (R + 1j*L*omega);
    arg_e = angle(e_DQ);
    abs_e = abs(e_DQ);
    xi = xi + arg_e;

    v_dq = V * exp(-1j*arg_e);
    i_dq = i_DQ * exp(-1j*arg_e);
    v_d = real(v_dq);
    v_q = imag(v_dq);
    i_d = real(i_dq);
    i_q = imag(i_dq);
    psi_f = abs_e/omega;
    T_m = psi_f * i_d - D*omega;

elseif floor(type/10) == 1  %10-19
 
    % ### state variables
    i_d = sym('i_d'); 
    i_q = sym('i_q');
    i_d_i = sym('i_d_i');
    i_q_i = sym('i_q_i');
    v_dc = sym('v_dc');
    v_dc_i = sym('v_dc_i');
    omega_pll_i = sym('omega_pll_i');

    % ### input and output
    v_d  = sym('v_d');
    v_q  = sym('v_q');
    P_dc = sym('P_dc');
    ang_r= sym('ang_r');
    % i_d = sym('i_d');
    % i_q = sym('i_q');
    % omega = sym('omega');

    % ### parameters
    C_dc = sym('C_dc');
    V_dc = sym('V_dc');
    kp_v_dc = sym('kp_v_dc');
    ki_v_dc = sym('ki_v_dc');
    kp_pll = sym('kp_pll');
    ki_pll = sym('ki_pll');
    tau_pll = sym('tau_pll');
    kp_i_dq = sym('kp_i_dq');
    ki_i_dq = sym('ki_i_dq');
    k_pf = sym('k_pf');
    L = sym('L');
    R = sym('R');

    % ### ancillary equations
    i_d_r = (V_dc - v_dc)*kp_v_dc + v_dc_i;
    %i_q_r = i_d_r * -k_pf;  %constant pf control, PQ node in power flow
    i_q_r = sym('i_q_r');    %constant q control, PQ/PV node in power flow
    e_d = (i_d - i_d_r)*kp_i_dq + i_d_i;
    e_q = (i_q - i_q_r)*kp_i_dq + i_q_i;
    e_ang = atan2(v_q,v_d) - ang_r;
      
    % ### output equations
    % empty

    % ### state equations
    dv_dc = (e_d*i_d + e_q*i_q - P_dc)/v_dc/C_dc;
    dv_dc_i = (V_dc - v_dc)*ki_v_dc;
    di_d_i = (i_d - i_d_r)*ki_i_dq;
    di_q_i = (i_q - i_q_r)*ki_i_dq;
    di_d = (v_d - R*i_d + omega * L*i_q - e_d)/L;
    di_q = (v_q - R*i_q - omega * L*i_d - e_q)/L;
    domega_pll_i = e_ang * ki_pll; 
    domega = (omega_pll_i + e_ang * kp_pll - omega)/tau_pll;

    % ### Jacobian
    x = [i_d i_q i_d_i i_q_i v_dc v_dc_i omega_pll_i omega].';
    f = [di_d di_q di_d_i di_q_i dv_dc dv_dc_i domega_pll_i domega].';
      
    u = [ang_r P_dc v_d v_q].';
    y = [omega v_dc i_d i_q].';

    As = jacobian(f,x);
    Bs = jacobian(f,u);
    Cs = jacobian(y,x);
    Ds = jacobian(y,u);

    % ### get parameters
    C_dc = para.C_dc;
    V_dc = para.V_dc;
    kp_v_dc = para.kp_v_dc;
    ki_v_dc = para.ki_v_dc;
    kp_pll = para.kp_pll;
    ki_pll = para.ki_pll;
    tau_pll = para.tau_pll;
    kp_i_dq = para.kp_i_dq;
    ki_i_dq = para.ki_i_dq;
    k_pf = para.k_pf;
    L = para.L;
    R = para.R;

    % ### get operating points
    P = flow(1);
    Q = flow(2);
    V = flow(3);
    xi = flow(4);
    omega = flow(5);

    i_d = P/V;
    i_q = -Q/V;
    v_d = V;
    v_q = 0;
    i_dq = i_d + 1j*i_q;
    v_dq = v_d + 1j*v_q;
    e_dq = v_dq - i_dq * (R + 1j*L*omega);
    e_d = real(e_dq);
    e_q = imag(e_dq);
    i_d_i = e_d;
    i_q_i = e_q;
    i_d_r = i_d;
    i_q_r = i_q;
    omega_pll_i = omega;
    v_dc_i = i_d;
    v_dc = V_dc;
    P_dc = e_d*i_d + e_q*i_q;
    ang_r= 0;
      
end

%% impedance transformation

% ### get numerical ABCD
An = double(subs(As));
Bn = double(subs(Bs));
Cn = double(subs(Cs));
Dn = double(subs(Ds));
Sn = ss(An,Bn,Cn,Dn);

% ### transfer function models
% embed frame dynamics
Kv = [-v_q ; v_d];
Ki = [-i_q ; i_d];
% integration for omega and unit gain for all
Se = ss(0,[1 0 0 0],[1;0;0;0;0],[zeros(1,4);eye(4)]);
Se = series(Sn,Se);
Se = feedback(Se,ss([],[],[],Kv),[3,4],1);
Se = series(Se,ss([],[],[],[[0;0;Ki],eye(4)]));
Txi = blkdiag(eye(2),[cos(xi) -sin(xi);sin(xi) cos(xi)]);
Sxi = ss([],[],[],Txi);
Sxiv= ss([],[],[],Txi^(-1));
Se = series(Sxiv,Se);
Se = series(Se,Sxi);
Gss = Se;
Gtf = tf(Se);

% ### symbolic models
I = eye(length(Gss.A));
Tj = blkdiag(eye(2),[1 1j;1 -1j]);
Gr = Gss.C *(s*I-Gss.A)^(-1)* Gss.B + Gss.D;   %real
Gc = Tj *Gr* Tj^(-1);                          %complex
%Gc = simplify(Gc);
%Gr = simplify(Gr);

end