% This function creates a descriptor state space model for devices
% connected to buses (such as synchronous generator, PLL-controlled VSI
% ...)

% Author(s): Yitong Li, Yunjie Gu

%% References:

% Y. Gu, Y. Li, Y. Zhu, T. C. Green, "Impedance-based whole-system modeling
% for a composite grid via embedding of frame dynamics." IEEE Transactions
% on Power Systems.

% Y. Li, Y. Gu, T. C. Green, "Interpreting frame tramsformations in ac
% systems as diagonolization of harmonic transfer function." IEEE
% Transctions on Circuit and Systems.


%% function
function [GmObj,GmDSS,DevicePara,DeviceEqui,DeviceDiscreDamping] ...
        = DeviceModel_Create(varargin) 

%% load arguments and common symbols
for n = 1:length(varargin)
    if(strcmpi(varargin{n},'Para'))
        Para = varargin{n+1};
    elseif(strcmpi(varargin{n},'Flow'))
        PowerFlow = varargin{n+1};
    elseif(strcmpi(varargin{n},'Type'))
        Type = varargin{n+1};
    elseif(strcmpi(varargin{n},'Freq'))
        w0 = 2*pi*varargin{n+1};
    elseif(strcmpi(varargin{n},'Ts'))
        Ts = varargin{n+1};
    end
end

% Set default values
try
    w0;
catch
    w0 = 2*pi*50;   %default base frequency
end

try
    Type;
catch
    Type = 0;
end

try 
    PowerFlow;
catch
    PowerFlow = [-1,0,1,0,w0];  %[P Q V xi omega]
    % note the frequency in flow can be different to 'freq' in parameter
    % the frequency in flow is steady-state frequency
    % 'freq' is the nominal frequency only used for default parameters
    % 'freq' is useless if 'para' and 'flow' are set by users
end

%% Create an object
switch floor(Type/10)
    
    % ### Synchronous generator
    case 0      % Type 0-9
        Device = Class_SynchronousMachine;
        Device.Para  = [Para.J;
                        Para.D;
                        Para.L;
                        Para.R;
                        Para.w0];       % (5)
                    
    % ### Grid-Following VSI
    case 1      % Type 10-19
        Device = Class_GridFollowingVSI;
        Device.Para = [Para.C_dc;
                       Para.V_dc;
                       Para.kp_v_dc;
                       Para.ki_v_dc;
                       Para.kp_pll; 	% (5)
                       Para.ki_pll;
                       Para.tau_pll;
                       Para.kp_i_dq;
                       Para.ki_i_dq;
                       Para.k_pf;    	% (10)
                       Para.L;
                       Para.R;
                       Para.w0;
                       Para.Gi_cd]; 	% (14)
                   
    % ### Grid-Forming VSI
    case 2  % Type 20-29
        Device = Class_GridFormingVSI;
        Device.Para = [Para.Lf;
                       Para.Rf;
                       Para.Cf;
                       Para.Lc;
                       Para.Rc;         % (5)
                       Para.Dw;
                       Para.Dv;
                       Para.Tf;
                       Para.P0;
                       Para.Q0;         % (10)
                       Para.w0;
                       Para.v_od0;
                       Para.v_oq0;
                       Para.kp_i_ldq;
                       Para.ki_i_ldq;   % (15)
                       Para.kp_v_odq;
                       Para.ki_v_odq;
                       Para.Gi_cd;
                       Para.Gv_cd;
                       Para.Fv;         % (20)
                       Para.Fi;
                       Para.Rov;
                       Para.Xov];       % (23)
                   
% 	% ### Passive load
%     case 9 % Type 90-99
%         Device = Class_PassiveLoad;
%         Device.Para = [Para.W0;
%                        Para.Connection;
%                        Para.R;
%                        Para.L];
    % ### Infinite Bus
    case 9 % Type 90-99
    Device = Class_SynchronousMachine;
    Device.Para  = [Para.J;
                    Para.D;
                    Para.L;
                    Para.R;
                    Para.w0];       % (5)
                   
    % ### Floating Bus
    case 10
        Device = Class_FloatingBus;
        Device.Para = [];
       
    % ### Otherwise
    otherwise
        error(['Error: device type']);
end

%% Calculate the linearized state space model
Device.DeviceType = Type;
Device.PowerFlow = PowerFlow;           % Power flow data
Device.SetString(Device);               % Set strings automatically
Device.Equilibrium(Device);             % Calculate the equilibrium
[x_e,u_e,y_e,xi] = Device.ReadEquilibrium(Device);
Device.Linearization(Device,x_e,u_e); 	% Linearize the model
Device.ConstructSS(Device);             % Construct the state space model

[~,ModelSS] = Device.ReadSS(Device);
[StateString,InputString,OutputString] = Device.ReadString(Device);

v_d = u_e(1);
v_q = u_e(2);
i_d = y_e(1);
i_q = y_e(2);

% Get the swing frame system model
Gm = ModelSS;   

% Output
DevicePara = Device.Para; 
DeviceEqui = {x_e,u_e,y_e,xi};

% Output the discretization damping resistance for simulation use
if floor(Type/10) <= 9
    Ak = ModelSS.A;
    Ck = ModelSS.C;
    Bk = ModelSS.B;
    Wk = inv(eye(length(Ak)) - Ts/2*Ak);
    MatrixY = Ts/2*Ck*Wk*Bk;
    MatrixY = MatrixY(1:2,1:2);
    MatrixR = inv(MatrixY);
    DeviceDiscreDamping = MatrixR(1,1);
else
    DeviceDiscreDamping = -1;
end

%% Impedance transformation: local swing frame dq -> local steady frame d'q'
if floor(Type/10) > 9
    Se = Gm;
else
    
% Effect of the frame perturbation on arbitrary signal udq:
% Complex vector dq frame:
% [ud'q'+] = [exp(j*epsilon), 0              ] * [udq+]
% [ud'q'-]   [0,              exp(-j*epsilon)]   [udq-]
% where "epsilon" is the angle that swing frame axes leads steady frame
% axes, which is also the angle that swing-frame signals lag the
% steady-frame signals.
%
% Linearized effect:
% Steady = Swing + Equilibrium * epsilon with epsilon = omega/s.
% Complex vector dq frame:
% [ud'q'+] = [udq+] + [ judq+0] * epsilon
% [ud'q'-]   [udq-]   [-judq-0]
% Transfer matrix dq frame:
% [ud'] = [ud] + [-uq0] * epsilon
% [uq']   [uq]   [ ud0]

% Get the steady-state operating points
V0 = [-v_q ; v_d];
I0 = [-i_q ; i_d];

% Integration for omega and unit gain for others
% Notes: New system has one more new output "w/s", which is added to the
% head of the original output vector. New system also has a new state
% "epsilon", which is added to the head of the original state vector.
for i = 1:length(OutputString)
    if strcmp(OutputString(i),'w')
        ind_w = i;
        break
    end
end
[~,~,ly1_w] = dss_GetDim(Gm);
Aw = 0;
Bw = zeros(1,ly1_w);   Bw(ind_w) = 1;
Cw = [1;zeros(ly1_w,1)];
Dw = [zeros(1,ly1_w);eye(ly1_w)];
Se = ss(Aw,Bw,Cw,Dw);
Se = series(Gm,Se);
StateString = [{'epsilon'},StateString];

% Embed the frame dynamics
Sfb = ss([],[],[],V0);
Se = feedback(Se,Sfb,[1,2],[1]);            % v = v' - V0 * w/s

[~,~,ly2_w] = dss_GetDim(Se);
Sff = ss([],[],[],[[I0;zeros((ly2_w-3),1)],eye(ly2_w-1)]);
Se = series(Se,Sff);                        % i' = i + I0 * w/s

end

%% Impedance transformation: local steady frame d'q' -> global steady frame D'Q'
if floor(Type/10) > 9
    % Device is a passive load
    Se = Gm;
else
    
% Effect of frame alignment on arbitrary signal udq:
% Complex vector dq frame:
% [uD'Q'+] = [exp(j*xi), 0         ] * [ud'q'+]
% [uD'Q'-]   [0,         exp(-j*xi)]   [ud'q'-]
%             \-------------------/
%                      Txi
% Transfer matrix dq frame:
% [uD'] = [cos(xi), -sin(xi)] * [ud']
% [uQ']   [sin(xi),  cos(xi)]   [uq']
%          \---------------/
%                 Txi
% where "xi" is the angle that local steady frame axes lead the global
% steady frame axes, which is also the local signals lag the global
% signals.
%
% Effect on the transfer function: similarity transform
% G_D'Q' = Txi*G_d'q'*inv(Txi)

% Get the dimension
[~,lu_xi,ly_xi] = dss_GetDim(Se);

% Transform matrix
Txi = [cos(xi),-sin(xi);
       sin(xi), cos(xi)];

% Connect Txi to system right-hand-side output: 
% Left multiplication of matrix, i.e., Txi*G;
Txi_left = blkdiag(Txi,eye(ly_xi-2));
Sxi_left = ss([],[],[],Txi_left);
Se = series(Se,Sxi_left);

% Connect inv(Txi) to system left-hand-side input: 
% Right multiplication of matrix, i.e., G*inv(Txi)
Txi_right = blkdiag(inv(Txi),eye(lu_xi-2));
Sxi_right= ss([],[],[],Txi_right);
Se = series(Sxi_right,Se);

end

%% Get the descritpor state space model
% Set the SS system back
Gm = Se;

% Change SS system to DSS system
An = Gm.A; Bn = Gm.B; Cn = Gm.C; Dn = Gm.D; 
En = eye(size(An));
GmDSS = dss(An,Bn,Cn,Dn,En);

%% Get the object model
% Create an object
GmObj = Class_Model_DSS;

% Load the model
GmObj.LoadModel(GmObj,GmDSS);

% Get the strings
GmObj.WriteString(GmObj,StateString,InputString,OutputString);

% Check dimension mismatch
obj_CheckDim(GmObj);

end