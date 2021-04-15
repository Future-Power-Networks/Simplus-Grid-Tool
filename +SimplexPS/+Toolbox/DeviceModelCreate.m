% This function creates the descriptor state space model for devices
% connected to buses (such as synchronous generators, voltage source
% inverters, ...)

% Author(s): Yitong Li, Yunjie Gu, Yue Zhu

%% References:

% Y. Gu, Y. Li, Y. Zhu, T. C. Green, "Impedance-based whole-system modeling
% for a composite grid via embedding of frame dynamics." IEEE Transactions
% on Power Systems, Early Access.

% Y. Li, Y. Gu, T. C. Green, "Interpreting frame tramsformations in ac
% systems as diagonolization of harmonic transfer function." IEEE
% Transctions on Circuit and Systems, 2020.

%% function
function [GmObj,GmDSS,DevicePara,DeviceEqui,DeviceDiscreDampingResistor,StateString,InputString,OutputString] ...
        = DeviceModelCreate(varargin) 

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
Flag_SwitchInOut = 0;   % Default: do not need to switch input and output
switch floor(Type/10)
    
    % ### Synchronous generator
    case 0      % Type 0-9
        Device = SimplexPS.Class.SynchronousMachine('DeviceType',Type);
        Device.Para  = [Para.J;
                        Para.D;
                        Para.L;
                        Para.R;
                        Para.w0];       % (5)
                    
    % ### Grid-following inverter
    case 1      % Type 10-19
        Device = SimplexPS.Class.GridFollowingVSI('DeviceType',Type);
        Device.Para = [
                            Para.Vdc;
                            Para.Cdc;
                            Para.wL;
                            Para.R;
                            Para.fvdc;
                            Para.fpll;
                            Para.fidq;
%                        Para.C_dc;
%                        Para.V_dc;
%                        Para.kp_v_dc;
%                        Para.ki_v_dc;
%                        Para.kp_pll; 	% (5)
%                        Para.ki_pll;
%                        Para.tau_pll;
%                        Para.kp_i_dq;
%                        Para.ki_i_dq;
%                        Para.k_pf;    	% (10)
%                        Para.L;
%                        Para.R;
%                        Para.w0;
%                        Para.Gi_cd
                       ]; 	% (14)
                   
    % ### Grid-forming inverter
    case 2  % Type 20-29
        if Type == 29
            Device = SimplexPS.Class.GridFormingVSI_Detuned('DeviceType',Type);
        else
            Device = SimplexPS.Class.GridFormingVSI('DeviceType',Type);
        end
        Device.Para = [
                        Para.wLf;
                        Para.Rf;
                        Para.wCf;
                        Para.wLc;
                        Para.Rc;
                        Para.Xov;
                        Para.Dw;
                        Para.fdroop;
                        Para.fvdc;
                        Para.fidq;
%                        Para.Lf;
%                        Para.Rf;
%                        Para.Cf;
%                        Para.Lc;
%                        Para.Rc;
%                        Para.Xov
%                        Para.Dw;
%                        Para.wf;
%                        Para.w_v_odq;
%                        Para.w_i_ldq;
%                        Para.w0                       
                       ];
    
    % ### Synchronous machine full model
    case 3  % Type 30-39
        if Type==30
            Device = SimplexPS.Class.SynchronousMachineFull_SM('DeviceType',Type);
        elseif Type==31
            Device = SimplexPS.Class.SynchronousMachineFull_SMAVR('DeviceType',Type);
        elseif Type==32
            Device = SimplexPS.Class.SynchronousMachineFull_SMAVRPSS('DeviceType',Type);
        elseif Type == 33
            Device = SimplexPS.Class.SynchronousMachineFull_SM_GOV('DeviceType',Type);
        elseif Type == 34
            Device = SimplexPS.Class.SynchronousMachineFull_SMAVRGOV('DeviceType',Type);
        elseif Type == 35
            Device = SimplexPS.Class.SynchronousMachineFull_SMAVRPSSGOV('DeviceType',Type);
        elseif Type == 38
            Device = SimplexPS.Class.SynchronousMachineFull_SMGOV_FS('DeviceType',Type);
        elseif Type == 39
            Device = SimplexPS.Class.SynchronousMachineFull_SMAVRPSSGOV_FS('DeviceType',Type);
        end
                Device.Para = [ Para.X;
                        Para.R;
                        Para.Xd; 
                        Para.Xd1; 
                        Para.Xd2; 
                        Para.Td1; 
                        Para.Td2;
                        Para.Xq;
                        Para.Xq1;
                        Para.Xq2;
                        Para.Tq1;
                        Para.Tq2;
                        Para.H;
                        Para.D;
                        %Para.AVRSel;
                        %Para.PSSSel;
                        %AVR parameters
                        Para.TR;
                        Para.KA;
                        Para.TA;
                        Para.VRmax;
                        Para.VRmin;
                        Para.KE;
                        Para.TE;
                        Para.E1;
                        Para.SE1;
                        Para.E2;
                        Para.SE2;
                        Para.KF;
                        Para.TF;
                        Para.KP;
                        Para.KI;
                        Para.KD;
                        Para.TD;
                        Para.KPSS;
                        Para.TW;
                        Para.T11;
                        Para.T12;
                        Para.T21;
                        Para.T22;
                        Para.T31;
                        Para.T32;
                        Para.VSSmax;
                        Para.VSSmin;
                        Para.Rgov;
                        Para.T1gov;
                        Para.T2gov;
                        Para.T3gov;
                        Para.Dtgov;
                        ];
    % ### Infinite bus
    case 9 % Type 90-99
        Device = SimplexPS.Class.InfiniteBus;
        Device.Para  = [];
        % Because the infinite bus is defined with "i" input and "v" output,
        % they need to be switched finally.
        Flag_SwitchInOut = 1;   
   
    % ### Floating bus
    case 10
        Device = SimplexPS.Class.FloatingBus;
        Device.Para = [];
        
    % ### Otherwise
    otherwise
        error(['Error: device type']);
end

%% Calculate the linearized state space model
Device.DeviceType = Type;                           % Device type
Device.Ts = Ts;                                     % Samping period
Device.PowerFlow = PowerFlow;                       % Power flow data
% if floot(Type/10) == 3
%     Device.PowerFlow=PowerFlow + [0 0 0 0 pi/2];
% end
Device.SetString(Device);                           % Set strings automatically
Device.SetEquilibrium(Device);                      % Calculate the equilibrium
[x_e,u_e,y_e,xi] = Device.GetEquilibrium(Device);   % Get the equilibrium
Device.SetSSLinearized(Device,x_e,u_e);             % Linearize the model

[~,ModelSS] = Device.GetSS(Device);                % Get the ss model
[StateString,InputString,OutputString] ...
    = Device.GetString(Device);                    % Get the string

% Get the swing frame system model
Gm = ModelSS;   

% Output
DevicePara = Device.Para; 
DeviceEqui = {x_e,u_e,y_e,xi};

% Output the discretization damping resistance for simulation use
if floor(Type/10) <= 5
    % CalcRv_old();
    Device.SetDynamicSS(Device,x_e,u_e);
    DeviceDiscreDampingResistor = Device.GetVirtualResistor(Device);
else
    DeviceDiscreDampingResistor = -1;
end

%% Check if the device needs to adjust its frame
if floor(Type/10) >= 9
    
else    
    
%% Frame transformation: local swing frame dq -> local steady frame d'q'
% Effect of the frame perturbation on arbitrary signal udq:
% Complex vector dq frame:
% [ud'q'+] = [exp(j*epsilon), 0              ] * [udq+]
% [ud'q'-]   [0,              exp(-j*epsilon)]   [udq-]
%             \_____________________________/
%                           \/
%            T_epsilon in complex vector form
% Transfer matrix dq frame:
% [ud'] = [cos(epsilon), -sin(epsilon)] * [ud]
% [uq']   [sin(epsilon),  cos(epsilon)]   [uq]
%          \_________________________/
%                      \/
%       T_epsilon in transfer matrix form
% where "epsilon" is the angle that swing frame axes leads steady frame
% axes, which is also the angle that swing-frame signals lag the
% steady-frame signals.
%
% Linearized:
% Steady = Swing + Equilibrium * epsilon with epsilon = omega/s.
% Complex vector dq frame:
% [ud'q'+] = [udq+] + [ judq+0] * epsilon
% [ud'q'-]   [udq-]   [-judq-0]
% Transfer matrix dq frame:
% [ud'] = [ud] + [-uq0] * epsilon
% [uq']   [uq]   [ ud0]

% Get the steady-state operating points
v_d = u_e(1);
v_q = u_e(2);
i_d = y_e(1);
i_q = y_e(2);

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
[~,~,ly1_w] = SimplexPS.SsGetDim(Gm);
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

[~,~,ly2_w] = SimplexPS.SsGetDim(Se);
Sff = ss([],[],[],[[I0;zeros((ly2_w-3),1)],eye(ly2_w-1)]);
Se = series(Se,Sff);                        % i' = i + I0 * w/s

%% Impedance transformation: local steady frame d'q' -> global steady frame D'Q'
% Effect of frame alignment on arbitrary signal udq:
% Complex vector dq frame:
% [uD'Q'+] = [exp(j*xi), 0         ] * [ud'q'+]
% [uD'Q'-]   [0,         exp(-j*xi)]   [ud'q'-]
%             \___________________/
%                      \/
%           T_xi in complex vector form
% Transfer matrix dq frame:
% [uD'] = [cos(xi), -sin(xi)] * [ud']
% [uQ']   [sin(xi),  cos(xi)]   [uq']
%          \_______________/
%                 \/
%     T_xi in transfer matrix form
% where "xi" is the steady-state angle that local steady frame axes lead
% the global steady frame axes, which is also the local signals lag the
% global signals.
%
% Effect on the transfer function: similarity transform
% G_D'Q' = T_xi * G_d'q' * inv(T_xi)

% Get the dimension
[~,lu_xi,ly_xi] = SimplexPS.SsGetDim(Se);

% Transformation matrix
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

% Set the SS system back
Gm = Se;

end

%% Get the descritpor state space model
% Change SS system to DSS system
An = Gm.A; Bn = Gm.B; Cn = Gm.C; Dn = Gm.D; 
En = eye(size(An));
GmDSS = dss(An,Bn,Cn,Dn,En);

%% Get the object model
% Create an object
GmObj = SimplexPS.Class.ModelBase;

% Load the model
GmObj.SetDSS(GmObj,GmDSS);

% Get the strings
GmObj.SetString(GmObj,StateString,InputString,OutputString);

% Switch input and output for required device
if Flag_SwitchInOut == 1
    GmObj = SimplexPS.ObjSwitchInOut(GmObj,2);
end

% Check dimension mismatch
SimplexPS.ObjCheckDim(GmObj);

end