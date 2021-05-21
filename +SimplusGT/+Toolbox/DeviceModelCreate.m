% This function creates the descriptor state space model for apparatuses
% connected to buses.

% Author(s): Yitong Li, Yunjie Gu

%% Notes:
%
% Note that frequency in power flow can be different to frequency in
% parameters. The frequency in power flow is steady-state frequency.


%% References:

% Y. Gu, Y. Li, Y. Zhu, T. C. Green, "Impedance-based whole-system modeling
% for a composite grid via embedding of frame dynamics." IEEE Transactions
% on Power Systems, Early Access.

% Y. Li, Y. Gu, T. C. Green, "Interpreting frame tramsformations in ac
% systems as diagonolization of harmonic transfer function." IEEE
% Transctions on Circuit and Systems, 2020.

%% function
function [GmObj,GmDSS,DevicePara,DeviceEqui,DiscreDampingResistor,OtherInputs,StateStr,InputStr,OutputStr] ...
        = DeviceModelCreate(DeviceBus,Type,PowerFlow,Para,Ts,ListBus) 

%% Create an object
SwInOutFlag = 0;   % Default: do not need to switch input and output
switch floor(Type/10)
    
    % Notes:
    %
    % The parameter order listed below will also influence the parameter
    % order in the system object, but will not influence the input order
    % from the user data.
    
    % =======================================
    % Ac apparatuses
    % =======================================
    % ### Synchronous generator
    case 0      % Type 0-9
        Device = SimplusGT.Class.SynchronousMachine('DeviceType',Type);
        Device.Para = [ Para.J;
                        Para.D;
                        Para.wL;
                        Para.R;
                        Para.w0];       % (5)
                    
    % ### Grid-following inverter
    case 1      % Type 10-19
        if Type~=19
            Device = SimplusGT.Class.GridFollowingVSI('DeviceType',Type);
        else
            Device = SimplusGT.Class.GridFollowingInverterStationary('DeviceType',Type);
        end
        Device.Para = [ Para.C_dc;
                        Para.V_dc;
                        Para.f_v_dc;
                        Para.f_pll;
                        Para.f_tau_pll;
                        Para.f_i_dq;
                        Para.wLf;
                        Para.R;
                        Para.w0];
                   
    % ### Grid-forming inverter
    case 2  % Type 20-29
        Device = SimplusGT.Class.GridFormingVSI('DeviceType',Type);
        Device.Para = [ Para.wLf;
                        Para.Rf;
                        Para.wCf;
                        Para.wLc;
                        Para.Rc;
                        Para.Xov;
                        Para.Dw;
                        Para.fdroop;
                        Para.fvdq;
                        Para.fidq;
                        Para.w0];
                   
    % ### Ac infinite bus
    case 9
        Device = SimplusGT.Class.InfiniteBusAc;
        Device.Para  = [];
        % Because the infinite bus is defined with "i" input and "v" output,
        % they need to be switched finally.
        SwInOutFlag = 1;
        SwInOutLength = 2;
   
    % ### Ac floating bus
    case 10
        Device = SimplusGT.Class.FloatingBusAc;
        Device.Para = [];
        
	% =======================================
    % Dc apparatuses
    % =======================================
    % ### Grid feeding buck converter
    case 101
    	Device = SimplusGT.Class.GridFeedingBuck('DeviceType',Type);
        Device.Para = [ Para.Vdc;
                        Para.Cdc;
                        Para.wL;
                        Para.R;
                        Para.fi;
                        Para.fvdc;
                        Para.w0];
  	% ### Dc infinite bus
    case 109
        Device = SimplusGT.Class.InfiniteBusDc;
        Device.Para  = [];
        SwInOutFlag = 1;
        SwInOutLength = 1;
   
    % ### Dc floating bus
    case 110
        Device = SimplusGT.Class.FloatingBusDc;
        Device.Para = [];
        
   	% =======================================
    % Interlinking apparatuses
    % =======================================
    case 200
        Device = SimplusGT.Class.InterlinkAcDc('DeviceType',Type);
        Device.Para = [ Para.C_dc;
                        Para.wL_ac;
                        Para.R_ac;
                        Para.wL_dc;
                        Para.R_dc;
                        Para.fidq;
                        Para.fvdc;
                        Para.fpll;
                        Para.w0];
    
    % ### Otherwise
    otherwise
        error(['Error: apparatus type']);
end

%% Calculate the linearized state space model
Device.DeviceType = Type;                           % Device type
Device.Ts = Ts;                                     % Samping period
Device.PowerFlow = PowerFlow;                       % Power flow data
Device.SetString(Device);                           % Set strings automatically
Device.SetEquilibrium(Device);                      % Calculate the equilibrium
[x_e,u_e,y_e,xi] = Device.GetEquilibrium(Device);   % Get the equilibrium
Device.SetSSLinearized(Device,x_e,u_e);             % Linearize the model

[~,ModelSS] = Device.GetSS(Device);                % Get the ss model
[StateStr,InputStr,OutputStr] ...
    = Device.GetString(Device);                    % Get the string

% Set ElecPortIOs and OtherInputs
if Type<1000
    Device.ElecPortIOs = [1,2];
    OtherInputs = u_e(3:end,:);     % dq frame ac apparatus
elseif 1000<=Type && Type<2000
    Device.ElecPortIOs = [1];
    OtherInputs = u_e(2:end,:);     % dc apparatus
elseif 2000<=Type && Type<3000
    Device.ElecPortIOs = [1,2,3];
    OtherInputs = u_e(4:end,:);     % ac-dc apparatus
else
    error(['Error']);
end

% Link the IO ports to bus number
InputStr = SimplusGT.AddNum2Str(InputStr,DeviceBus);
OutputStr = SimplusGT.AddNum2Str(OutputStr,DeviceBus);

% For 2-bus apparatus, adjust electrical port strings
if length(DeviceBus)==2  % A multi-bus apparatus 
    InputStr{1} = ['v_d',num2str(DeviceBus(1))];
    InputStr{2} = ['v_q',num2str(DeviceBus(1))];
   	OutputStr{1} = ['i_d',num2str(DeviceBus(1))];
    OutputStr{2} = ['i_q',num2str(DeviceBus(1))];
    
  	InputStr{3} = ['v',num2str(DeviceBus(2))];
    OutputStr{3} = ['i',num2str(DeviceBus(2))];
elseif length(DeviceBus) == 1
else
    error(['Error: Each apparatus can only be connected to one or two buses.']);
end

% Get the swing frame system model
Gm = ModelSS;   

% Output
DevicePara = Device.Para; 
DeviceEqui = {x_e,u_e,y_e,xi};

% Output the discretization damping resistance for simulation use
if Type<90 || (1000<=Type && Type<1090) || (2000<=Type && Type<2090)
    % CalcRv_old();
    
    Device.SetDynamicSS(Device,x_e,u_e);
    DiscreDampingResistor = Device.GetVirtualResistor(Device);
else
    DiscreDampingResistor = -1;
end

%% Check if the apparatus needs to adjust its frame
if (Type>=90 && Type<1000)
    % No need for frame dynamics embedding
elseif (Type>=1000 && Type<2000)
    % No need for frame dynamics embedding
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
% Notes: 
% New system has one more new output "w/s", which is added to the head of
% the original output vector. New system also has a new state "epsilon",
% which is added to the head of the original state vector. This means that
% the w has to be one of the outputs of the original model.
for i = 1:length(OutputStr)
    w_port = SimplusGT.AddNum2Str({'w'},DeviceBus);
    w_port = w_port{1};
    if strcmpi(OutputStr(i),w_port)
        ind_w = i;
        break
    end
end
if isempty(ind_w)
    error(['Error: w has to be in the output of the original apparatus model.']);
end
[~,~,ly1_w] = SimplusGT.SsGetDim(Gm);
Aw = 0;
Bw = zeros(1,ly1_w);   Bw(ind_w) = 1;
Cw = [1;zeros(ly1_w,1)];
Dw = [zeros(1,ly1_w);eye(ly1_w)];
Se = ss(Aw,Bw,Cw,Dw);
Se = series(Gm,Se);
StateStr = [{'epsilon'},StateStr];

% Embed the frame dynamics
Sfb = ss([],[],[],V0);
Se = feedback(Se,Sfb,[1,2],[1]);            % v = v' - V0 * w/s

[~,~,ly2_w] = SimplusGT.SsGetDim(Se);
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
[~,lu_xi,ly_xi] = SimplusGT.SsGetDim(Se);

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
GmObj = SimplusGT.Class.ModelBase;

% Load the model
GmObj.SetDSS(GmObj,GmDSS);

% Get the strings
GmObj.SetString(GmObj,StateStr,InputStr,OutputStr);

% Switch input and output for required apparatus
if SwInOutFlag == 1
    GmObj = SimplusGT.ObjSwitchInOut(GmObj,SwInOutLength);
end

% Check dimension mismatch
SimplusGT.ObjCheckDim(GmObj);
[StateStr,InputStr,OutputStr] = GmObj.GetString(GmObj);

end