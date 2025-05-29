% This function creates the descriptor state space model for apparatuses
% connected to buses.

% Author(s): Yitong Li, Yunjie Gu

%% Notes:
%
% Note that frequency in power flow can be different to frequency in
% parameters. The frequency in power flow is steady-state frequency.


%% References:

% Yitong Li, Yunjie Gu, Timothy C. Green, "Mapping of dynamics between
% mechanical and electrical ports in SG-IBR composite grids" arXiv,
% preprint.

% Yunjie Gu, Yitong Li, Yue Zhu, Timothy C. Green, "Impedance-based
% whole-system modeling for a composite grid via embedding of frame
% dynamics." IEEE Transactions on Power Systems, 2021.

% Yitong Li, Yunjie Gu, Timothy C. Green, "Interpreting frame
% tramsformations in ac systems as diagonolization of harmonic transfer
% function." IEEE Transctions on Circuit and Systems, 2020.

%% function
function [GmObj,GmDSS,ApparatusPara,ApparatusEqui,DiscreDampingResistor,OtherInputs,StateStr,InputStr,OutputStr] ...
        = ApparatusModelCreate(ApparatusBus,Type,PowerFlow,Para,Ts,ListBus) 

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
        Apparatus = SimplusGT.Class.SynchronousMachine('ApparatusType',Type);
        Apparatus.Para = [ Para.J;
                        Para.D;
                        Para.wL;
                        Para.R;
                        Para.w0];       % (5)
                    
    % ### Grid-following inverter
    case 1      % Type 10-19
        if Type~=19
            Apparatus = SimplusGT.Class.GridFollowingVSI('ApparatusType',Type);
        else
            Apparatus = SimplusGT.Class.GridFollowingInverterStationary('ApparatusType',Type);
        end
        Apparatus.Para = [ Para.C_dc;
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
        Apparatus = SimplusGT.Class.GridFormingVSI('ApparatusType',Type);
        Apparatus.Para = [ Para.wLf;
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
        Apparatus = SimplusGT.Class.InfiniteBusAc;
        Apparatus.Para  = [];
        % Because the infinite bus is defined with "i" input and "v" output,
        % they need to be switched finally.
        SwInOutFlag = 1;
        SwInOutLength = 2;
   
    % ### Ac floating bus
    case 10
        Apparatus = SimplusGT.Class.FloatingBusAc;
        Apparatus.Para = [];
        
	% =======================================
    % Dc apparatuses
    % =======================================
    % ### Grid feeding buck converter
    case 101
    	Apparatus = SimplusGT.Class.GridFeedingBuck('ApparatusType',Type);
        Apparatus.Para = [ Para.Vdc;
                        Para.Cdc;
                        Para.wL;
                        Para.R;
                        Para.fi;
                        Para.fvdc;
                        Para.w0];
  	% ### Dc infinite bus
    case 109
        Apparatus = SimplusGT.Class.InfiniteBusDc;
        Apparatus.Para  = [];
        SwInOutFlag = 1;
        SwInOutLength = 1;
   
    % ### Dc floating bus
    case 110
        Apparatus = SimplusGT.Class.FloatingBusDc;
        Apparatus.Para = [];
        
   	% =======================================
    % Interlinking apparatuses
    % =======================================
    case 200
        Apparatus = SimplusGT.Class.InterlinkAcDc('ApparatusType',Type);
        Apparatus.Para = [ Para.C_dc;
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
Apparatus.ApparatusType = Type;                           % Apparatus type
Apparatus.Ts = Ts;                                     % Samping period
Apparatus.PowerFlow = PowerFlow;                       % Power flow data
Apparatus.SetString(Apparatus);                           % Set strings automatically
Apparatus.SetEquilibrium(Apparatus);                      % Calculate the equilibrium
[x_e,u_e,y_e,xi] = Apparatus.GetEquilibrium(Apparatus);   % Get the equilibrium
Apparatus.SetSSLinearized(Apparatus,x_e,u_e);             % Linearize the model

[~,ModelSS] = Apparatus.GetSS(Apparatus);                % Get the ss model
[StateStr,InputStr,OutputStr] ...
    = Apparatus.GetString(Apparatus);                    % Get the string

% Set ElecPortIOs and OtherInputs
if Type<1000
    Apparatus.ElecPortIOs = [1,2];
    OtherInputs = u_e(3:end,:);     % dq frame ac apparatus
elseif 1000<=Type && Type<2000
    Apparatus.ElecPortIOs = [1];
    OtherInputs = u_e(2:end,:);     % dc apparatus
elseif 2000<=Type && Type<3000
    Apparatus.ElecPortIOs = [1,2,3];
    OtherInputs = u_e(4:end,:);     % ac-dc apparatus
else
    error(['Error']);
end

% Add bus number index to all IO ports
InputStr = SimplusGT.AddNum2Str(InputStr,ApparatusBus);
OutputStr = SimplusGT.AddNum2Str(OutputStr,ApparatusBus);
if ~isempty(StateStr)
StateStr = SimplusGT.AddNum2Str(StateStr,ApparatusBus);
end
% Notes: 
% For interlink, the port number would be like vd2-3, w2-3, etc. The
% connection between Gm and Zbus directly dependents on the IO bus number
% index, so for interlink, the electrical port number should be adjusted
% next.
if length(ApparatusBus) == 1    % A normal apparauts
elseif length(ApparatusBus)==2  % A multi-bus apparatus, e.g., interlink converter
    % Ac port
    InputStr{1} = ['v_d',num2str(ApparatusBus(1))];
    InputStr{2} = ['v_q',num2str(ApparatusBus(1))];
   	OutputStr{1} = ['i_d',num2str(ApparatusBus(1))];
    OutputStr{2} = ['i_q',num2str(ApparatusBus(1))];
    
    % Dc port
  	InputStr{3} = ['v',num2str(ApparatusBus(2))];
    OutputStr{3} = ['i',num2str(ApparatusBus(2))];
else
    error(['Error: Each apparatus can only be connected to one or two buses.']);
end

% Get the swing frame system model
Gm = ModelSS;   

% Output
ApparatusPara = Apparatus.Para; 
ApparatusEqui = {x_e,u_e,y_e,xi};

% Output the discretization damping resistance for simulation use
if Type<90 || (1000<=Type && Type<1090) || (2000<=Type && Type<2090)
    % CalcRv_old();
    
    Apparatus.SetDynamicSS(Apparatus,x_e,u_e);
    DiscreDampingResistor = Apparatus.GetVirtualResistor(Apparatus);
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
    w_port = SimplusGT.AddNum2Str({'w'},ApparatusBus);
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
StrEpsilon = SimplusGT.AddNum2Str({'epsilon'},ApparatusBus);
StateStr = [StrEpsilon,StateStr];

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

%% Get the descriptor state space model
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