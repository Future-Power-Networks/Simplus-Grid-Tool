% This is the main function for analying the power grids by
% communication theory

%% Enable settings
% Enable inner loop
Enable_VoltageNode_InnerLoop    = 1;    % 1/0: star-delta conversion for flux inductance of voltage node 
Enable_CurrentNode_InnerLoop    = 1;    % 1/0: inner-current loop impedance of current node

% PLL settings
Enable_vq_PLL                   = 1;    % 1/0: change Q-PLL to vq-PLL
Enable_PLL_LPF                  = 1;    % 1/0: if the PLL with an additional 100Hz LPF
w_PLL_LPF                       = 2*pi*100;     % Bandwidth of the PLL LPF

% Enable plot
Enable_Plot_Eigenvalue          = 1;    % 1/0: Plot eigenvalues.

%% Fault settings
NfaultBus = 37;
Tfault = 1/60*3;

%% Update power flow
[V,I] = SimplusGT.Synchron.UpdateVI(PowerFlowNew);

%% Nodal admittance matrix
% Symbolic
s = sym('s');

% Calculate nodal admittance matrix
fprintf('Calculate nodal admittance matrix...\n')
Ybus = SimplusGT.Synchron.YbusCalcSym(ListLineNew,Wbase,'albe');       % We get the alpha/beta frame nodal admittance matrix

% Convert Ybus to double
dW = 1e-10*(1+Wbase);
Ybus_ = double(subs(Ybus,'s',1i*(Wbase+dW)));   % Used for calculating derivative numerically
Ybus = double(subs(Ybus,'s',1i*Wbase));
% Nnotes:
% Ybus should statisfy: I = Ybus*V

% Get Ybus during fault
ListLineFault = SimplusGT.Synchron.GetListLineFault(NfaultBus,ListLineNew,ListBus);
YbusFault = SimplusGT.Synchron.YbusCalcSym(ListLineFault,Wbase,'albe');
YbusFault = double(subs(YbusFault,'s',1i*Wbase));

%% Reorder the Data
fprintf('Reorder the data...\n')
% Notes:
%
% In this code, the bus/node should be orderred in this sequence:
% [all voltage nodes, all current nodes, all floating bus nodes], i.e.,
% [v bus, ..., v bus, i bus, ..., i bus, f bus, ... f bus].
% Hence, in this section, we re-order the data obtained from excel first,
% to make sure that this required sequence can be obtained. Noting that,
% the device data, the power flow data, and the network line data should
% all be re-orderred.

% [ApparatusSourceType,ApparatusParaNew,IndexEbus,IndexVbus,IndexIbus,IndexFbus,OrderOld2New,Ybus,Ybus_,ExistVbus,ExistIbus,ExistFbus,V,I] = SimplusGT.Synchron.FunReorderData(ListBus,ApparatusType,V,I,Ybus,Ybus_,Para);
SimplusGT.Synchron.ReorderData();

%% The Influence of node type and their parameters on Ybus
fprintf('Consider the influence of node type on node admittance matrix...\n')

% Notes:
% If using s-domain Ybus calculation and using vectors as input, then, the
% system admittance matrix has to be symmetric. Fortunately, the passive
% component, and the inner loops of inverters, are indeed symmtric in
% complex dq and alpha/beta frame.

% Find the node index
SimplusGT.Synchron.FindNodeIndex();

% Handle voltage, current, and floating nodes
% [V,I,J,D,kp_pll,ki_pll,Ybus,Ybus_,N_Node] = SimplusGT.Synchron.FunHandleNode(N_Bus,V,I,ApparatusParaNew,Ybus,Ybus_,NumVbus1st,NumIbus1st,NumFbus1st,ExistVbus,ExistIbus,ExistFbus,Wbase,dW,Enable_VoltageNode_InnerLoop,Enable_CurrentNode_InnerLoop);
SimplusGT.Synchron.HandleNode();

%% Network matrix
fprintf('Calculate network matrix: hybrid admittance/impedance matrix, or equivalently channel gain...\n')

% Convert the nodol admittance matrix to hybrid admittance/impedance matrix
Gbus = SimplusGT.Synchron.HybridMatrixYZ(Ybus,NumIbus1st);
% GbusVI  = Gbus;
% GbusVIF = SimplusGT.Synchron.HybridMatrixYZ(YbusVIF,NumIbus1st);

% For numerically calculating GbusPrime later
Gbus_ = SimplusGT.Synchron.HybridMatrixYZ(Ybus_,NumIbus1st);

% Notes:
% It should be ensured that the buses are listed in the form like this:
% [Vbus1, Vbus2, Vbus3, Ibus4, Ibus5, ...]
    
Gbus = -Gbus;  	% Change the power direction to load convention.
              	% Noting that this operation is different from Ybus
                % = -Ybus if the system has current nodes. The current
                % direction is not important actually.
                
% For numerically calculating GbusPrime
Gbus_ = -Gbus_;

% Get G_prime
% Notes: It is calculaed by numerical method
GbusPrime = (Gbus_ - Gbus)/(1i*dW);         	% Consider

% Get fault Gbus
GbusFault = SimplusGT.Synchron.HybridMatrixYZ(YbusFault,NumIbus1st);
GbusFault = - GbusFault;

%% 
fprintf('Calculate network matrix: complex power...\n')
% Update input and output so that they correspond to the hybrid
% admittance/impedance matrix, i.e., Output = -Gbus*Input
Input = [V(1:NumIbus1st-1);
         I(NumIbus1st:end)];
Output = [I(1:NumIbus1st-1);
          V(NumIbus1st:end)];
      
% Normalize the current node because of PLL
InputNormalized = Input;        % Initialize
for i = 1:length(Input)
    if ApparatusSourceType(i) == 2
        % Notes:
        % If current source, then normalize the current of it, in order to
        % match the actual feedback signal, i.e., voltage rather than
        % power, of the PLL inverter later.
        %
        % This vq-PLL effect is considerred into the S matrix next, rather
        % than the T or H matrix. This effect only needs to be considerred
        % once.
        if Input(i) == 0
            InputNormalized(i,1) = 0;
        else
            InputNormalized(i,1) = Input(i)/abs(Input(i));
        end
    end
end

% Get S matrix
if Enable_vq_PLL
    S = conj(InputNormalized)*transpose(Input);
else
    S = conj(Input)*transpose(Input);
end

%% 
fprintf('Calculate network matrix: mu, GAMMA, and gamma...\n')
% Get mu
for i = 1:N_Node
    if ApparatusSourceType(i) == 1          % Voltage node
        mu(i) = 0;         % W = P
    elseif ApparatusSourceType(i) == 2      % Current node
        mu(i) = pi/2;      % W = Q
        if Enable_vq_PLL
            theta_i = angle(-I(i));
            theta_v = angle(V(i));
            mu(i) = pi/2 - (theta_i-theta_v);      % The Q direction is changed to vq direction.
        end
    else
        error(['Error']);
    end
end

% Get GAMMA and gamma
for m = 1:N_Node
    for n = 1:N_Node
        GAMMA(m,n) = abs(Gbus(m,n)*S(m,n));
        gamma(m,n) = pi/2 + mu(m) - angle(Gbus(m,n));
    end
end

% Get GAMMA and gamma during fault
for m = 1:N_Node
    for n = 1:N_Node
        GAMMAFault(m,n) = abs(GbusFault(m,n)*S(m,n));
        gammaFault(m,n) = pi/2 + mu(m) - angle(GbusFault(m,n));
    end
end

%%
fprintf('Calculate network matrix: inertia, damping...\n')
% Initialize
Hmat = eye(N_Node);
Hinv = inv(Hmat);
Dmat = eye(N_Node);

% Update voltage node
if ExistVbus == 1
for i = 1:(NumIbus1st-1)
    % The inertia of a SG is J
    Hmat(i,i) = J{i};
    Dmat(i,i) = D{i};
    
  	Hinv(i,i) = 1/Hmat(i,i);
    Hinv(i,i) = double(Hinv(i,i));
end
end

% Update current node
if ExistIbus == 1
for i = NumIbus1st:N_Node
    % The inertia of an inverter is ki_pll.
    if Enable_PLL_LPF == 0
        Hmat(i,i) = 1/ki_pll{i};                        % PI format
        Dmat(i,i) = kp_pll{i}/ki_pll{i};
    else                                    
        Hmat(i,i) = 1/(w_PLL_LPF*kp_pll{i});         	% LPF format
        Dmat(i,i) = 1/kp_pll{i};
        ki_pll{i} = 0;
    end
    Hinv(i,i) = 1/Hmat(i,i);
    Hinv(i,i) = double(Hinv(i,i));
end
end

% Power reference
for i = 1:N_Node
    if ApparatusSourceType(i) == 1          % Voltage node
        Wref(i) = -PowerFlowNew{i}(1);
    elseif ApparatusSourceType(i) == 2      % Current node
        Wref(i) = 0;
 	end
end

%%
if 0
SimplusGT.Synchron.SmallSignalAnalysis();
end

%%
if 1
fprintf('Run time-domain simulation...\n')
SimplusGT.Synchron.Simulation(Input,GAMMA,gamma,GAMMAFault,gammaFault,Hmat,Dmat,Wref,Wbase,N_Node,FigN);
end

%% Clear useless variables

clear('V','I');
clear('Ybus','Ybus_','YbusFault','Gbus','Gbus_','GbusFault','GbusPrime');
clear('J','D','kp_pll','ki_pll');
clear('ExistVbus','ExistIbus','ExistFbus');
clear('IndexEbus','IndexVbus','IndexIbus','IndexFbus');
clear('feedin','feedoutL1','feedoutL2');
clear('GAMMA','gamma','GAMMAFault','gammaFault','mu');
clear('K','ang_K','FreqShift','ang_FreqShift');
clear('ApparatusParaNew','ApparatusSourceType');
clear('N_Node');